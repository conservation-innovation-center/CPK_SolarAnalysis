#' Standardize continuous data to z-score
#' 
#' @param x numeric vector
#' @return vector of same length as \code{x} containing z-transformed values
#' @examples 
#' v <- runif(100, 10, 12)
#' calc_z(v)
calc_z <- function(x){
  mn <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  z <- (x-mn)/s
  return(z)
}

#' Standardize continuous data to (0,1)
#' 
#' @param x numeric vector
#' @return vector of same length as \code{x} containing standardized values
#' @examples 
#' v <- runif(100, 10, 12)
#' standardize(v)
standardize <- function(x){
  mx <-max(x, na.rm = TRUE)
  mn <-min(x, na.rm = TRUE)-0.0001
  z <- ((x-mn)/(mx-mn))
  return(z)
}

logistic <- function(x){
  y <- exp(x)/(exp(x) + 1)
  return(y)
}

logit <- function(x){
  y <- log(x/(1-x))
  return(y)
}

plot_post_hist <- function(jags_out, param){
  # get index of desired parameter in jags simulation object
  index = grep(param, names(jags_out$sims.list))
  hist(out$sims.list[[index]], breaks = 50, col = 'grey', freq = F, xlab = param)
  abline(v = quantile(out$sims.list[[index]], prob=c(0.025, 0.975)), col = 'red', lwd = 2)
}

#' Create a dataset that can be used by rjags and jagsUI
#' 
#' @param dataframe a data frame containing all variables (must be a dataframe NOT A TIBBLE)
#' @param random a chr vector identifying the column(s) in dataframe containing random variables
#' @param continuous a chr vector identifying the column(s) in dataframe containing continuous predictor variables
#' @param categorical a chr vector identifying the column(s) in dataframe containing categorical predictor variables
#' @param responses a chr vector identifying the column(s) in dataframe containing response variables
#' @return data list that can be provided to the 'data' argument of jagsUI::jags()
#' @example
#' dat <- make_jags_data(
#' dataframe = cleanDF,
#' random = c('statef'),
#' continuous = c('impervious16', 'open16', 'tree_cover16', 'slope'),
#' categorical = c('tmax', 'statef'),
#' response = c('solar', 'ttd'))
make_jags_data <- function(dataframe, random = NULL, continuous, categorical, responses){
  # first create a list of our continuous variables
  cont <- list()
  for(var in continuous){
    cont <- append(cont, list(dataframe[,c(var)]))
  }
  names(cont) <- continuous
  
  # then create a list of our categorical variables
  cate <- list()
  for(var in categorical){
    cate <- append(cate, list(dataframe[,c(var)]))
  }
  names(cate) <- categorical
  
  # list of our response variables
  resp <- list()
  for(var in responses){
    resp <- append(resp, list(dataframe[,c(var)]))
  }
  names(resp) <- responses
  
  # list of our constants
  const <- list(
    N = nrow(dataframe)
  )
  
  if(length(random)>0){
    ## TO DO: Make this dynamic. Currently assumes we will only ever use state as a random variable
    for(var in random){
      const <- append(const, list(S = length(unique(dataframe[,c(var)]))))
    }
  }
  
  # constrict all continuous variables to [0-1]
  # then combine with categorical variables
  dat <- map(cont, standardize)%>%
    append(cate)%>%
    append(const)%>%
    append(resp)
  
  return(dat) 
}

#' write a binomial jags model 
#' 
#' Given the names of response and predictor variables, create text representing a simple
#' binomial model in jags model format and write to a specified file
#' 
#' @param filename string specifying destination file
#' @param response string identifying name of response variable in jags data list
#' @param predictors character vector containing names of predictor varibles in jags data list
#' @returns function to generate initial values
write_jags_binomial <- function(filename, response, predictors){
  file.create(filename)
  responseString <- sprintf('%s[i] ~ dbern(psi[i])', response)
  # create a string of predictor variables following beta_var*var[i]+...
  predictorString <- paste('logit(alpha)', paste('beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  rep('dnorm(0, 0.0001', length(predictors))
  priorsString <- paste('alpha ~ dnorm(0, 0.0001)', paste('beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    '), sep = '\n    ')
  inits <- list(rnorm)
  # for(var in predictors) {
  #   prior <- sprintf('beta_%s ~ dnorm(0, 0.0001)', var)
  #   priorsString <- paste(priorsString, prior, sep = '\n    ')
  #   # add predictor to our list of inits
  # 
  # }
  baseString <-"
  model{
    #Likelihood
    # binomial data model
    for (i in 1:N){
      %s
      logit(psi[i]) <- %s
    }
    
    #Priors
    %s
  }"
  finalString = sprintf(baseString, responseString, predictorString, priorsString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ls <- list(alpha = rnorm(1,0,1))
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[2:length(ls)] <- paste('beta', predictors, sep = "_")
    return(ls)
  }
  
  return(inits)
}

#' write a binomial jags model with random effects
#' 
#' Given the names of response and predictor variables, create text representing a simple
#' binomial model in jags model format and write to a specified file
#' 
#' @param filename string specifying destination file
#' @param response string identifying name of response variable in jags data list
#' @param predictors character vector containing names of predictor varibles in jags data list
#' @param random string indentifying the name of variable used for random effects on intercept
#' @returns function to generate initial values
#' 
#' @details generates a text file specifying a jags model in the following format: \code{
#' model{
#' #Likelihood
#' alpha_shp2 <- alpha_shp1*((1/mu) - 1)
#' alpha_shp1 <- (((1-mu)/v)-1/mu )* pow(mu, 2)
#' # for each level of random
#' for (s in 1:S){
#'  # alpha[s] ~ dnorm(alpha_mu, alpha_tau)
#'  alpha[s] ~ dbeta(alpha_shp1, alpha_shp2)
#'  
#' }
#' # binomial data model
#' for (i in 1:N){
#'  solar[i] ~ dbern(psi[i])
#'  logit(psi[i]) <- logit(alpha[s[i]]) + \strong{beta_predictors*predictors[i]}
#' }
#'
#' #Priors
#' # alpha_mu ~ dnorm(0, 0.0001)
#' # alpha_tau <- pow(sigma, -2)
#' # sigma ~ dunif(0, 1000)
#' mu ~ dunif(0, 1)
#' v ~ dunif(0, 0.2)

#' beta_predictors[i]~dnorm(0, 0.0001)
#' ...
#' }
#'}
write_jags_binomial_raneff <- function(filename, response, predictors, random){
  file.create(filename)
  responseString <- sprintf('%s[i] ~ dbern(psi[i])', response)
  alphaString <- sprintf('logit(alpha[%s[i]])', random)
  # create a string of predictor variables following beta_var*var[i]+...
  predictorString <- paste(alphaString, paste('beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  priorsString <- paste('beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')

  # for(var in predictors) {
  #   prior <- sprintf('beta_%s ~ dnorm(0, 0.0001)', var)
  #   priorsString <- paste(priorsString, prior, sep = '\n    ')
  #   # add predictor to our list of inits
  # 
  # }
  baseString <-"
  model{
    #Likelihood
    alpha_shp2 <- alpha_shp1*((1/mu) - 1)
    alpha_shp1 <- (((1-mu)/v)-1/mu )* pow(mu, 2)
    # for each level of random
    for (s in 1:S){
      # alpha[s] ~ dnorm(alpha_mu, alpha_tau)
      alpha[s] ~ dbeta(alpha_shp1, alpha_shp2)
  
    }
    # binomial data model
    for (i in 1:N){
      %s
      logit(psi[i]) <- %s
    }
    
    #Priors
    # alpha_mu ~ dnorm(0, 0.0001)
    # alpha_tau <- pow(sigma, -2)
    # sigma ~ dunif(0, 1000)
    mu ~ dunif(0, 1)
    v ~ dunif(0, 0.2)
    
    %s
  }"
  finalString = sprintf(baseString, responseString, predictorString, priorsString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ls <- list(
      mu = 0.5,
      v = 0.125
      # alpha_mu = rnorm(1,0,1),
      # sigma = runif(0, 1000)
      )
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[3:length(ls)] <- paste('beta', predictors, sep = "_")
    return(ls)
  }
  
  return(inits)
}

#' write a joint binomial weibull jags model with random effects
#' 
#' @description 
#' Given the names of binary and continuous response variables, and predictor variables
#'  create text representing a simple binomial model in jags model format and write to a specified file
#' 
#' @param filename string specifying destination file
#' @param response string identifying name of response variable in jags data list
#' @param predictors character vector containing names of predictor varibles in jags data list
#' @param random string indentifying the name of variable used for random effects on intercept
#' @returns function to generate initial values
#' 
#' @details 
#' As a side effect, the funciton generates a text file specifying a jags model using the template: \code{
#' model{
#' #Likelihood
#'
#' # hyperparameters for random effects on 
#' shape_shape <- pow(shape_mu, 2)/shape_v
#' shape_rate <- shape_mu/shape_v
#'
#' rate_alpha_shape <- pow(rate_alpha_mu, 2)/rate_alpha_v
#' rate_alpha_rate <- rate_alpha_mu/rate_alpha_v

#' psi_alpha_shp2 <- psi_alpha_shp1*((1/psi_alpha_mu) - 1)
#' psi_alpha_shp1 <- (((1-psi_alpha_mu)/psi_alpha_v)-1/psi_alpha_mu )* pow(psi_alpha_mu, 2)
#'
#' # random intercept per state on binomial intercept and weibul shape
#' for (s in 1:S){
#'   psi_alpha[s] ~ dbeta(psi_alpha_shp1, psi_alpha_shp2)
#'  # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
#'  # BUGS uses shape and rate parameterization of gamma
#'   shape[s] ~ dgamma(shape_shape, shape_rate)
#'  # the rate or scale parameter is [0, Inf] variability in ttd data
#'   rate_alpha[s] ~ dgamma(rate_alpha_shape, rate_alpha_rate)
#' }
#'
#' # binomial data model
#' for (i in 1:N){
#'  # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
#'   \code{bin_response}[i] ~ dbern(psi[i])
#'   # probability of ever being developed linear fxn of covariates with state-specific intercept
#'   logit(psi[i]) <- logit(psi_alpha[\code{random}[i]]) + \strong{psi_beta_predictors*predictors[i]}
#'  
#'  # time to detection is a weibull process with state-specific shape and rate determined by covariates
#'   \code{cont_response}[i] ~ dweib(shape[\code{random}[i]], rate[i])
#'   log(rate[i]) <- log(rate_alpha[\code{random}[i]]) + \strong{rate_beta_predictors*predictors[i]}
#'}
#'
#' #Priors
#' # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
#' # to simulate gamma(0.0001, 0.0001) with no state effect
#' shape_mu ~ dunif(0, 5)
#' shape_v ~ dunif(0, 1000) 
#' rate_alpha_mu ~ dunif(0, 5)
#' rate_alpha_v ~ dunif(0, 1000)
#' psi_alpha_mu ~ dunif(0, 1) 
#' psi_alpha_v ~ dunif(0, 0.2)
#' psi_beta_predictors[i] ~ dnorm(0, 0.0001)
#' ...
#' psi_beta_predictors[n] ~ dnorm(0, 0.0001)
#' rate_beta_predictors[i] ~ dnorm(0, 0.0001)
#' ...
#' rate_beta_predictors[n] ~ dnorm(0, 0.0001)
#' }
write_jags_weibull_raneff <- function(filename, bin_response, cont_response, predictors, random){
  file.create(filename)
  binResponseString <- sprintf('%s[i] ~ dbern(psi[i])', bin_response)
  binAlphaString <- sprintf('logit(psi_alpha[%s[i]])', random)
  # create a string of predictor variables following beta_var*var[i]+...
  binPredictorString <- paste(binAlphaString, paste('psi_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  contResponseString <- sprintf('%s[i] ~ dweib(shape[%s[i]], rate[i])', cont_response, random)
  contAlphaString <- sprintf('log(rate_alpha[%s[i]])', random)
  contPredictorString <- paste(contAlphaString, paste('rate_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')

  binPriorsString <- paste('psi_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  contPriorString <- paste('rate_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  
  baseString <- "
  model{
    #Likelihood
    
    # hyperparameters for random effects on 
    shape_shape <- pow(shape_mu, 2)/shape_v
    shape_rate <- shape_mu/shape_v
    
    rate_alpha_shape <- pow(rate_alpha_mu, 2)/rate_alpha_v
    rate_alpha_rate <- rate_alpha_mu/rate_alpha_v
    
    # psi_alpha_shp2 <- psi_alpha_shp1*((1/psi_alpha_mu) - 1)
    # psi_alpha_shp1 <- (((1-psi_alpha_mu)/psi_alpha_v)-1/psi_alpha_mu )* pow(psi_alpha_mu, 2)
    
    # random intercept per state on binomial intercept and weibul shape
    for (s in 1:S){
      # psi_alpha[s] ~ dbeta(psi_alpha_shp1, psi_alpha_shp2)
      # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
      # BUGS uses shape and rate parameterization of gamma
      shape[s] ~ dgamma(shape_shape, shape_rate)
      # the rate or scale parameter is [0, Inf] variability in ttd data
      rate_alpha[s] ~ dgamma(rate_alpha_shape, rate_alpha_rate)
    }
    
    # binomial data model
    for (i in 1:N){
      # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
      # %s
      # probability of ever being developed linear fxn of covariates with state-specific intercept
      # logit(psi[i]) <- %s
      
      # time to detection is a weibull process with state-specific shape and rate determined by covariates
      %s
      log(rate[i]) <- %s
      
      # Expected data under current model
      ttd_sim[i] ~ dweib(shape[statei[i]], rate[i])
      # solar_exp[i] ~ dbern(psi[i])
    }
    
    #Priors
    # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
    # to simulate gamma(0.0001, 0.0001) with no state effect
    shape_mu ~ dunif(0, 5)
    shape_v ~ dunif(0, 1000) 
    rate_alpha_mu ~ dunif(0, 5)
    rate_alpha_v ~ dunif(0, 1000)
    # psi_alpha_mu ~ dunif(0, 1) 
    # psi_alpha_v ~ dunif(0, 0.2)
    %s
  }"
  
  finalString = sprintf(
    baseString,
    binResponseString,
    binPredictorString,
    contResponseString,
    contPredictorString,
    #binPriorsString,
    contPriorString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ls <- list(
      # we will chose mean and variance for gamma that corresponds to a 'flat' gamma(0.001, 0.001)
      shape_mu= 1,
      shape_v= 1000,
      rate_alpha_mu = 1,
      rate_alpha_v = 1000
      # we will use mean and variance for a 'flat' beta dist beta(0.5, 0.5)
      # psi_alpha_mu = 0.5, 
      # psi_alpha_v = 0.125
    )
    for(var in rep(predictors)){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[5:(4+length(predictors))] <- paste('rate_beta', predictors, sep = "_")
    # names(ls)[(5+length(predictors)):length(ls)] <- paste('rate_beta', predictors, sep = "_")
    
    return(ls)
  }
  
  return(inits)
  
}

#' calculate the expected response of a binomial model
#' 
#' @description 
#' uses the mean posterior estimates for all estimated parameters in a binomial jags
#' model to estimate the expected response given the required set of predictor variables
#' 
#' @param dataset a jagsUI structured dataset containing predictor variables
#' @param out a jagsUI model output object produced using \code{\link{write_jags_binomial}} or \code{\ling{write_jags_binomial_raneff}}
#' @param predictors character vector containing names of predictor variables. must be supplied in same order as in model creation
#' @param random string containing name of random variable
#' @return vector containing expected responses at each observation based on predictor variable values
#' 
#' @example 
#' continuous <- c('impervious16', 'open16', 'tree_cover16', 'slope') 
#' 
#' dat <- make_jags_data(
#'   dataframe = cleanDF,
#'   random = c('statei'),
#'   continuous = continuous,
#'   categorical = c('tmax', 'statei'),
#'   response = c('solar', 'ttd', 'd')
#' )
#' 
#' inits <- write_jags_binomial_raneff(modfile, 'solar', continuous, 'statei')
#' 
#' out <- jags(data = dat, inits = inits, ...)
#' 
#' mu <- compute_mean_response(dat, out, continuous, 'statei')
#' 
compute_binomial_response <- function(dataset, out, predictors, random = NULL){
  df <- as.data.frame(dataset[c(predictors, random)])
  means <- as.data.frame(out$mean[paste('beta', predictors, sep = "_")])
  pred_func <- function(x){
    if(random){
      int <- logit(out$mean$alpha[df[x, random]])
    }
    else{
      int <- logit(out$mean$alpha)
    }
    y <- logistic(int + sum(df[x, predictors]*means))
    return(y)
  }
  mu <- map_dbl(1:nrow(df), pred_func)
  return(mu)
}

#' Expected response of Weibull model
#'
#' @param dataset a jagsUI structured dataset containing n predictor variables
#' @param params either the sims.list or mean objects from jagsUI out created by \code{\link{write_jags_binomial}} or \code{\ling{write_jags_binomial_raneff}}
#' @param predictors character vector containing names of predictor variables. must be supplied in same order as in model creation
#'
#' @return matrix of expected values of same dimensions as param
#'
#' @examples
#' continuous <- c('impervious16', 'open16', 'tree_cover16', 'slope') 
#' 
#' dat <- make_jags_data(
#'   dataframe = cleanDF,
#'   random = c('statei'),
#'   continuous = continuous,
#'   categorical = c('tmax', 'statei'),
#'   response = c('solar', 'ttd', 'd')
#' )
#' 
#' inits <- write_jags_weibull_raneff(modfile, 'solar', continuous, 'statei')
#' 
#' out <- jags(data = dat, inits = inits, ...)
#' # compute expecte value for each observation using mean parameter posteriors
#' ttd_exp <- compute_weibull_response(dat, out$mean, continuous)
#' 
#' # compute expected value for each observation at each MCMC iteration
#' ttd_exp <- compute_weibull_response(dat, out$sims.list, continuous)
compute_weibull_response <- function(dataset, params, predictors){
  pred_df <- as.data.frame(dataset[predictors])
  betas <- as.data.frame(params[grepl('rate_beta', names(params))])
  # simRates is a niter x nobs matrix. for mean values, niter is 1
  simRates <- t(as.matrix(pred_df) %*% t(as.matrix(betas)))
  if(dim(params$shape)[1] == dim(simRates)[1]){
    fx <- function(x){
      i <- dat$statei[x]
      rate <- exp(simRates[,x] + log(params$rate_alpha[,i]))
      shape <- params$shape[,i]
      mu <- meanWeibullJags(shape, rate)
      return(mu)
    }
  }
  else{
    fx <- function(x){
      i <- dat$statei[x]
      rate <- exp(simRates[,x] + log(params$rate_alpha[i]))
      shape <- params$shape[i]
      mu <- meanWeibullJags(shape, rate)
      return(mu)
    }
  }
  expected <- vapply(1:nrow(pred_df), fx, FUN.VALUE = numeric(nrow(simRates)))
  return(expected)
}


#' JAGS Weibull Distribution
#' 
#' @description 
#' Density function for the Weibull distribution with parameters \code{shape} and \code{rate}
#' as parameterized by JAGS.  
#' \eqn{f(x|v,L) = vLx^(v-1)exp(-Lx^v)}
#'
#' @param x vector of values
#' @param shape shape parameter, v
#' @param rate rate parameter, Lambda
#'
#' @return density of the specified Weibull distribution at x
dweibullJags <- function(x, shape, rate){
  y <- shape*rate*(x^(shape-1))*exp(-rate*(x^shape))
  return(y)
}

#' JAGS Weibull Distribution
#' 
#' @description 
#' Quantile function for the Weibull distribution with parameters \code{shape} and \code{rate}
#' as parameterized by JAGS.  
#' \eqn{F(x|v,L) = 1-exp(-Lx^v)}
#'
#' @param x vector of quantiles
#' @param shape shape parameter, v
#' @param rate rate parameter, Lambda
#'
#' @return quantile of the specified Weibull distribution at x
qweibullJags <- function(x, shape, rate){
  y <- 1-(exp(-rate*(x^shape)))
  return(y)
}

hweibullJags <- function(x, shape, rate){
  y  <- rate*shape*(x^(shape-1))
  return(y)
}

meanWeibullJags <- function(shape, rate){
  y <- (rate^(-1/shape))*gamma(1+(1/shape))
  return(y)
}

write_jags_nbin_raneff <- function(filename, bin_response, cont_response, predictors, random){
  file.create(filename)
  binResponseString <- sprintf('%s[i] ~ dbern(psi)', bin_response)
  # binAlphaString <- sprintf('logit(psi_alpha[%s[i]])', random)
  # create a string of predictor variables following beta_var*var[i]+...
  # binPredictorString <- paste(binAlphaString, paste('psi_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  contResponseString <- sprintf('%s[i] ~ dnegbin(p[%s[i]], 1)', cont_response, random)
  contAlphaString <- sprintf('logit(p_alpha[%s[i]])', random)
  contPredictorString <- paste(contAlphaString, paste('p_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  # binPriorsString <- paste('psi_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  contPriorString <- paste('p_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  
  baseString <- "
  model{
    #Likelihood
    
    # Distribution Parameters for Random Efects on
    # 1. Negative Binomial probability p ~ dbeta(shape1, shape2)
    p_alpha_shp2 <- p_alpha_shp1*((1/p_alpha_mu) - 1)
    p_alpha_shp1 <- (((1-p_alpha_mu)/p_alpha_v)-1/p_alpha_mu )* pow(p_alpha_mu, 2)
    
    # 2. Bernouli probability p ~ dbeta(shape1, shape2)
    # psi_alpha_shp2 <- psi_alpha_shp1*((1/psi_alpha_mu) - 1)
    # psi_alpha_shp1 <- (((1-psi_alpha_mu)/psi_alpha_v)-1/psi_alpha_mu )* pow(psi_alpha_mu, 2)
    
    # random intercept per state on binomial intercept and neg binomial intercept
    for (s in 1:S){
      # psi_alpha[s] ~ dbeta(psi_alpha_shp1, psi_alpha_shp2)
  
      p_alpha[s] ~ dbeta(p_alpha_shp1, p_alpha_shp2)
    }
    
    # binomial data model
    for (i in 1:N){
      # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
      %s
      
      # time to detection is a negative binomial process with state-specific shape and rate determined by covariates
      %s
      logit(p[i]) <- %s
      
      # model for censoring observed arrays due to not seeing into the future
      # whether we see an array is a bernouli process determined by
      # theta is 0 if site will never be developed (i.e. z[i] = 0) 
      #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
      
      d[i]~dbern(theta[i])
      theta[i] <-solar[i]*step(ttd[i] - tmax[i]) + (1-solar[i])
      
      # Expected data under current model
      ttd_sim[i] ~ dnegbin(p[statei[i]], 1)
      ttd_exp[i] <- p[statei[i]]/pow((1-p[statei[i]]), 2)
      chi2[i] <- pow(ttd[i] - ttd_exp[i], 2)/ttd_exp[i]
      chi2_sim[i] <- pow(ttd_sim[i] - ttd_exp[i], 2)/ttd_exp[i]
    }
    
    fit <- sum(chi2[])
    fit_sim <- sum(chi2_sim[])
    bpv <- step(fit_sim - fit)
    
    #Priors
    # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
    # to simulate gamma(0.0001, 0.0001) with no state effect
    psi ~ dunif(0, 1)
    p_alpha_mu ~ dunif(0, 1)
    p_alpha_v ~ dunif(0, 0.2)
    # psi_alpha_mu ~ dunif(0, 1) 
    # psi_alpha_v ~ dunif(0, 0.2)
    %s
  }
  "
  finalString = sprintf(
    baseString,
    binResponseString,
    # binPredictorString,
    contResponseString,
    contPredictorString,
    # binPriorsString,
    contPriorString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    zst <- rep(1, dat$N)
    ttdst <- dat$tmax+1 # creating some fake times greater than tmax by 1
    ttdst[!is.na(dat$ttd)] <- NA # this strictly overwrites the value of unobserved nodes!!!
    ls <- list(
      # z = zst,
      ttd = ttdst,
      # we will use mean and variance for a 'flat' beta dist beta(0.5, 0.5)
      # psi_alpha_mu = 0.5,
      # psi_alpha_v = 0.125,
      psi = dunif(1),
      p_alpha_mu = 0.5,
      p_alpha_v = 0.125
    )
    n <- length(ls)
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('p_beta', predictors, sep = "_")
    # names(ls)[(7+length(predictors)):length(ls)] <- paste('p_beta', predictors, sep = "_")
    
    return(ls)
  }
  
  return(inits)
  
}
