### BINOMIAL ###
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
  responseString <- sprintf('%s[i] ~ dbin(p[i], tmax[i])T(0,1)', response)
  alphaString <- 'logit(alpha)'
  # create a string of predictor variables following beta_var*var[i]+...
  predictorString <- paste(alphaString, paste('beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  priorsString <- paste('beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  
  baseString <-"
  model{
    #Likelihood
    # binomial data model
    for (i in 1:N){
      %s
      logit(p[i]) <- %s
    
    # calculate expected values
    solar_exp[i] <- tmax[i]*p[i]
    solar_sim[i] ~ dbin(p[i], tmax[i])
    chi2[i] <- pow(solar[i] - solar_exp[i], 2)/solar_exp[i]
    chi2_sim[i] <- pow(solar_sim[i] - solar_exp[i], 2)/solar_exp[i]
    }
    
    fit <- mean(chi2[])
    fit_sim <- mean(chi2_sim[])
    bpv <- step(fit - fit_sim)
    
    #Priors
    alpha ~ dunif(0, 1)
    %s
  }"
  finalString = sprintf(baseString, responseString, predictorString, priorsString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ls <- list(
      alpha = runif(1,0,1)
    )
    n <- length(ls)
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('beta', predictors, sep = "_")
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
  responseString <- sprintf('%s[i] ~ dbern(p[i])', response)
  alphaString <- sprintf('logit(alpha[%s[i]])', random)
  # alphaString <- 'logit(alpha0)'
  # create a string of predictor variables following beta_var*var[i]+...
  predictorString <- paste(alphaString, paste('beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  priorsString <- paste('beta_', predictors, '~dnorm(0, 0.001)', sep = "", collapse = '\n    ')

  baseString <-"
  model{
    #Likelihood
    alpha_shp1 <- (((1-alpha_mu)/alpha_v) - (1/alpha_mu))*pow(alpha_mu, 2)
    alpha_shp2 <- (((1-alpha_mu)/alpha_v) - (1/alpha_mu))*alpha_mu*(1-alpha_mu)
    # for each level of random
    for (s in 1:S){
      alpha[s] ~ dbeta(alpha_shp1, alpha_shp2)
    }
    # binomial data model
    for (i in 1:N){
      %s
      logit(p[i]) <- %s
    
    # calculate expected values
    solar_exp[i] <- p[i]
    solar_sim[i] ~ dbern(p[i])
    chi2[i] <- pow(solar[i] - solar_exp[i], 2)/solar_exp[i]
    chi2_sim[i] <- pow(solar_sim[i] - solar_exp[i], 2)/solar_exp[i]
    }
    
    fit <- mean(chi2[])
    fit_sim <- mean(chi2_sim[])
    bpv <- step(fit - fit_sim)
    
    #Priors
    alpha_mu ~ dunif(0, 1)
    alpha_v ~ dunif(0, 0.2)
    %s
  }"
  finalString = sprintf(baseString, responseString, predictorString, priorsString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ls <- list(
      alpha_mu = 0.5,
      alpha_v = 0.08
      # alpha_mu = rnorm(1,0,1),
      # sigma = runif(0, 1000)
    )
    n <- length(ls)
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('beta', predictors, sep = "_")
    return(ls)
  }
  
  return(inits)
}

#' write a jags binomial model for variable importance evaluation
#' 
#' Given the names of response and predictor variables, create text representing a simple
#' binomial model in jags model format and write to a specified file
#' 
#' @param filename string specifying destination file
#' @param response string identifying name of response variable in jags data list
#' @param predictors character vector containing names of predictor varibles in jags data list
#' @param random string indentifying the name of variable used for random effects on intercept
#' @param out jags output containing mean and sd values for variables
#' @returns function to generate initial values
#' 
write_jags_binomial_varEval <- function(filename, response, predictors, random, out){
  file.create(filename)
  means <- unlist(out$mean[paste('beta', predictors, sep = '_')])
  sds <- unlist(out$sd[paste('beta', predictors, sep = '_')])
  
  responseString <- sprintf('%s[i] ~ dbin(p[i], tmax[i])T(0,1)', response)
  alphaString <- sprintf('logit(alpha[%s[i]])', random)
  # create a string of predictor variables following v_var*beta_var*var[i]+...
  predictorString <- paste(alphaString, paste('w_', predictors, '*beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  #
  wPrior <- paste('w_', predictors, '~dbern(0.5)', sep = "", collapse = '\n    ')
  # set our beta priors to the posterior distributions from the full model
  # beta_slope ~ dnorm(beta_slope_mu, beta_slope_tau)
  betaPriors <- paste('beta_', predictors, ' ~ dnorm(beta_', predictors,'_mu, beta_', predictors, '_tau)', sep = "", collapse = '\n    ')
  
  # beta_slope_mu <- w_slope*mean
  mus <- paste('beta_', predictors, '_mu <- w_', predictors, '*', means, sep = "", collapse = '\n    ')
  # beta_slope_tau <- w_slope*tau + (1-w_slope)*0.1
  taus <- paste('beta_', predictors, '_tau <- w_', predictors, '*', 1/(sds^2), '+ (1 - w_', predictors, ')*0.1', sep = "", collapse = '\n    ')
  
  baseString <-"
  model{
    #Likelihood
    alpha_shp1 <- (((1-alpha_mu)/alpha_v) - (1/alpha_mu))*pow(alpha_mu, 2)
    alpha_shp2 <- (((1-alpha_mu)/alpha_v) - (1/alpha_mu))*alpha_mu*(1-alpha_mu)
    # for each level of random
    for (s in 1:S){
      alpha[s] ~ dbeta(alpha_shp1, alpha_shp2)
    }
    # binomial data model
    for (i in 1:N){
      %s
      logit(p[i]) <- %s
    }
    
    #Priors
    alpha_mu ~ dunif(0, 1)
    alpha_v ~ dunif(0, 0.2)
    # w_alpha ~ dbern(0.5)
    # alpha0 ~ dnorm(0, 0.01)
    %s
    %s
    %s
    %s
  }"
  finalString = sprintf(baseString, responseString, predictorString, wPrior, betaPriors, mus, taus)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ls <- list(
      alpha_mu = 0.5,
      alpha_v = 0.08
    )
    n <- length(ls)
    for(var in 1:length(predictors)){
      ls <- append(ls, rnorm(1, means[var], sds[var]))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('beta', predictors, sep = "_")
    return(ls)
  }
  
  return(inits)
}

## NEGATIVE BINOMIAL
write_jags_binNbin <- function(filename, trunc, cens, cont_response, predictors, random){
  file.create(filename)
  # binResponseString <- sprintf('%s[i] ~ dbern(psi)', bin_response)
  binAlphaString <- 'logit(psi_alpha)'
  # create a string of predictor variables following beta_var*var[i]+...
  # binPredictorString <- paste(binAlphaString, paste('psi_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  contResponseString <- sprintf('%s[i] ~ dnegbin(p[i]*z[i]+0.00001, 1)T(%s[i],)', cont_response, random, trunc)
  contAlphaString <- 'logit(p_alpha)'
  contPredictorString <- paste(contAlphaString, paste('p_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  # censorString <- sprintf('%s[i]~dinterval(%s[i], tmax[i])', cens, cont_response)
  censorString1 <- sprintf('%s[i] ~ dbern(theta[i])', cens)
  censorString2 <- sprintf('theta[i] <- z[i] * step(%s[i] - tmax[i]) + (1 - z[i])', cont_response)
  
  # binPriorsString <- paste('psi_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  contPriorString <- paste('p_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  
  baseString <- "
  model{
    #Likelihood
    
    # Distribution Parameters for Random Efects on
    # 1. Negative Binomial probability p ~ dbeta(shape1, shape2)
    psi_shp2 <- psi_shp1*((1/psi_mu) - 1)
    psi_shp1 <- (((1-psi_mu)/psi_v)-1/psi_mu )* pow(psi_mu, 2)
    
    # random intercept per state on binomial intercept and neg binomial intercept
    for (s in 1:S){
      psi[s] ~ dbeta(psi_shp1, psi_shp2)
    }
    
    # for every observation
    for (i in 1:N){
      # model for suitability of a site to solar development, z[i]
      z[i] ~ dbern(psi[i])
      
      # binPredictorString
      
      # time to detection is a negative binomial process ignoring censoring
      # contResponseString
      %s
      # contPredictorString
      logit(p[i]) <- %s
      
      # model for censoring observed arrays due to not seeing into the future
      
      # OPTION 1:
      # whether we see an array is a bernouli process determined by
      # theta is 0 if site will never be developed (i.e. z[i] = 0) 
      #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
      
      %s
      %s
      
      # OPTION 2:
      # use jags dinterval to model censoring...?
      
      # Expected data under current model
      ttd_sim[i] ~ dnegbin(p[i], 1)T(l[i],)
      ttd_exp[i] <- (1-p[i])/p[i]
      chi2[i] <- pow(ttd[i] - ttd_exp[i], 2)/ttd_exp[i]
      chi2_sim[i] <- pow(ttd_sim[i] - ttd_exp[i], 2)/ttd_exp[i]
    }
    
    fit <- mean(chi2[])
    fit_sim <- mean(chi2_sim[])
    bpv <- step(fit_sim - fit)
    
    #Priors
    # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
    # to simulate gamma(0.0001, 0.0001) with no state effect
    psi_mu ~ dunif(0, 1)
    psi_v ~ dunif(0, 0.2)
    # psi_alpha_mu ~ dunif(0, 1) 
    # psi_alpha_v ~ dunif(0, 0.2)
    p_alpha ~ dnorm(0, 0.0001)
    # contPriorsString
    %s
  }
  "
  finalString = sprintf(
    baseString,
    # binResponseString,
    binPredictorString,
    contResponseString,
    # contPredictorString,
    # binPriorsString,
    censorString1,
    censorString2,
    # binPriorsString
    contPriorString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ttdst <- dat$tmax+1 # creating some fake times greater than tmax by 1
    ttdst[!is.na(dat$ttd)] <- NA # this strictly overwrites the value of unobserved nodes!!!
    zst <- rep(1, length(dat$ttd))
    ls <- list(
      z = zst,
      ttd = ttdst,
      # we will use mean and variance for a 'flat' beta dist beta(0.5, 0.5)
      # psi_alpha_mu = 0.5,
      # psi_alpha_v = 0.125,
      # psi = runif(1),
      # p_mu = 0.5,
      # p_v = 0.125
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

write_jags_nbin <- function(filename, trunc, cens, cont_response, predictors, random){
  file.create(filename)
  # binResponseString <- sprintf('%s[i] ~ dbern(psi)', bin_response)
  # binAlphaString <- sprintf('logit(psi_alpha[%s[i]])', random)
  # create a string of predictor variables following beta_var*var[i]+...
  # binPredictorString <- paste(binAlphaString, paste('psi_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  contResponseString <- sprintf('%s[i] ~ dnegbin(p[i], 1)T(%s[i],)', cont_response, trunc)
  contAlphaString <- sprintf('logit(p_alpha[%s[i]])', random)
  contPredictorString <- paste(contAlphaString, paste('p_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  censorString <- sprintf('%s[i]~dinterval(%s[i], tmax[i])', cens, cont_response)
  # binPriorsString <- paste('psi_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  contPriorString <- paste('p_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  
  baseString <- "
  model{
    #Likelihood
    
    # Distribution Parameters for Random Efects on
    # 1. Negative Binomial probability p ~ dbeta(shape1, shape2)
    p_alpha_shp2 <- p_alpha_shp1*((1/p_alpha_mu) - 1)
    p_alpha_shp1 <- (((1-p_alpha_mu)/p_alpha_v)-1/p_alpha_mu )* pow(p_alpha_mu, 2)
    
    # random intercept per state on binomial intercept and neg binomial intercept
    for (s in 1:S){
      p_alpha[s] ~ dbeta(p_alpha_shp1, p_alpha_shp2)
    }
    
    # for every observation
    for (i in 1:N){
      # time to detection is a negative binomial process with state-specific shape and rate determined by covariates
      %s
      logit(p[i]) <- %s
      
      # model for censoring observed arrays due to not seeing into the future
      # whether we see an array is a bernouli process determined by
      # theta is 0 if site will never be developed (i.e. z[i] = 0) 
      #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
      
      %s
      
      # Expected data under current model
      ttd_sim[i] ~ dnegbin(p[i], 1)
      ttd_exp[i] <- (1-p[i])/p[i]
      chi2[i] <- pow((ttd[i] - ttd_exp[i]), 2)/ttd_exp[i]
      chi2_sim[i] <- pow((ttd_sim[i] - ttd_exp[i]), 2)/ttd_exp[i]
    }
    
    fit <- mean(chi2[])
    fit_sim <- mean(chi2_sim[])
    bpv <- step(fit_sim - fit)
    
    #Priors
    # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
    # to simulate gamma(0.0001, 0.0001) with no state effect
    p_alpha_mu ~ dunif(0, 1)
    p_alpha_v ~ dunif(0, 0.2)
    # psi_alpha_mu ~ dunif(0, 1) 
    # psi_alpha_v ~ dunif(0, 0.2)
    # psi ~ dunif(0,1)
    %s
  }
  "
  finalString = sprintf(
    baseString,
    # binResponseString,
    # binPredictorString,
    contResponseString,
    contPredictorString,
    # binPriorsString,
    censorString,
    contPriorString)
  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ttdst <- dat$tmax+1 # creating some fake times greater than tmax by 1
    ttdst[!is.na(dat$ttd)] <- NA # this strictly overwrites the value of unobserved nodes!!!
    ls <- list(
      # z = zst,
      ttd = ttdst,
      # we will use mean and variance for a 'flat' beta dist beta(0.5, 0.5)
      # psi_alpha_mu = 0.5,
      # psi_alpha_v = 0.125,
      # psi = runif(1),
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

### WEIBULL ###

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
#' @param jagsDat list containing data for jags model output by \code{make_jags_dat}
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
write_jags_weibull <- function(filename, trunc, censor, response, predictors, random, jagsDat){
  file.create(filename)
  responseString <- sprintf('%s[i] ~ dweib(shape[%s[i]], rate[i])T(%s[i],)', response, random, trunc)
  alphaString <- sprintf('log(rate_alpha[%s[i]])', random)
  predictorString <- paste(alphaString, paste('rate_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
  
  censorString <- sprintf('%s[i] ~ dinterval(%s[i], tmax[i])', censor, response)
  
  simString <- sprintf('%s_sim[i] ~ dweib(shape[%s[i]], rate[i])T(%s[i],)', response, random, trunc)
  expString <- sprintf('%s_exp[i] <- pow(rate[i], -1/shape[%s[i]])*exp(loggam(1+ (1/shape[%s[i]])))', response, random, random)
  # binPriorsString <- paste('psi_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  priorString <- paste('rate_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
  
  baseString <- "
  model{
    #Likelihood
    # hyperparameters for random effects on 
    shape_shape <- pow(shape_mu, 2)/shape_v
    shape_rate <- shape_mu/shape_v
    
    rate_alpha_shape <- pow(rate_alpha_mu, 2)/rate_alpha_v
    rate_alpha_rate <- rate_alpha_mu/rate_alpha_v

    # random intercept per state on binomial intercept and weibul shape
    for (s in 1:S){
      # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
      # BUGS uses shape and rate parameterization of gamma
      shape[s] ~ dgamma(shape_shape, shape_rate)
      # the rate or scale parameter is [0, Inf] variability in ttd data
      rate_alpha[s] ~ dgamma(rate_alpha_shape, rate_alpha_rate)
    }
    
    # for each observation
    for (i in 1:N){
      # time to detection is a weibull process with state-specific shape and rate determined by covariates
      #responseString
      %s
      #predictorString
      log(rate[i]) <- %s
      
      # account for right censoring of ttd data - we cant look infinitely into the future
      #censorString
      %s

      # Expected data under current model
      #simString
      %s
      #expString
      %s
      chi2[i] <- pow((ttd[i] - ttd_exp[i]), 2)/ttd_exp[i]
      chi2_sim[i] <- pow((ttd_sim[i] - ttd_exp[i]), 2)/ttd_exp[i]
    }
    
    fit <- mean(chi2[])
    fit_sim <- mean(chi2_sim[])
    bpv <- step(fit_sim - fit)
    
    #Priors
    # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
    # to simulate gamma(0.0001, 0.0001) with no state effect
    shape_mu ~ dunif(0, 5)
    shape_v ~ dunif(0, 1000) 
    rate_alpha_mu ~ dunif(0, 5)
    rate_alpha_v ~ dunif(0, 1000)
    #priorString
    %s
  }"
  
  finalString = sprintf(
    baseString,
    responseString,
    predictorString,
    censorString,
    simString,
    expString,
    priorString)

  cat(finalString, sep = "")
  # write our model string to specified file
  cat(finalString, file = filename)
  
  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ttdst <- jagsDat$tmax+1 # creating some fake times greater than tmax by 1
    # ttdst <- weib_dat$tmax+1 # creating some fake times greater than tmax by 1
    ttdst[!is.na(jagsDat[[response]])] <- NA
    # ttdst[!is.na(weib_dat$ttd)] <- NA # this strictly overwrites the value of unobserved nodes!!!
    ls <- list(
      # we will chose mean and variance for gamma that corresponds to a 'flat' gamma(0.001, 0.001)
      ttd = ttdst,
      shape_mu= 2,
      shape_v= 1000,
      rate_alpha_mu = 2,
      rate_alpha_v = 1000
    )
    n <- length(ls)
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('rate_beta', predictors, sep = "_")
    # names(ls)[(5+length(predictors)):length(ls)] <- paste('rate_beta', predictors, sep = "_")
    
    return(ls)
  }
  
  return(inits)
  
}


# write_jags_binWeib <- function(filename, trunc, bin_response, cont_response, predictors, random, jagsDat){
#   file.create(filename)
#   binResponseString <- sprintf('%s[i] ~ dbern(psi)', bin_response)
#   # binAlphaString <- sprintf('logit(psi[i]) <- logit(psi_alpha[%s[i]])', random)
#   # create a string of predictor variables following beta_var*var[i]+...
#   # binPredictorString <- paste(binAlphaString, paste('psi_beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
#   
#   contResponseString <- sprintf('%s[i] ~ dweib(shape[%s[i]], rate[i])', cont_response, random)
#   rateString <- sprintf('rate[i] <- %s[i]*lambda[i]', bin_response)
#   contAlphaString <- sprintf('log(alpha[%s[i]])', random)
#   contPredictorString <- paste(contAlphaString, paste('beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + "), sep = ' + ')
# 
#   # censorString <- sprintf('%s[i] ~ dinterval(%s[i], tmax[i])', censor, cont_response)
#   
#   # binPriorsString <- paste('psi_beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
#   priorString <- paste('beta_', predictors, '~dnorm(0, 0.0001)', sep = "", collapse = '\n    ')
#   
#   baseString <- "
#   model{
#     #Likelihood
#     # hyperparameters for random effects on 
#     shape_shape <- pow(shape_mu, 2)/shape_v
#     shape_rate <- shape_mu/shape_v
#     
#     alpha_shape <- pow(alpha_mu, 2)/alpha_v
#     alpha_rate <- alpha_mu/alpha_v
#     
#     # psi_alpha_shp1 <- (((1-psi_alpha_mu)/psi_alpha_v) - (1/psi_alpha_mu))*pow(psi_alpha_mu, 2)
#     # psi_alpha_shp2 <- (((1-psi_alpha_mu)/psi_alpha_v) - (1/psi_alpha_mu))*psi_alpha_mu*(1-psi_alpha_mu)
# 
#     
#     # random intercept per state on binomial intercept and weibul shape
#     for (s in 1:S){
#       alpha[s] ~ dgamma(alpha_shape, alpha_rate)
#       # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
#       # BUGS uses shape and rate parameterization of gamma
#       shape[s] ~ dgamma(shape_shape, shape_rate)
#       # the rate or scale parameter is [0, Inf] variability in ttd data
#     }
#     
#     # for each observation
#     for (i in 1:N){
#       # the observation of a solar panel is a binomial process
#       # binResponseString
#       
#       solar[i] ~ dbern(psi)
#       
#       # censoring
#       d[i] ~ dbern(theta[i])
#       theta[i] <- (1-z[i]) + step(ttd[i] - tmax[i])*z[i]
#       
#       # time to detection is a weibull process with state-specific shape and rate determined by covariates
#       # contResponseString
#       ttd[i] ~ dweib(shape[statei[i]], rate[i])
#       # rateString
#       rate[i] <- z[i]*lambda[i]
#       # contPredictorString
#       log(lambda[i]) <- log(alpha[statei[i]]) + beta_slope*slope[i] + beta_tree_cover16*tree_cover16[i] + beta_road_dist[i]*road_dist[i]
#       
#       # Expected data under current model
#       # Expected value for joing weibull is lambda * psi
#       ttd_sim[i] ~ dweib(shape[statei[i]], lambda[i]*psi)
#       ttd_exp[i] <- pow(lambda[i]*psi, -1/shape[statei[i]])*exp(loggam(1 + (1/shape[statei[i]])))
#       chi2[i] <- pow(ttd[i] - ttd_exp[i], 2)/ttd_exp[i]
#       chi2_sim[i] <- pow(ttd_sim[i] - ttd_exp[i], 2)/ttd_exp[i]
#     }
#     
#     fit <- mean(chi2[])
#     fit_sim <- mean(chi2_sim[])
#     bpv <- step(fit_sim - fit)
#     
#     #Priors
#     # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
#     # to simulate gamma(0.0001, 0.0001) with no state effect
#     shape_mu ~ dunif(0, 5)
#     shape_v ~ dunif(0, 1000) 
#     alpha_mu ~ dunif(0, 5)
#     alpha_v ~ dunif(0, 1000)
#     psi ~ dunif(0,1)
#     # psi_alpha_mu ~ dunif(0, 1)
#     # psi_alpha_v ~ dunif(0, 0.2)
#     # priorString
#     beta_slope~dnorm(0, 0.0001)
#     beta_tree_cover16~dnorm(0, 0.0001)
#     beta_road_dist~dnorm(0, 0.0001)
#   }"
#   
#   # finalString = sprintf(
#   #   baseString,
#   #   binResponseString,
#   #   contResponseString,
#   #   rateString,
#   #   contPredictorString,
#   #   priorString)
#   
#   cat(baseString, sep = "")
#   # write our model string to specified file
#   cat(baseString, file = filename)
#   
#   # dynamically create inits function that will generate initial values for beta parameters for each predictor
#   inits <- function(){
#     ttdst <- jagsDat[[cont_response]]+1 # create some fake zeros where we didn't detect solar
#     ttdst[!is.na(jagsDat[[cont_response]])] <- NA # this strictly overwrites the value of unobserved nodes!!!
#     zst <- rep(1, length(jagsDat[[1]]))
#     ls <- list(
#       # we will chose mean and variance for gamma that corresponds to a 'flat' gamma(0.001, 0.001)
#       z = zst,
#       ttd = ttdst,
#       shape_mu= 2,
#       shape_v= 1000,
#       alpha_mu = 2,
#       alpha_v = 1000,
#       psi = runif(1,0,1)
#       # we will use mean and variance for a 'flat' beta dist beta(0.5, 0.5)
#       # psi_alpha_mu = 0.5,
#       # psi_alpha_v = 0.125
#     )
#     n <- length(ls)
#     for(var in predictors){
#       ls <- append(ls, rnorm(1,0,1))
#     }
#     names(ls)[(n+1):(n+length(predictors))] <- paste('beta', predictors, sep = "_")
#     # names(ls)[(5+length(predictors)):length(ls)] <- paste('rate_beta', predictors, sep = "_")
#     
#     return(ls)
#   }
#   
#   return(inits)
#   
# }

write_jags_binWeib <- function(filename, write = FALSE, jagsDat, predictors, bin_response, cont_response, censor, trunc, random = NULL){

  binString <- sprintf('%s[i] ~ dbern(psi)', bin_response)
  contBetaString <- paste('beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + ")
  betaPriors <- paste('beta_', predictors, '~dnorm(0, 0.01)', sep = "", collapse = '\n  ')
  
  if(is.null(random)){
    contString <- sprintf('%s[i] ~ dweib(shape, lambda[i])T(%s[i], )', cont_response, trunc)
    contAlphaString <- 'log(alpha)'
    hyperParams <- "#"
    randShape <- "#"
    otherPriors <- "
    shape ~ dgamma(0.0001, 0.0001)
    alpha ~ dgamma(0.0001, 0.0001)
    psi ~ dunif(0,1)"
    expected <- "ttd_sim[i] ~ dweib(shape, lambda[i])T(l[i], )
    ttd_exp[i] <- pow(lambda[i], -1/shape)*exp(loggam(1 + (1/shape)))"
  }else{
    contString <- sprintf('%s[i] ~ dweib(shape[%s[i]], lambda[i])T(%s[i], )', cont_response, random, trunc)
    contAlphaString <- sprintf('log(alpha[%s[i]])', random)
    hyperParams <- "
    shape_shape <- pow(shape_mu, 2)/shape_v
    shape_rate <- shape_mu/shape_v
    
    alpha_shape <- pow(alpha_mu, 2)/alpha_v
    alpha_rate <- alpha_mu/alpha_v"
    randShape <- "
    for (s in 1:S){
      shape[s] ~ dgamma(shape_shape, shape_rate)
      alpha[s] ~ dgamma(alpha_shape, alpha_rate)
    }"
    otherPriors <- "
    shape_mu ~ dunif(0,5)
    shape_v ~ dunif(0, 1000)
    alpha_mu ~ dunif(0,5)
    alpha_v ~ dunif(0, 1000)
    psi ~ dunif(0,1)"
    expected <- sprintf(
    "ttd_sim[i] ~ dweib(shape[%s[i]], lambda[i])T(l[i], )
    ttd_exp[i] <- pow(lambda[i], -1/shape[%s[i]])*exp(loggam(1 + (1/shape[%s[i]])))",
    random,
    random,
    random)
  }
  
  contPredictorString <- sprintf('%s + %s', contAlphaString, contBetaString)
  priorString <- paste(otherPriors, betaPriors, sep = '\n  ')
  baseString <- "
  model {
    # hyperparameters for random effects on 
    %s
    # mean shape from a preliminary run without random effects was 3.2
    %s
    #LIKELIHOOD
    for (i in 1:N){
      # likelihood solar
      # binString
      %s
      
      # likelihood time to detection
      # contString
      %s
      # contPredictorString
      log(lambda[i]) <- %s
      
      # censoring
      %s[i] ~ dbern(theta[i])
      theta[i] <- (1-%s[i]) * step(ttd[i] - tmax[i]) + (%s[i])
      
      # expected values
      %s
      chi2[i] <- pow(ttd[i] - ttd_exp[i], 2)/ttd_exp[i]
      chi2_sim[i] <- pow(ttd_sim[i] - ttd_exp[i], 2)/ttd_exp[i] 
    }
  
  fit <- mean(chi2[])
  fit_sim <- mean(chi2_sim[])

  # PRIORS
  %s
  }
  "
  string <- sprintf(
    baseString,
    hyperParams,
    randShape,
    binString,
    contString,
    contPredictorString,
    censor,
    bin_response, bin_response,
    expected,
    priorString
  )
  
  cat(string, sep = "")
  if(write){
    file.create(filename)
    cat(string, file = filename)
  }

  inits <- function(){
    ttdst <- jagsDat$tmax + 1
    ttdst[jagsDat[[censor]] == 0] <- NA
    if(is.null(random)){
      ls <- list(
        ttd = ttdst,
        alpha = 0.003,
        shape = 3.5,
        psi = runif(1)
      )
    }else{
      ls <- list(
        ttd = ttdst,
        alpha_mu = 0.003,
        shape_mu = 3.5,
        psi = runif(1)
      )
    }
    n <- length(ls)
    for(var in predictors){
      ls <- append(ls, rnorm(1,0,1))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('beta', predictors, sep = "_")
    return(ls)
  
  }
  return(inits)
}

#' write a jags joint binomial-weibull model for variable importance evaluation
#' 
#' Given the names of response and predictor variables, create text representing a simple
#' binomial model in jags model format and write to a specified file
#' 
#' @param filename string specifying destination file
#' @param jagsDat list of input data structured for jags model
#' @param bin_response string identifying name of binary response variable in jags data list
#' @param cont_response string identifying name of continuous response variable in jags data list
#' @param predictors character vector containing names of predictor variables in jags data list
#' @param censor string identifying name of binary censor indicator
#' @param trunc string identifying name of variable holding left truncation value
#' @param random string identifying the name of variable used for random effects on intercept
#' @param out jags output containing mean and sd values for variables
#' @returns function to generate initial values
write_jags_binWeib_varEval <- function(filename, write = FALSE, jagsDat, predictors, bin_response, cont_response, censor, trunc, random, out){
  means <- unlist(out$mean[paste('beta', predictors, sep = '_')])
  sds <- unlist(out$sd[paste('beta', predictors, sep = '_')])
  
  binString <- sprintf('%s[i] ~ dbern(psi)', bin_response)
  contString <- sprintf('%s[i] ~ dweib(shape[%s[i]], lambda[i])T(%s[i], )', cont_response, random, trunc)
  alphaString <- sprintf('log(alpha[%s[i]])', random)
  # alphaString <- 'w_alpha*log(alpha[%s[i]]) + (1-w_alpha)*log(alpha))'
  betaString <- paste('w_', predictors, '*beta_', predictors, '*', predictors, '[i]', sep = '', collapse = " + ")
  predictorString <- sprintf('%s + %s', alphaString, betaString)
  wPriors <- paste('w_', predictors, '~dbern(0.5)', sep = "", collapse = '\n    ')
  # set our beta priors to the posterior distributions from the full model
  # beta_slope ~ dnorm(beta_slope_mu, beta_slope_tau)
  betaPriors <- paste('beta_', predictors, ' ~ dnorm(beta_', predictors,'_mu, beta_', predictors, '_tau)', sep = "", collapse = '\n    ')
  
  # beta_slope_mu <- w_slope*mean
  mus <- paste('beta_', predictors, '_mu <- (1-w_', predictors, ')*', means, sep = "", collapse = '\n    ')
  # mus <- paste('beta_', predictors, '_mu <- ', means, sep = "", collapse = '\n    ')
  # beta_slope_tau <- w_slope*tau + (1-w_slope)*0.1
  taus <- paste('beta_', predictors, '_tau <- (1 - w_', predictors, ')*', 1/((2*sds)^2), '+ w_', predictors, '*0.05', sep = "", collapse = '\n    ')
  # taus <- paste('beta_', predictors, '_tau <- ', 1/(sds^2), sep = "", collapse = '\n    ')
  baseString <-"
  model {
    # hyperparameters for random effects on 
    
    # mean shape from a preliminary run without random effects was 3.2
    shape_shape <- pow(shape_mu, 2)/shape_v
    shape_rate <- shape_mu/shape_v
    
    alpha_shape <- pow(alpha_mu, 2)/alpha_v
    alpha_rate <- alpha_mu/alpha_v
    
    for (s in 1:S){
      shape[s] ~ dgamma(shape_shape, shape_rate)
      alpha[s] ~ dgamma(alpha_shape, alpha_rate)
    }
    for (i in 1:N){
      # likelihood solar
      # solar[i] ~ dbern(psi)
      %s
      
      # likelihood time to detection
      # ttd[i] ~ dweib(shape[statei[i]], lambda[i])T(l[i], )
      %s
      log(lambda[i]) <- %s
      
      # censoring
      %s[i] ~ dbern(theta[i])
      theta[i] <- (1-%s[i]) * step(ttd[i] - tmax[i]) + %s[i]

    }

  # shape ~ dgamma(0.0001, 0.0001)  
  shape_mu ~ dnorm(3.8, 1.5625)
  shape_v ~ dunif(0, 1000)
  # alpha ~ dgamma(0.0001, 0.0001)
  alpha_mu ~ dnorm(0.003, 40000)
  alpha_v ~ dunif(0, 1000)
  psi ~ dunif(0,1)
  %s
  %s
  %s
  %s
  }
  "
  finalString = sprintf(
    baseString,
    binString,
    contString,
    predictorString,
    censor,
    bin_response, bin_response,
    wPriors,
    betaPriors,
    mus,
    taus)
  
  cat(finalString, sep = "")
  # write our model string to specified file
  if(write){
    file.create(filename)
    cat(finalString, file = filename)
  }

  # dynamically create inits function that will generate initial values for beta parameters for each predictor
  inits <- function(){
    ttdst <- jagsDat$tmax + 1
    ttdst[jagsDat[[censor]] == 0] <- NA
    ls <- list(
      ttd = ttdst,
      alpha_mu = out$mean$alpha_mu,
      psi = runif(1),
      shape_mu = out$mean$shape_mu,
      alpha = out$mean$alpha,
      shape = out$mean$shape
    )
    n <- length(ls)
    for(var in 1:length(predictors)){
      ls <- append(ls, rnorm(1, means[var], sds[var]))
    }
    names(ls)[(n+1):(n+length(predictors))] <- paste('beta', predictors, sep = "_")
    return(ls)
  }
  
  return(inits)
}
