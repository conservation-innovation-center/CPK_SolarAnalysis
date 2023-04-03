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



#` Calculate the shape and rate parameters of a Weibull distribution
#'
#' @description 
#' Given coefficient estimates from a jagsUI model specifying a linear
#' relationship between covariates #’ #’and the shape and/or rate parameters of 
#' a Weibull distribution, calculate the shape and rate at a series of observed
#' covariate values.
#' 
#' @param dataset JAGS input data from \code{make_jags_data}
#' @param params list of parameter estimates from jagsUI output
#' @param predictors vector of parameter names
#' 
#' @return data frame with columns 'shape' and 'rate' and nrows equal to length
#' of elements in \code{dataset}
#' 
#' @require dplyr
compute_weibull_params <- function(dataset, params, predictors){
  pred_df <- as.data.frame(dataset[predictors])
  betas <- as.data.frame(params[grepl('beta', names(params))])
  
  # simRates is a niter x nobs matrix. for mean values, niter is 1
  simRates <- t(as.matrix(pred_df) %*% t(as.matrix(betas)))
  
  # Define a function to calculate shape and rate parameters given
  # beta coefficients and predictor variables
  # if we are calculating the parameters at every iteration...
  if(dim(params$shape)[1] == dim(simRates)[1]){
    fx <- function(x){
      i <- dataset$statei[x]
      rate <- exp(simRates[,x] + log(params$alpha[,i]))
      shape <- params$shape[,i]
      return(data.frame(shape = shape, rate = rate))
    }
  }
  #...or using mean values from the MCMC chain
  else{
    fx <- function(x){
      i <- dataset$statei[x]
      rate <- exp(simRates[,x] + log(params$alpha[i]))
      shape <- params$shape[i]
      return(data.frame(shape = shape, rate = rate))
    }
  }
  # apply our function to each observation
  parameters <- lapply(1:nrow(pred_df), fx)%>%
    bind_rows()
  return(parameters)
}

#' Variance of Weibull distribution (JAGS)
#'
#' @param shape a vector of shape parameters
#' @param rate a vector of rate parameters as defined in JAGS
#'
#' @return vector of variance estimates of same length as shape
compute_weibull_variance <- function(shape, rate){
  scale <- jagsWeib2Rweib(shape, rate)
  variance <- varWeibull(shape, scale)
  return(variance)
}

#' Expected response of Weibull model
#'
#'
#' @param shape a vector of shape parameters
#' @param rate a vector of rate parameters as defined in JAGS
#'
#' @return vector of expected values of same length as shape
compute_weibull_response <- function(shape, rate){
  expected <- meanWeibullJags(shape, rate)
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
  y <- rate*shape*(x^(shape-1))*(exp(-rate*(x^shape)))
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

#' Weibull hazard function (JAGS)
#' 
#' @description
#' Returns the result of the hazard function for a Weibull distribution with parameters \code{shape} and \code{lambda}
#' as parameterized in JAGS
#' 
#' @param x vector of quantiles
#' @param shape shape parameter, v
#' @param lambda rate parameter, Lambds
#' 
#' @return vector of hazard rate of the specified Weibull distribution at x
hweibullJags <- function(x, shape, lambda){
  y  <- lambda*shape*(x^(shape-1))
  return(y)
}

meanWeibullJags <- function(shape, lambda){
  y <- (lambda^(-1/shape))*gamma(1+(1/shape))
  return(y)
}

#' Variance of the Weibull distribution
#' @description 
#' Cacluates the variance of a Weibull distribution as parameterized in R
#' with \code{shape} and \code{scale}
#' @param shape numeric vector of shape parameters for the weibull distribution
#' @param scale numeric vector of scale parameters fo rthe weibull distribution
#' @details 
#' The Weibull distribution with \code{shape} and \code{scale} parameters has
#' density f(x) = (shape/scale)(x/scale)^(shape-1)exp(-(x/scale)^shape) as 
#' calculated by \code{dweibull(shape, scale)}.
#' @return vector of variance values
varWeibull <- function(shape, scale){
  var <- (scale^2)*(gamma(1+(2/shape)) - (gamma(1+(1/shape))^2))
  return(var)
}

beta_mean <- function(alpha, beta){
  return(alpha/(alpha+beta))
}

beta_var <- function(alpha, beta){
  return((alpha*beta)/(((alpha+beta)^2) *(alpha + beta + 1)))
}
beta_shapes <-function(mean, var){
  psi_shp1 <- (((1-mean)/var) - (1/mean))*(mean^2)
  psi_shp2 <- (((1-mean)/var) - (1/mean))*mean*(1-mean)
  return(list(shp1 = psi_shp1, shp2 = psi_shp2))
}

compute_nbin_response <- function(dataset, params, predictors){
  pred_df <- as.data.frame(dataset[predictors])
  betas <- as.data.frame(params[grepl('p_beta', names(params))])
  # simRates is a nobs x niter. for mean values, niter is 1
  simRates <- t(as.matrix(pred_df) %*% t(as.matrix(betas)))
  fx <- function(x){
    i <- dat$statei[x]
    p <- logistic(simRates[,x] + params$p_alpha[i])
    mu <- p/(1-p)
    return(mu)
  }
  expected <- vapply(1:nrow(pred_df), fx, FUN.VALUE = numeric(nrow(simRates)))
  return(expected)
}

jagsWeib2Rweib <- function(shape, lambda){
  scale <- 1/(lambda^(1/shape))
  return(scale)
}

gamma_params <- function(mean, var){
  shape <- (mean^2)/var
  rate <- mean/var
  return(c('shape' = shape, 'rate' = rate))
}
