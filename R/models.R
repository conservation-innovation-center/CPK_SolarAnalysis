## TODO: FUnctionize each class of model

# Binomial
model{
  #Likelihood
  # binomial data model
  for (i in 1:N){
    solar[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(alpha) + beta1*impervious[i] + beta2*tree[i] + beta3*open[i]
  }
  
  #Priors
  alpha ~ dbeta(1, 1)
  beta1 ~ dnorm(0, 0.00001)
  beta2 ~ dnorm(0, 0.0001)
  beta3 ~ dnorm(0, 0.0001)
}

# Binomial w/Random state effect
model{
  #Likelihood
  # for each level of state
  for (s in 1:S){
    # two ways we can do this
    # 1. do a normal dist on a typical intercept
    alpha[s] ~ dnorm(alpha_mu, alpha_tau)
    # 2. do a beta dist on an intercept that will be alpha transformed as a 'baseline' probability
    alpha[s] ~ dbeta(alpha_shp1, alpha_shp2)
    # seems sketchy to set priors on beta hyperparameters. so use moment matching to achieve a desired dist
    #alpha_shp1 <- mu*alpha_shp2/(1+mu)
    alhpa_shp2 <- alpha_shp1*((1/mu) - 1)
    alpha_shp1 <- (((1-mu)/v)-1/mu )* pow(mu, 2)
  }
  # binomial data model
  for (i in 1:N){
    solar[i] ~ dbern(psi[i])
    # 1 from above
    logit(psi[i]) <- alpha[state[i]] + beta1*impervious[i] + beta2*tree[i] + beta3*open[i]
    # 2 from above
    logit(psi[i]) <- logit(alpha[state[i]]) + beta1*impervious[i] + beta2*tree[i] + beta3*open[i]
  }
  
  #Priors
  # 1 from above
  mu ~ dunif(0, 1)
  v ~ dgamma(1, 20)
  # 2 from above
  alpha_mu ~ dnorm(0, 0.0001)
  alpha_tau <- pow(sigma, -2)
  sigma ~ dunif(0, 1000)
  beta1 ~ dnorm(0, 0.00001)
  beta2 ~ dnorm(0, 0.0001)
  beta3 ~ dnorm(0, 0.0001)
}

# Weibull
# Note the response variable here is ttd NOT binary solar
# BUGS dweibull uses scale and rate parameters
model{
  #Likelihood
  # binomial data model
  for (i in 1:N){
    # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
    z[i] ~ dbern(psi[i])
    # probability of ever being developed linear fxn of covariates with state-specific intercept
    logit(psi[i]) <- psi_alpha + psi_beta1*impervious[i] + psi_beta2*tree[i] + psi_beta3*open[i]
    
    # time to detection is a weibull process with state-specific shape and rate determined by covariates
    ttd[i] ~ dweib(shape, rate[i])
    log(rate[i]) <- rate_alpha + rate_beta1*impervious[i] + rate_beta2*tree[i] + rate_beta3*open[i]
    
  }
  
  #Priors
  shape ~ dgamma(0.001, 0.001)
  psi_alpha ~ dnorm(0, 0.0001)
  psi_beta1 ~ dnorm(0, 0.0001)
  psi_beta2 ~ dnorm(0, 0.0001)
  psi_beta3 ~ dnorm(0, 0.0001)
  rate_alpha ~ dnorm(0, 0.0001)
  rate_beta1 ~ dnorm(0, 0.0001)
  rate_beta2 ~ dnorm(0, 0.0001)
  rate_beta3 ~ dnorm(0, 0.0001)
  
}

# WEibull with state random effects
model{
  #Likelihood
  # random intercept per state on binomial intercept and weibul shape
  for (s in 1:S){
    psi_alpha[s] ~ dnorm(psi_alpha_mu, psi_alpha_tau)
    # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
    # BUGS uses shape and rate parameterization of gamma
    # abov
    shape[s] ~ dgamma(shape_shape, shape_rate)
    shape_shape <- pow(mu, 2)/v
    shape_rate <- mu/v
    # the rate or scale parameter is [0, Inf] variability in ttd data
    rate_alpha[s] ~ dnorm(rate_alpha_mu, rate_alpha_tau)
  }
  
  # binomial data model
  for (i in 1:N){
    # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
    z[i] ~ dbern(psi[i])
    # probability of ever being developed linear fxn of covariates with state-specific intercept
    logit(psi[i]) <- psi_alpha[state[i]] + psi_beta1*impervious[i] + psi_beta2*tree[i] + psi_beta3*open[i]
    
    # time to detection is a weibull process with state-specific shape and rate determined by covariates
    ttd[i] ~ dweib(shape[state[i]], rate[i])
    log(rate[i]) <- rate_alpha[state[i]] + rate_beta1*impervious[i] + rate_beta2*tree[i] + rate_beta3*open[i]
  }
  
  #Priors
  # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
  # to simulate gamma(0.0001, 0.0001) with no state effect
  mu ~ dunif(0, 5)
  v ~ dunif(0, 1000) 
  rate_alpha_mu ~ dnorm(0, 0.0001)
  rate_alpha_tau <- pow(sigma, -2)
  psi_alpha_mu ~ dnorm(0, 0.0001) 
  psi_alpha_tau <- pow(sigma, -2) 
  sigma ~ dunif(0, 1000)
  rate_beta1 ~ dnorm(0, 0.0001)
  rate_beta2 ~ dnorm(0, 0.0001)
  rate_beta3 ~ dnorm(0, 0.0001)
  psi_beta1 ~ dnorm(0, 0.0001)
  psi_beta2 ~ dnorm(0, 0.0001)
  psi_beta3 ~ dnorm(0, 0.0001)
}

# Censored Weibull
# Note the response variable here is ttd NOT binary solar
model{
  #Likelihood
  # binomial data model
  for (i in 1:N){
    # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
    z[i] ~ dbern(psi[i])
    # probability of ever being developed linear fxn of covariates with state-specific intercept
    logit(psi[i]) <- psi_alpha + psi_beta1*impervious[i] + psi_beta2*tree[i] + psi_beta3*open[i]
    
    # time to detection is a weibull process with state-specific shape and rate determined by covariates
    ttd[i] ~ dweib(shape, rate[i])
    log(rate[i]) <- rate_alpha + rate_beta1*impervious[i] + rate_beta2*tree[i] + rate_beta3*open[i]
    
    # model for censoring observed arrays due to not seeing into the future
    # whether we see an array is a bernouli process determined by
    solar[i] ~ dbern(theta[i])
    # theta is 0 if site will never be developed (i.e. z[i] = 0) 
    #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
    theta[i] <- z[i]*step(ttd[i] - Tmax[i]) + (1 - z[i])
    
  }
  
  #Priors
  psi_alpha ~ dnorm(0, 0.0001)
  psi_beta1 ~ dnorm(0, 0.0001)
  psi_beta2 ~ dnorm(0, 0.0001)
  psi_beta3 ~ dnorm(0, 0.0001)
  rate_alpha ~ dnorm(0, 0.0001)
  rate_beta1 ~ dnorm(0, 0.0001)
  rate_beta2 ~ dnorm(0, 0.0001)
  rate_beta3 ~ dnorm(0, 0.0001)
  
}

# Censored Weibull with state random effects
model{
  #Likelihood
  # random intercept per state on binomial intercept and weibul shape
  for (s in 1:S){
    psi_alpha[s] ~ dnorm(psi_alpha_mu, psi_alpha_tau)
    # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
    # BUGS uses shape and rate parameterization of gamma
    # abov
    shape[s] ~ dgamma(shape_shape, shape_rate)
    shape_shape <- pow(mu, 2)/v
    shape_rate <- mu/v
    # the rate or scale parameter is [0, Inf] variability in ttd data
    rate_alpha[s] ~ dnorm(rate_alpha_mu, rate_alpha_tau)
  }
  
  # binomial data model
  for (i in 1:N){
    # zi = 'true' site occupancy - whether it will ever be developed. NOT what we observed
    z[i] ~ dbern(psi[i])
    # probability of ever being developed linear fxn of covariates with state-specific intercept
    logit(psi[i]) <- psi_alpha[state[i]] + psi_beta1*impervious[i] + psi_beta2*tree[i] + psi_beta3*open[i]
    
    # time to detection is a weibull process with state-specific shape and rate determined by covariates
    ttd[i] ~ dweib(shape[state[i]], rate[i])
    log(rate[i]) <- rate_alpha[state[i]] + rate_beta1*impervious[i] + rate_beta2*tree[i] + rate_beta3*open[i]
    
    # model for censoring observed arrays due to not seeing into the future
    # whether we see an array is a bernouli process determined by
    solar[i] ~ dbern(theta[i])
    # theta is 0 if site will never be developed (i.e. z[i] = 0) 
    #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
    theta[i] <- z[i]*step(ttd[i] - Tmax[i]) + (1 - z[i])
    
  }
  
  #Priors
  # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
  # to simulate gamma(0.0001, 0.0001) with no state effect
  mu ~ dunif(0, 5)
  v ~ dunif(0, 1000) 
  rate_alpha_mu ~ dnorm(0, 0.0001)
  rate_alpha_tau <- pow(sigma, -2)
  psi_alpha_mu ~ dnorm(0, 0.0001) 
  psi_alpha_tau <- pow(sigma, -2) 
  sigma ~ dunif(0, 1000)
  rate_beta1 ~ dnorm(0, 0.0001)
  rate_beta2 ~ dnorm(0, 0.0001)
  rate_beta3 ~ dnorm(0, 0.0001)
  psi_beta1 ~ dnorm(0, 0.0001)
  psi_beta2 ~ dnorm(0, 0.0001)
  psi_beta3 ~ dnorm(0, 0.0001)
}