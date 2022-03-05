library(jagsUI)
library(tidyverse)

## THIS CURRENTLY RAISES AN ERROR ABOUT Z BEING UNDEFINABLE
## QUESTIONS:
## Do we need to estimate z? It comes from the AHM book, but they had repeated
## observations at each site and z estimated true occupancy. We had

#Analyze Longley economic data in JAGS
#Number employed as a function of GNP
######################################
## 1. Collect and Package Data ##
######################################
#Load data 
cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
head(cleanDF)
#Separate data objects


#Input data objects must be numeric, and must be
#scalars, vectors, matrices, or arrays.

# first create a list of our continuous variables
cont <- list(
  impervious = cleanDF$impervious16,
  open = cleanDF$open16,
  tree = cleanDF$tree_cover16#,
  #ag = cleanDF$cultivated16
  #slope = df$slope,
  #road = df$road_dist,
  #line = df$line_dist,
  #housing = df$housing,
  #area = df$Shape_Area
)

# then create a list of our categorical variables
cat <- list(
  state = as.numeric(cleanDF$statef),
  Tmax = cleanDF$tmax
)

# list of our response variables
resp <- list(
  solar = cleanDF$solar,
  d = 1-cleanDF$solar, # for censored models, d is indicator of potential censorship
  ttd = cleanDF$ttd
)

# list of our constants
const <- list(
  N = nrow(cleanDF),
  S = length(unique(cleanDF$statef))
)

# constrict all continuous variables to [0-1]
# then combine with categorical variables
dat <- map(cont, standardize)%>%
  append(cat)%>%
  append(const)%>%
  append(resp)

######################################
## 2. Write model file ##
######################################
#Write a model in the BUGS language
#Generate model file directly in R
#(could also read in existing model file)
#Identify filepath of model file
modfile <- 'jags/weibull/test'
file.create(modfile)
#Write model to file
writeLines("
model{
  #Likelihood

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
    d[i] ~ dbern(theta[i])
    # theta is 0 if site will never be developed (i.e. z[i] = 0) 
    #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
    theta[i] <- z[i]*step(ttd[i] - Tmax[i]) + (1 - z[i])
  }
  # random intercept per state on binomial intercept and weibul shape
  for (s in 1:S){
    psi_alpha[s] ~ dnorm(psi_alpha_mu, psi_alpha_tau)
    # the shape parameter for weibull is [0, Inf] & indicates increasing, decreasing, or steady risk
    # BUGS uses shape and rate parameterization of gamma
    shape[s] ~ dgamma(shape_shape, shape_rate)
    # the rate or scale parameter is [0, Inf] variability in ttd data
    rate_alpha[s] ~ dnorm(rate_alpha_mu, rate_alpha_tau)
  }
  
  shape_shape <- pow(mu, 2)/v
  shape_rate <- mu/v
  
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
", con=modfile)
######################################
## 3. Initialize Parameters ##
######################################
#Best to generate initial values using function
ttdst <- rep(cleanDF$tmax+1) # creating some fake times greater than tmax by 1
ttdst[!is.na(cleanDF$ttd)] <- NA

inits <- function(){
  list(
    mu=runif(0, 5),
    v=runif(0, 1000),
    z = rep(0,dat$N),
    ttd = ttdst,
    rate_alpha_mu=1,
    sigma = 1000,
    rate_beta1=rnorm(1,0,1),
    rate_beta2=rnorm(1,0,1),
    rate_beta3=rnorm(1,0,1),
    psi_beta1=rnorm(1,0,1),
    psi_beta2=rnorm(1,0,1),
    psi_beta3=rnorm(1,0,1)
  )
}

#In many cases, JAGS can pick initial values automatically;
#you can leave argument inits=NULL to allow this.

######################################
## 4. Set parameters to monitor ##
######################################
#Choose parameters you want to save output for
#Only parameters in this list will appear in output object
#(deviance is added automatically if DIC=TRUE)
#List must be specified as a character vector
params <- c(
  'z',
  map_chr(1:length(cont), function(x){sprintf('psi_beta%s',x)}),
  map_chr(1:length(cont), function(x){sprintf('rate_beta%s', x)}),
  map_chr(1:dat$S, function(x){sprintf('psi_alpha[%s]',x)}),
  map_chr(1:dat$S, function(x){sprintf('rate_alpha[%s]', x)})
)
######################################
## 5. Run Analysis ##
######################################
#Call jags function; specify number of chains, number of adaptive iterations,
#the length of the burn-in period, total iterations, and the thin rate.
out <- jags(data = dat,
            inits = inits,
            parameters.to.save = params,
            model.file = modfile,
            n.chains = 3,
            n.adapt = 100,
            n.iter = 1000,
            n.burnin = 500,
            n.thin = 2)
#Arguments will be passed to JAGS; you will see progress bars
#and other information
#Examine output summary
out
#Look at output object elements
names(out)
#Plot traces and posterior densities
plot(out)
#Plot traces
traceplot(out)
#Update model another 1000 iterations
out <- update(out,n.iter = 1000)

## Follow up analyses

## Evaluate residuals

# first we compute the posterior mean response...?
mu <- out$mean$alpha + logistic(out$mean$beta1*dat$impervious + out$mean$beta2*dat$tree + out$mean$beta3*dat$open)


