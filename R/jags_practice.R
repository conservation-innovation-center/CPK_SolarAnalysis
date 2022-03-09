library(AHMbook)
library(tidyverse)
library(fitdistrplus)
test_dat <- simOccttd()
data('ttdPeregrine')
# data set is a dataframe with 70 obersvations each at 1 of 38 sites
str(ttdPeregrine)

nobs <- length(ttdPeregrine$SiteNumber)
d <- as.numeric(is.na(ttdPeregrine$ttd)) # indicator for potential censoring. 1 = censored
tod <- calc_z(ttdPeregrine$MinOfDay)

perData <- list(
  M=max(ttdPeregrine$SiteNumber), # the number of sites = 38
  site = ttdPeregrine$SiteNumber, # uniue identifier per site
  tod = tod, # standardized time of day of detection
  male = as.numeric(ttdPeregrine$sex)-1, # indicator for male individual
  ttd = ttdPeregrine$ttd, # time to detection - NA where no observation
  d = d, # indicator for censored observation. 1 where ttd is NA
  nobs = nobs,
  Tmax = ttdPeregrine$Tmax
)

cat(file = 'jags/weibulltest.txt',
    "
    model{
    
    #Priors
    psi ~ dunif(0,1)
    lambda.int~dgamma(0.001, 0.001)
    alpha1~dnorm(0, 0.001)
    alpha2~dnorm(0, 0.001)
    shape~dgamma(0.001, 0.001)
    
    #Likelihood
    for (i in 1:M){ #for each site
      z[i]~dbern(psi)
    }
    
    for (i in 1:nobs){
      ttd[i] ~ dweib(shape, lambda[i])
      log(lambda[i]) <- log(lambda.int) + alpha1*tod[i]+alpha2*pow(tod[i], 2)
      
      d[i]~dbern(theta[i])
      theta[i] <-z[site[i]]*step(ttd[i] - Tmax[i]) + (1-z[site[i]])
    }
    }
    ")

# these inits are needed b/c z is not part of data. we are estimating it
# with this data we start by assuming all sites occupied (z = 1) and all
# non-detections are due to censoring (where d = 0)
zst <- rep(1, perData$M)
ttdst <- rep(perData$Tmax+1) # creating some fake times greater than tmax by 1
ttdst[perData$d==0] <- NA # this strictly overwrites the value of unobserved nodes!!!
inits <- function(){
  list(
    z=zst,
    ttd = ttdst,
    psi = runif(1),
    lambda.int = runif(1), # in the book this was runif(2) bc there was a lambda.int[1] and lambda.int[2]
    alpha1=rnorm(1),
    alpha2=rnorm(1),
    shape = runif(1)
  )
}
params <- c('psi', 'lambda.int', 'alpha1', 'alpha2', 'z', 'shape')
out <- jags(
  data = perData,
  inits = inits,
  parameters.to.save = params,
  model.file = 'jags/weibulltest.txt',
  n.chains = 3,
  n.adapt = 100,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 2)
test_dat[is.na(test_dat$ttd),]
