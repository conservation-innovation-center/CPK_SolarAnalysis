library(AHMbook)
library(tidyverse)
library(fitdistrplus)
test_dat <- simOccttd()

test_dat <- ttdPeregrine
# data set is a dataframe with 70 obersvations each at 1 of 38 sites
str(test_dat)

nobs <- length(test_dat$SiteNumber)
d <- as.numeric(is.na(test_dat$ttd)) # indicator for potential censoring. 1 = censored
tod <- calc_z(test_dat$MinOfDay)

win.data <- list(
  M=max(test_dat$SiteNumber), # the number of sites = 38
  site = test_dat$SiteNumber, # uniue identifier per site
  tod = tod, # standardized time of day of detection
  male = as.numeric(test_dat$sex)-1, # indicator for male individual
  ttd = test_dat$ttd, # time to detection - NA where no observation
  d = d, # indicator for censored observation. 1 where ttd is NA
  nobs = nobs,
  Tmax = test_dat$Tmax
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
    sexratio~dunif(0,1)
    
    #Likelihood
    for (i in 1:M){ #for each site
      z[i]~dbern(psi)
    }
    
    for (i in 1:nobs){
      ttd[i] ~ dweib(shape, lambda[i])
      log(lambda[i]) <- log(lambda.int) + alpha1*tod[i]+aplha2*pow(tod[i], 2)
      
      d[i]~dbern(theta[i])
      theta[i] <-z[site[i]]*step(ttd[i] - Tmax[i]) + (1+z[site[i]])
    }
    ")
zst <- rep(1, win.data$M)
ttdst <- rep(win.data$Tmax+1) # creating some fake times greater than tmax by 1
ttdst[win.data$d==0] <- NA #
inits <- function(){
  list(
    z=zst,
    ttd = ttdst,
    psi = runif(1),
    lambda.int = runif(1), # in the book this was runif(2) bc there was a lambda.int[1] and lambda.int[2]
    alpha=rnorm(1),
    alpha2=rnorm(1),
    shape = runif(1)
  )
}
params <- c('psi', 'lambda.int', 'alpha1', 'alpha2', 'z', 'shape')
out <- jags(
  data = win.data,
  inits = inits,
  parameters.to.save = params,
  model.file = 'jags/weibulltest.txt',
  n.chains = 3,
  n.adapt = 100,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 2)
test_dat[is.na(test_dat$ttd),]
