library(jagsUI)
library(tidyverse)

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
modfile <- 'jags/binomial/test'
file.create(modfile)
#Write model to file
writeLines("
model{
  #Likelihood
  # for each level of state
  for (s in 1:S){
    # two ways we can do this
    # 1. do a normal dist on a typical intercept
    alpha[s] ~ dnorm(alpha_mu, alpha_tau)

  }
  # binomial data model
  for (i in 1:N){
    solar[i] ~ dbern(psi[i])
    # 1 from above
    logit(psi[i]) <- alpha[state[i]] + beta1*impervious[i] + beta2*tree[i] + beta3*open[i]
  }
  
  #Priors
  alpha_mu ~ dnorm(0, 0.0001)
  alpha_tau <- pow(sigma, -2)
  sigma ~ dunif(0, 1000)
  beta1 ~ dnorm(0, 0.00001)
  beta2 ~ dnorm(0, 0.0001)
  beta3 ~ dnorm(0, 0.0001)
}
", con=modfile)
######################################
## 3. Initialize Parameters ##
######################################
#Best to generate initial values using function
inits <- function(){
  list(
    alpha_mu =rnorm(1,0,1),
    sigma = runif(0, 1000),
    beta1=rnorm(1,0,1),
    beta2=rnorm(1,0,1),
    beta3=rnorm(1,0,1)
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
  map_chr(1:dat$S, function(x){sprintf('alpha[%s]', x)}),
  map_chr(1:length(cont), function(x){sprintf('beta%s', x)})
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


