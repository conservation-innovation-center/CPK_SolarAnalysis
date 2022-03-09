source('R/functions.R')
library(tidyverse)
library(jagsUI)

# Load data 
cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
head(cleanDF)

# specify our continuous predictor variables and model file
continuous <- c('impervious16', 'open16', 'tree_cover16', 'slope')
modfile <- 'jags/weibull/test'

# Prepare data for jags
# Input data objects must be numeric, and must be
# scalars, vectors, matrices, or arrays.
dat <- make_jags_data(
  dataframe = cleanDF,
  random = c('statei'),
  continuous = continuous,
  categorical = c('tmax', 'statei'),
  response = c('solar', 'ttd', 'd')
)

# check the data structure to make sure things look as expected
str(dat)

# create the jags model, init, and params objects
inits <- write_jags_weibull_raneff(modfile, 'solar', 'ttd', continuous, 'statei')

# do some manual modification to our monitored parameters
params <- names(inits())
params <- params[grepl('beta', params)]%>%
# params <- append(params, paste('psi_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))
  append(paste('rate_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(paste('shape[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append('ttd_exp')

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

# first we compute the posterior mean response...?
pred_df <- as.data.frame(dat[continuous])
rateBetaMeans <- as.data.frame(out$mean[grepl('rate_beta', names(out$mean))])
rateBetas <- as.data.frame(out$sims.list[grepl('rate_beta', names(out$sims.list))])

simRates <- t(as.matrix(pred_df) %*% t(as.matrix(rateBetas)))

ttd_mu <- vapply(1:nrow(pred_df), function(x){
  i <- dat$statei[x]
  mu <- exp(simRates[,i] + log(out$sims.list$rate_alpha[,i]))
  return(mu)
},
FUN.VALUE = numeric(nrow(simRates))
)

  
# GOODNESS OF FIT
# We will compute a posterior p-value using a chi-square stat on ttd
# For each iteration get the chi2 of observed values and a randomly generated data point
# first need to compute expected value

# TODO: Create a generic function to compute expected value of weibull distribution given covariates
ttd_exp <- compute_weibull_response(out$sims.list$rate_alpha, as.data.frame(out$sims.list[grepl('rate_beta', names(out$sims.list))]), pred_df, out$sims.list$shape)

function(alpha, betas, predictors, shape){
  int <- log(alpha)
  rate <- exp(int + sum(betas*predictors))
  y <- meanWeibullJags(shape, rate)
  return(y)
}
