source('R/functions.R')
library(tidyverse)
library(jagsUI)

# Load data 
cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
head(cleanDF)

# specify our continuous predictor variables and model file
continuous <- c('impervious16', 'open16', 'tree_cover16', 'slope')
modfile <- 'jags/binomial/test'

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
inits <- write_jags_binomial_raneff(modfile, 'solar', continuous, 'statei')

# do some manual modification to our monitored parameters
params <- names(inits())
params <- params[grepl('beta', params)]
params <- append(params, paste('alpha[', 1:length(unique(dat$statei)), ']', sep = ""))

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

mu <- compute_binomial_response(dat, out, continuous, 'statei')
residuals <- dat$solar - mu
