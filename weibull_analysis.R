source('R/functions.R')
library(tidyverse)
library(jagsUI)
library(plotly)

## SETUP
# Load data 
cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
head(cleanDF)

# specify our continuous predictor variables and model file
continuous <- c('impervious16',
                'open16',
                'tree_cover16',
                'cultivated16', 
                'slope',
                'GAP_Sts',
                'line_dist',
                'road_dist',
                'POPULATION',
                'lat')

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

# specify training hyperparameters
nthin = 2
nchain = 3
niter = 5000
nburn = 500
nadapt = 100

### WEIBULL ANALYSIS ###
# create model file
modfile <- 'jags/weibull/full'
file.create(modfile)

# create the jags model, init, and params objects
inits <- write_jags_weibull_raneff(modfile, 'solar', 'ttd', continuous, 'statei')

# do some manual modification to our monitored parameters
params <- names(inits())
params <- params[grepl('beta', params)]%>%
# params <- append(params, paste('psi_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))
  append(paste('rate_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(paste('shape[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append('ttd_sim')

#Call jags function; specify number of chains, number of adaptive iterations,
#the length of the burn-in period, total iterations, and the thin rate.
out <- jags(data = dat,
            inits = inits,
            parameters.to.save = params,
            model.file = modfile,
            n.chains = nchain,
            n.adapt = nadapt,
            n.iter = niter,
            n.burnin = nburn,
            n.thin = nthin)

# GOODNESS OF FIT
# We will compute a posterior p-value using a chi-square stat on ttd
# For each iteration get the chi2 of observed values and a randomly generated data point

# this is a 750 x 5127 matrix of expected ttd at each site for each iteration
ttd_exp <- compute_weibull_response(dat, out$sims.list, continuous)
# quick check to see if we're getting reasonable results
hist(ttd_exp, breaks = 100, xlim = c(0,20))

# compute chi-square of observed data
chi2 <- t((apply(ttd_exp, 1, function(x){((dat$ttd - x)^2)/x})))
fit <- apply(chi2, 1, function(x){mean(x, na.rm = TRUE)})
chi2Sim <- ((out$sims.list$ttd_sim - ttd_exp)^2)/ttd_exp
fitSim <- apply(chi2Sim, 1, function(x){mean(x, na.rm = TRUE)})

# compare the chisquare distributions of the observed and simulated data
hist(chi2, col = 'white', breaks = 100)
hist(chi2Sim, col = 'grey', breaks = 1000, add = TRUE)
# plot_ly(alpha = 0.6)%>%
#   add_histogram(
#     x = chi2,
#     marker = list(
#       color = 'blue'
#     )
#   )%>%
#   add_histogram(
#     x = chi2Sim,
#     marker = list(
#       color = 'orange'
#     )
#   )%>%
#   layout(
#     barmode = 'overlay'
#   )
# 
# rm(chi2, chi2New)

plot(fit, fitSim, xlim = c(0,1), ylim = c(0,1), xlab = 'Fit actual', ylab = 'Fit simulated')
abline(0,1,lwd=2)

# RESIDUALS ANALYSIS
ttd_exp <- compute_weibull_response(dat, out$mean, continuous)

resid <- data.frame(resid = dat$ttd - ttd_exp, obs=dat$ttd, exp = ttd_exp)

# plot residuals by actual value
plot(resid$obs, resid$resid, xlab = 'observed value', ylab = 'residual')
plot(resid$resid)

#### NEGATIVE BINOMIAL ANALYSIS

modfile <- 'jags/negbin/full'

# create the jags model, init, and params objects
inits <- write_jags_nbin_raneff(modfile, 'solar', 'ttd', continuous, 'statei')

# do some manual modification to our monitored parameters
params <- names(inits())
params <- params[grepl('beta', params)]%>%
  # append(paste('psi_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(paste('p_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(c('fit', 'fit_sim', 'bpv', 'psi'))

#Call jags function; specify number of chains, number of adaptive iterations,
#the length of the burn-in period, total iterations, and the thin rate.
out <- jags(data = dat,
            inits = inits,
            parameters.to.save = params,
            model.file = modfile,
            n.chains = nchain,
            n.adapt = nadapt,
            n.iter = niter,
            n.burnin = nburn,
            n.thin = nthin)
