source('R/functions.R')
source('R/models.R')
library(tidyverse)
library(jagsUI)
library(plotly)
library(sf)

## SETUP
# Load data 
cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
head(cleanDF)

# transform predictor variables based on histograms
cleanDF$slope <- sqrt(cleanDF$slope)
cleanDF$line_dist <- log(cleanDF$line_dist+0.1)
cleanDF$road_dist <- log(cleanDF$road_dist+0.1)
cleanDF$GAP_Sts <- log(6-cleanDF$GAP_Sts)
cleanDF$pdensity <- log(cleanDF$pdensity + 0.1)

# Summary stats
filter(cleanDF, year < 2022)%>%group_by(statef)%>%summarize(sd = sd(area/100000), mn = mean(area/100000))
glm(data = cleanDF[cleanDF$year < 2022,], formula = I(area/100000) ~ state -1, family = gaussian(link = 'log'))
# Temporal analysis
timeDF <- readRDS(file = 'www/data/timelinedat.rds')
timeDF[timeDF$area == 0, c(2,4,5)] <- NA
cleanTime <- timeDF[!is.na(timeDF$area),]

timeDat <- list(
  area = cleanTime$std_area*100,
  N = nrow(cleanTime),
  year = as.numeric(cleanTime$year) - 2016,
  state = as.numeric(as.factor(cleanTime$State)),
  S = 5
)

# specify training hyperparameters
testMCMCparams <- list(
  nthin = 2,
  nchain = 3,
  niter = 1000,
  nburn = 500,
  nadapt = 100
)

MCMCparams <- list(
  nthin = 2,
  nchain = 3,
  niter = 10000,
  nburn = 2000,
  nadapt = 500
)

## TEMPORAL ANALYSIS
modfile <- 'jags/normal/full_random'

timeDat <- list(
  area = cleanTime$std_area*100,
  N = nrow(cleanTime),
  year = as.numeric(cleanTime$year) - 2016,
  state = as.numeric(as.factor(cleanTime$State)),
  S = 5
)

params <- paste('alpha[', 1:timeDat$S, ']', sep = "")%>%
  append(paste('beta[', 1:timeDat$S, ']', sep = ""))%>%
  append(c('area_exp', 'fit', 'fit_sim'))

outTime <- jags(data = timeDat,
                parameters.to.save = params,
                model.file = modfile,
                n.chains = MCMCparams$nchain,
                n.adapt = MCMCparams$nadapt,
                n.iter = MCMCparams$niter,
                n.burnin = MCMCparams$nburn,
                n.thin = MCMCparams$nthin)

plot(outTime$sims.list$fit, outTime$sims.list$fit_sim)

# pairwise test for posterior overlap
sum(outTime$sims.list$beta[,2] > outTime$sims.list$beta[,5])/12000

# BIODIVERSITY ANALYSIS
modfile <- 'jags/normal/biodiversity'

bioDat <- list(
  richness = round(biod_long$richness[biod_long$taxon == 'all']),
  N = nrow(biod_long[biod_long$taxon == 'all',]),
  solar = biod_long$solar[biod_long$taxon == 'all'],
  state = as.numeric(as.factor(biod_long$state[biod_long$taxon == 'all'])),
  S = 5
)

params <- paste('alpha[', 1:bioDat$S, ']', sep = "")%>%
  append(paste('beta[', 1:bioDat$S, ']', sep = ""))%>%
  append(c('richness_exp', 'fit', 'fit_sim'))

outBio <- jags(data = bioDat,
                parameters.to.save = params,
                model.file = modfile,
                n.chains = MCMCparams$nchain,
                n.adapt = MCMCparams$nadapt,
                n.iter = MCMCparams$niter,
                n.burnin = MCMCparams$nburn,
                n.thin = MCMCparams$nthin)

plot(outBio$sims.list$fit, outBio$sims.list$fit_sim)
abline(0,1,lwd=2)

# continue running the model until convergence
outBio<- update(
  object = outBio,
  n.iter = MCMCparams$niter,
  n.thin = MCMCparams$nthin
)

plot(outBio$sims.list$fit, outBio$sims.list$fit_sim)
abline(0,1,lwd=2)

# specify our continuous predictor variables and model file
test <- c('slope','tree_cover16','road_dist')
full <- c('impervious16',
          'open16',
          'tree_cover16',
          'cultivated16',
          'ssurgo',
          'slope', #square root(x)
          'GAP_Sts', #log(6-x)
          'line_dist', #square root(x)
          'road_dist', #square root(x)
          'pdensity',
          'income',
          'lat')

names(full) <- full

## DATA PREP - Prepare data for jags
# Input data objects must be numeric, and must be
# scalars, vectors, matrices, or arrays.
dat <- make_jags_data(
  dataframe = cleanDF[cleanDF$year > 2016,],
  random = c('statei'),
  continuous = full,
  categorical = c('tmax', 'statei'),
  response = c('ttd', 'd', 'l', 'solar')
)

# weib_dat <- dat
# weib_dat$l <- dat$l +1
# weib_dat$ttd <- dat$ttd +1
# weib_dat$tmax <- dat$tmax +1

test_dat <- make_jags_data(
  dataframe = cleanDF[cleanDF$year > 2016,],
  random = c('statei'),
  continuous = test,
  categorical = c('tmax', 'statei'),
  response = c('ttd', 'd', 'l', 'solar')
)

# check the data structure to make sure things look as expected
str(dat)

### BERNOULLI ANALYSIS ###
modfile <- 'jags/bernoulli/full_random'
file.create(modfile)
inits <- write_jags_binomial_raneff(modfile, 'solar', full, 'statei')
params <- names(inits())
params <- params[grepl('beta', params)]%>%
  append(paste('alpha[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(c('solar_exp', 'fit', 'fit_sim'))

outBern <- jags(data = dat,
               inits = inits,
               parameters.to.save = params,
               model.file = modfile,
               n.chains = MCMCparams$nchain,
               n.adapt = MCMCparams$nadapt,
               n.iter = MCMCparams$niter,
               n.burnin = MCMCparams$nburn,
               n.thin = MCMCparams$nthin)

plot(outBern$sims.list$fit, outBern$sims.list$fit_sim, xlim = c(0,20), ylim = c(0,20))
abline(0,1,lwd=2)

# continue running the model until convergence
outBern<- update(
  object = outBern,
  n.iter = MCMCparams$niter,
  n.thin = MCMCparams$nthin
)

saveRDS(outBern, file= 'data/output/bernoulli/full_random.rds')

resid <- data.frame(obs = dat$solar, exp = outBin$mean$solar_exp, n = dat$tmax)%>%
  mutate(resid = obs - exp)%>%
  mutate(pearson = resid/sqrt(exp*(n- exp)/n))

plot(resid$obs, resid$pearson)

outBin <- readRDS(file = 'data/output/binomial/binomial_full.rds')

## Variable selection
## we will use the approach of fitting a model with indicator variables for each predictor.
## Predictor priors are set to the posterior from the full model run

outBern <- readRDS(file = 'data/output/bernoulli/full_random.rds')

modfile <- 'jags/bernoulli/full_variableEval'
file.create(modfile)
inits <- write_jags_binomial_varEval(modfile, 'solar', full, 'statei', outBern)
params <- names(inits())[grepl('beta', full)]%>%
  append(c(paste('w', full, sep = '_')))

evalBern <- jags(data = dat,
               inits = inits,
               parameters.to.save = params,
               model.file = modfile,
               n.chains = MCMCparams$nchain,
               n.adapt = MCMCparams$nadapt,
               n.iter = MCMCparams$niter,
               n.burnin = MCMCparams$nburn,
               n.thin = MCMCparams$nthin)

saveRDS(evalBern, file = 'data/output/bernoulli/full_eval.rds')

mod <- paste(evalBern$sims.list$w_impervious16,
             evalBern$sims.list$w_open16,
             evalBern$sims.list$w_tree_cover16,
             evalBern$sims.list$w_cultivated16,
             evalBern$sims.list$w_ssurgo,
             evalBern$sims.list$w_slope,
             evalBern$sims.list$w_GAP_Sts,
             evalBern$sims.list$w_line_dist,
             evalBern$sims.list$w_road_dist,
             evalBern$sims.list$w_POPULATION,
             evalBern$sims.list$w_lat)

as.data.frame(table(mod))%>%
  mutate(prob = Freq/sum(Freq))%>%
  arrange(Freq)

table(mod)/sum(table(mod))
## Without Random Effects ##
modfile <- 'jags/binomial/full_norandom'
file.create(modfile)
inits <- write_jags_binomial(modfile, 'solar', full)
params <- names(inits())
params <- params[grepl('beta', params)]%>%
  append(c('p', 'solar_exp', 'fit', 'fit_sim'))

outBin <- jags(data = dat,
               inits = inits,
               parameters.to.save = params,
               model.file = modfile,
               n.chains = MCMCparams$nchain,
               n.adapt = MCMCparams$nadapt,
               n.iter = MCMCparams$niter,
               n.burnin = MCMCparams$nburn,
               n.thin = MCMCparams$nthin)

plot(outBin$sims.list$fit, outBin$sims.list$fit_sim)
abline(0,1,lwd=2)
resid <- data.frame(obs = dat$solar, exp = outBin$mean$solar_exp, n = dat$tmax)%>%
  mutate(resid = obs - exp)%>%
  mutate(pearson = resid/sqrt(exp*(n- exp)/n))

plot(resid$obs, resid$pearson)


### BIN_WEIB MIXTURE
weib_dat <- test_dat
weib_dat$l <- test_dat$l +1
weib_dat$ttd <- test_dat$ttd +1
weib_dat$tmax <- test_dat$tmax +1

modfile <- 'jags/binWeib/test_truncated_censored_no'
file.create(modfile)
inits <- write_jags_binWeib(
  filename = modfile,
  write = TRUE,
  jagsDat = test_dat,
  predictors = test,
  bin_response = 'solar',
  cont_response = 'ttd',
  censor = 'd',
  trunc = 'l')#,
  # random = 'statei')

params <- c("alpha", "shape" ,"beta_slope","beta_tree_cover16","beta_road_dist") 
#Call jags function; specify number of chains, number of adaptive iterations,
#the length of the burn-in period, total iterations, and the thin rate.

# alpha and shape parameters were having a hard time converging,
# so we run a test and use the posteriors for these to inform priors for a run
# with all predictors

testout <- jags(data = test_dat,
            inits = inits,
            parameters.to.save = params,
            model.file = modfile,
            n.chains = MCMCparams$nchain,
            n.adapt = MCMCparams$nadapt,
            n.iter = MCMCparams$niter,
            n.burnin = MCMCparams$nburn,
            n.thin = MCMCparams$nthin)

alpha_mu <- testout$mean$alpha # 0.003332244
alpha_sd <- testout$sd$alpha # 0.0006388908
shape_mu <- testout$mean$shape # 3.9604
shape_sd <- testout$sd$shape # 0.1287315

modfile <- 'jags/binWeib/full_truncated_censored_no'
# manually modify this file to set priors on alpha and shape equal to
# gamma distributions with equivalent mean and variance
alpha_gamma_shape <- alpha_mu^2/(5*alpha_sd)^2
alpha_gamma_rate <- alpha_mu/(5*alpha_sd)^2
shape_gamma_shape <- shape_mu^2/(5*shape_sd)^2
shape_gamma_rate <- shape_mu/(5*shape_sd)^2

file.create(modfile)
inits <- write_jags_binWeib(
  filename = modfile,
  write = FALSE,
  jagsDat = dat,
  predictors = full,
  bin_response = 'solar',
  cont_response = 'ttd',
  censor = 'd',
  trunc = 'l')#,
# random = 'statei')

# do some manual modification to our monitored parameters
params <- names(inits())
params <- params[grepl('beta', params)]%>%
  # params <- append(params, paste('psi_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))
  # append(paste('alpha[', 1:length(unique(weib_dat$statei)), ']', sep = ""))%>%
  # append(paste('shape[', 1:length(unique(weib_dat$statei)), ']', sep = ""))%>%
  append(c('ttd_exp', 'fit', 'fit_sim', 'shape', 'alpha'))

out <- jags(data = dat,
                inits = inits,
                parameters.to.save = params,
                model.file = modfile,
                n.chains = MCMCparams$nchain,
                n.adapt = MCMCparams$nadapt,
                n.iter = MCMCparams$niter,
                n.burnin = MCMCparams$nburn,
                n.thin = MCMCparams$nthin)

out <- update(
  object = out,
  n.iter = MCMCparams$niter,
  n.thin = MCMCparams$nthin
)
# convergence!!

# file naming scheme {dataset}_{truncated}_{censored}_{random?}
saveRDS(out, file = 'data/output/binWeib/full_truncated_censored_no.rds')

# GOODNESS OF FIT
out <- readRDS(file = 'data/output/binWeib/full_truncated_censored_no.rds')
# We will compute a posterior p-value using a chi-square stat on ttd
# For each iteration get the chi2 of observed values and a randomly generated data point
hist(out$mean$ttd_exp, breaks = 500, xlim = c(0,20))
plot(out$sims.list$fit, out$sims.list$fit_sim, xlim = c(0, 20), ylim = c(0,20))
abline(0,1,lwd=2)

# RESIDUALS ANALYSIS
params <- compute_weibull_params(weib_dat, out$mean, full)
ttd_exp <- compute_weibull_response(weib_dat, out$mean, full)
ttd_var <- compute_weibull_variance(weib_dat, out$mean, full)
ttd_sim <- purrr::map2_dbl(params$shape, params$rate, function(x,y){
  rweibull(1, x, jagsWeib2Rweib(x, y))
}
  )

# to calculated standardized residuals, we divide each residual by the sqrt of the variance
# of the weibull distribution

resid <- data.frame(obs=weib_dat$ttd,
                    exp = ttd_exp,
                    var = ttd_var,
                    sim = ttd_sim)%>%
  mutate(resid = obs - exp,
         resid_sim = sim - exp,
         stresid = resid/sqrt(var),
         stresid_sim = resid-sim/sqrt(var))

# plot residuals by actual value
plot_ly(x = resid$obs, y = resid$resid, type = 'scatter', mode = 'markers')%>%
  layout(
    xaxis = list(
      title = 'Observed ttd'
    ),
    yaxis = list(
      title = 'Studentized residual'
    )
  )


plot(resid$obs, resid$stresid, xlab = 'observed value', ylab = 'residual')
hist(resid$stresid, breaks = 50)
qqplot(resid$stresid, resid$stresid_sim)
abline(0,1)

## ADD RANDOM EFFECTS
# now that we have converged beta estimates for all parameters, let's try to get
# random effects on alpha and shape
out <- readRDS(file = 'data/output/binWeib/full_truncated_censored_no.rds')
modfile <- 'jags/binWeib/full_truncated_censored_both'
# manually modify this file to set priors on betas equal to posterior from model with no random effects
file.create(modfile)
inits <- write_jags_binWeib(
  filename = modfile,
  write = FALSE,
  jagsDat = dat,
  predictors = full,
  bin_response = 'solar',
  cont_response = 'ttd',
  censor = 'd',
  trunc = 'l',
  random = 'statei')

# do some manual modification to our monitored parameters
params <- names(inits())
params <- params[grepl('beta', params)]%>%
  # params <- append(params, paste('psi_alpha[', 1:length(unique(dat$statei)), ']', sep = ""))
  append(paste('alpha[', 1:length(unique(weib_dat$statei)), ']', sep = ""))%>%
  append(paste('shape[', 1:length(unique(weib_dat$statei)), ']', sep = ""))%>%
  append(c('alpha_mu', 'shape_mu'))

outBinWeib <- jags(data = dat,
            inits = inits,
            parameters.to.save = params,
            model.file = modfile,
            n.chains = MCMCparams$nchain,
            n.adapt = MCMCparams$nadapt,
            n.iter = MCMCparams$niter,
            n.burnin = MCMCparams$nburn,
            n.thin = MCMCparams$nthin)

# continue running the model until convergence
outBinWeib<- update(
  object = outBinWeib,
  n.iter = MCMCparams$niter,
  n.thin = MCMCparams$nthin
)

saveRDS(outBinWeib, file = 'data/output/binWeib/fulllog_truncated_censored_both.rds')
outBinWeib <- readRDS(file = 'data/output/binWeib/fulllog_truncated_censored_both.rds')
# RESIDUALS ANALYSIS
params <- compute_weibull_params(dat, outBinWeib$mean, full)
ttd_exp <- compute_weibull_response(params$shape, params$rate)
ttd_var <- compute_weibull_variance(params$shape, params$rate)
ttd_sim <- purrr::map2_dbl(params$shape, params$rate, function(x,y){
  rweibull(1, x, jagsWeib2Rweib(x, y))
}
)

# to calculated standardized residuals, we divide each residual by the sqrt of the variance
# of the weibull distribution

resid <- data.frame(obs=dat$ttd,
                    shape = params$shape,
                    rate = params$rate,
                    exp = ttd_exp-1,
                    var = ttd_var,
                    sim = ttd_sim)%>%
  mutate(resid = obs - exp,
         resid_sim = sim - exp,
         stresid = resid/sqrt(var),
         stresid_sim = resid-sim/sqrt(var),
         fit = resid/(exp^2),
         fit_sim = resid_sim/(exp^2)
         )

plot(resid$obs, resid$stresid, ylim = c(-4, 4), xlab = 'Observed value', ylab = 'Standardized residual')
plot(resid$fit, resid$fit_sim, xlim = c(0, 2), ylim = c(0, 2), xlab = 'Observed data', ylab = 'Simulated data')
abline(0,1)

varProbs <- map_dfc(full, function(var){
  beta <- paste('beta', var, sep='_')
  p <- sum(outBinWeib$sims.list[[beta]] > 0)/length(outBinWeib$sims.list[[beta]])
  u <- outBinWeib$mean[[beta]]
  sd <- outBinWeib$sd[[beta]]
  return(c(p, u, sd))
})

### VARIABLE EVALUATION
outBinWeib <- readRDS(file = 'data/output/binWeib/fulllog_truncated_censored_both.rds')
modfile <- 'jags/binWeib/fulllog_truncated_censored_both_GVS'
file.create(modfile)
inits <- write_jags_binWeib_varEval(
  filename = modfile,
  write = FALSE,
  jagsDat = dat,
  predictors = full,
  bin_response = 'solar',
  cont_response = 'ttd',
  censor = 'd',
  trunc = 'l',
  random = 'statei',
  out = outBinWeib)

params <- names(inits())
params <- params[grepl('beta', params)]%>%
  append(c(paste('w', full, sep = '_')))%>%
  append(paste('alpha[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(paste('shape[', 1:length(unique(dat$statei)), ']', sep = ""))%>%
  append(c('alpha_mu', 'shape_mu', 'alpha_v', 'shape_v'))

gvsBinWeib <- jags(data = dat,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = modfile,
                 n.chains = MCMCparams$nchain,
                 n.adapt = MCMCparams$nadapt,
                 n.iter = MCMCparams$niter,
                 n.burnin = MCMCparams$nburn,
                 n.thin = MCMCparams$nthin)

saveRDS(gvsBinWeib, file = 'data/output/binweib/fulllog_truncated_censored_both_GVS.rds')

gvsBinWeib<- update(
  object = gvsBinWeib,
  n.iter = MCMCparams$niter,
  n.thin = MCMCparams$nthin
)

gvsBinWeib <- readRDS(file = 'data/output/binweib/fulllog_truncated_censored_both_GVS.rds' )

as.data.frame(table(mod))%>%
  mutate(prob = Freq/sum(Freq))%>%
  arrange(Freq)

gvsProbs <- map_dfc(full, function(var){
  param <- paste('w', var, sep = "_")
  beta <- paste('beta', var, sep='_')
  p <- sum(gvsBinWeib$sims.list[[param]] == 1)/length(gvsBinWeib$sims.list[[param]])
  u <- mean(gvsBinWeib$sims.list[[beta]][gvsBinWeib$sims.list[[param]] == 1], na.rm = TRUE)
  sd <- sd(gvsBinWeib$sims.list[[beta]][gvsBinWeib$sims.list[[param]] == 1], na.rm = TRUE)
  return(c(p, u, sd))
})

# Kuo & Mallick variable selection
modfile <- 'jags/binWeib/fulllog_truncated_censored_both_KM'
file.create(modfile)
inits <- write_jags_binWeib_varEval(
  filename = modfile,
  write = FALSE,
  jagsDat = dat,
  predictors = full,
  bin_response = 'solar',
  cont_response = 'ttd',
  censor = 'd',
  trunc = 'l',
  random = 'statei',
  out = outBinWeib)

params <- names(inits())
params <- params[grepl('beta', params)]%>%
  append(c(paste('w', full, sep = '_')))%>%
  append(paste('alpha[', 1:length(unique(weib_dat$statei)), ']', sep = ""))%>%
  append(paste('shape[', 1:length(unique(weib_dat$statei)), ']', sep = ""))%>%
  append(c('alpha_mu', 'shape_mu', 'alpha_v', 'shape_v'))

kmBinWeib <- jags(data = dat,
                   inits = inits,
                   parameters.to.save = params,
                   model.file = modfile,
                   n.chains = MCMCparams$nchain,
                   n.adapt = MCMCparams$nadapt,
                   n.iter = MCMCparams$niter,
                   n.burnin = MCMCparams$nburn,
                   n.thin = MCMCparams$nthin)

saveRDS(kmBinWeib, file = 'data/output/binweib/fulllog_truncated_censored_both_KM.rds')
kmProbs <- map_dfc(full, function(var){
  param <- paste('w', var, sep = "_")
  beta <- paste('beta', var, sep='_')
  p <- sum(kmBinWeib$sims.list[[param]] == 1)/length(kmBinWeib$sims.list[[param]])
  u <- mean(kmBinWeib$sims.list[[beta]][kmBinWeib$sims.list[[param]] == 1], na.rm = TRUE)
  sd <- sd(kmBinWeib$sims.list[[beta]][kmBinWeib$sims.list[[param]] == 1], na.rm = TRUE)
  return(c(p, u, sd))
})
# RESIDUALS
predictors <- c('open16', 'cultivated16', 'ssurgo', 'slope', 'GAP_Sts', 'road_dist', 'pdensity', 'income', 'lat')
pred_df <- as.data.frame(weib_dat[predictors])
betas <- as.data.frame(varProbs[2, predictors])
simRates <- t(as.matrix(pred_df) %*% t(as.matrix(betas)))
parameters <- lapply(1:nrow(pred_df), function(x){
  i <- weib_dat$statei[x]
  rate <- exp(simRates[,x] + log(kmBinWeib$mean$alpha[i]))
  shape <- kmBinWeib$mean$shape[i]
  return(data.frame(shape = shape, rate = rate))
})%>%bind_rows()

ttd_exp <- compute_weibull_response(parameters$shape, parameters$rate)-1
ttd_var <- compute_weibull_variance(parameters$shape, parameters$rate)
ttd_sim <- purrr::map2_dbl(parameters$shape, parameters$rate, function(x,y){
  rweibull(1, x, jagsWeib2Rweib(x, y))
}
)

# to calculated standardized residuals, we divide each residual by the sqrt of the variance
# of the weibull distribution

resid <- data.frame(obs=weib_dat$ttd,
                    exp = ttd_exp,
                    var = ttd_var,
                    sim = ttd_sim)%>%
  mutate(resid = obs - exp,
         resid_sim = sim - exp,
         stresid = resid/sqrt(var),
         stresid_sim = resid-sim/sqrt(var))

plot(resid$obs, resid$stresid)
##SPATIAL
library(sf)
library(spatial)
arrays <- st_read('data/solar_analysis_data.geojson')
bounds <- st_bbox(arrays)
pts <- st_centroid(arrays[arrays$year < 2022,])

areas <- c('New York' = 122049149763, 'Delaware' = 5045925646, 'Maryland' = 25151100280, 'Virginia'= 102257717110, 'Pennsylvania' = 115884442321)

calc_nni <- function(state){
  pts <- st_centroid(arrays[arrays$year < 2022 & arrays$state == state,])
  pp <- map2_dfr(pts$index, pts$geometry, function(x,y){
    return(data.frame(index = x, x = y[1], y = y[2]))}
  )
  dm <- dist(pp[,c('x', 'y')])
  dmat <- as.matrix(dm)
  dmat[dmat==0] <- NA
  nn <- apply(dmat,1,function(x){min(x, na.rm = TRUE)})
  enn <- sqrt(areas[state]/length(nn))*0.5
  return(nn/enn)
}

nnis <- lapply(1:length(areas), function(x){
  state <- names(areas)[x]
  nns <- calc_nni(state)
  df <- data.frame(state = rep(state, length(nns)), nearest = nns)
  return(df)
}
)%>%bind_rows()
glm(data = nnis, formula = nearest~state-1, family = gaussian(link = 'log'))

pp <- map2_dfr(pts$index, pts$geometry, function(x,y){
  return(data.frame(index = x, x = y[1], y = y[2]))}
)
dm <- dist(pp[,c('x', 'y')])
dmat <- as.matrix(dm)
dmat[dmat==0] <- NA
nn <- apply(dmat,1,function(x){min(x, na.rm = TRUE)})


pp <- ppp(pp$x, pp$y, c(bounds[1], bounds[3]), c(bounds[2], bounds[4]))
kest <- Kest(pp, rmax = 50000)
gest <- Gest(pp, r = seq(0, 5000, 50))
