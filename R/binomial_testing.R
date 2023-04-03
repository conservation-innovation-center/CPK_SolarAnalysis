test_dat <- make_jags_data(
  dataframe = cleanDF[cleanDF$year > 2016,],
  random = c('statei'),
  continuous = test,
  categorical = c('tmax', 'statei'),
  response = c('ttd', 'd', 'l')
)

filename <- 'jags/binNegbin/test_dcensoring'
file.create(filename)

baseString <- "
  model{
    #Likelihood
    
    # Distribution Parameters for Random Efects on
    # 1. Negative Binomial probability p ~ dbeta(shape1, shape2)
    # https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
    p_shp1 <- (((1-p_mu)/p_v) - (1/p_mu))*pow(p_mu, 2)
    p_shp2 <- (((1-p_mu)/p_v) - (1/p_mu))*p_mu*(1-p_mu)
    # https://cchecastaldo.github.io/BayesianShortCourse/content/lectures/MomentMatching.pdf
    # psi_shp1 <- (pow(psi_mu, 2) - pow(psi_mu, 3) - (psi_mu*psi_v))/psi_v
    # psi_shp2 <- (psi_mu - (2*pow(psi_mu, 2)) + pow(psi_mu, 3) - psi_v + (psi_mu*psi_v))/psi_v
    
    # random intercept per state on binomial intercept and neg binomial intercept
    for (s in 1:S){
      p[s] ~ dbeta(p_shp1, p_shp2)
    }
    
    # for every observation
    for (i in 1:N){
      # model for suitability of a site to solar development, z[i]
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- logit(psi_alpha) + psi_beta_slope*slope[i] + psi_beta_tree_cover16*tree_cover16[i] + psi_beta_road_dist*road_dist[i]
      
      # time to detection is a negative binomial process ignoring censoring
      # contResponseString
      ttd[i] ~ dnegbin(p[statei[i]]*z[i]+0.00001, 1)T(l[i],)
      # contPredictorString
      
      # model for censoring observed arrays due to not seeing into the future
      
      # OPTION 1:
      # whether we see an array is a bernouli process determined by
      # theta is 0 if site will never be developed (i.e. z[i] = 0) 
      #  or will be developed but not detected yet (i.e. z[i] = 1, ttd[i] > Tmax[i])
      
      d[i] ~ dinterval(ttd[i], tmax[i])
      # d[i] ~ dbern(theta[i])
      # theta[i] <- (z[i] * step(ttd[i] - tmax[i])) + (1 - z[i])
      
      # OPTION 2:
      # use jags dinterval to model censoring...?
      
      # Expected data under current model
      ttd_sim[i] ~ dnegbin(p[statei[i]]*z[i]+0.00001, 1)T(l[i],)
      ttd_exp[i] <- (1-(p[statei[i]]*z[i]+0.00001))/(p[statei[i]]*z[i]+0.00001)
      chi2[i] <- pow(ttd[i] - ttd_exp[i], 2)/ttd_exp[i]
      chi2_sim[i] <- pow(ttd_sim[i] - ttd_exp[i], 2)/ttd_exp[i]
    }
    
    fit <- mean(chi2[])
    fit_sim <- mean(chi2_sim[])
    bpv <- step(fit_sim - fit)
    
    #Priors
    # we want the mean of the gamma dist on weibull shape to be 1 and variance 1000
    # to simulate gamma(0.0001, 0.0001) with no state effect
    p_mu ~ dunif(0, 1)
    p_v ~ dunif(0, 0.2)
    # psi_alpha_mu ~ dunif(0, 1) 
    # psi_alpha_v ~ dunif(0, 0.2)
    psi_alpha ~ dunif(0,1)
    psi_beta_slope ~ dnorm(0, 0.0001)
    psi_beta_tree_cover16 ~ dnorm(0, 0.0001)
    psi_beta_road_dist ~ dnorm(0, 0.0001)
  }
"

# write our model string to specified file
cat(baseString, file = filename)

# dynamically create inits function that will generate initial values for beta parameters for each predictor
inits <- function(){
  ttdst <- test_dat$tmax+1 # creating some fake times greater than tmax by 1
  ttdst[!is.na(test_dat$ttd)] <- NA # this strictly overwrites the value of unobserved nodes!!!
  zst <- rep(1, length(test_dat$ttd))
  ls <- list(
    z = zst,
    ttd = ttdst,
    # psi_mu = 0.5,
    # psi_v = 0.08,
    # p_alpha = 0.5,
    psi_beta_slope = rnorm(1,0,1),
    psi_beta_tree_cover16 = rnorm(1,0,1),
    psi_beta_road_dist = rnorm(1,0,1)
  )
  return(ls)
}

params <- names(inits())
params <- c('d', 'z', 'ttd', 'p', 'psi')


out <- jags(data = test_dat,
            inits = inits,
            parameters.to.save = params,
            model.file = filename,
            n.chains = testMCMparams$nchain,
            n.adapt = testMCMparams$nadapt,
            n.iter = testMCMparams$niter,
            n.burnin = testMCMparams$nburn,
            n.thin = testMCMparams$nthin)
