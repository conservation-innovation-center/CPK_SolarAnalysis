#' Standardize continuous data to z-score
#' 
#' @param x numeric vector
#' @return vector of same length as \code{x} containing z-transformed values
#' @examples 
#' v <- runif(100, 10, 12)
#' calc_z(v)
calc_z <- function(x){
  mn <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  z <- (x-mn)/s
  return(z)
}

#' Standardize continuous data to (0,1)
#' 
#' @param x numeric vector
#' @return vector of same length as \code{x} containing standardized values
#' @examples 
#' v <- runif(100, 10, 12)
#' standardize(v)
standardize <- function(x){
  mx <-max(x, na.rm = TRUE)
  mn <-min(x, na.rm = TRUE)
  z <- (x-mn)/(mx-mn)
  return(z)
}

logistic <- function(x){
  y <- exp(x)/(exp(x) + 1)
  return(y)
}

plot_post_hist <- function(jags_out, param){
  # get index of desired parameter in jags simulation object
  index = grep(param, names(jags_out$sims.list))
  hist(out$sims.list[[index]], breaks = 50, col = 'grey', freq = F, xlab = param)
  abline(v = quantile(out$sims.list[[index]], prob=c(0.025, 0.975)), col = 'red', lwd = 2)
}

make_jags_model <- function(class, random, continuous){
  if class == 'binomial'{
    
  }
}