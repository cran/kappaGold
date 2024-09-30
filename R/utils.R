# mkuhn, 2024-09-11
# helper/utility functions


#' delete-1 jackknife estimator
#' 
#' Quick simple jackknife routine to estimate bias and standard error of an
#' estimator.
#' @param est estimator function
#' @param idx maximal index vector for data of estimator
#' @return list with jackknife information, bias and SE
#' @references https://de.wikipedia.org/wiki/Jackknife-Methode
victorinox <- function(est, idx) {
  nx <- length(idx)
  estd <- numeric(nx)
  
  for (i in seq_along(idx)) {
    estd[i] <- est(idx[-i])
  }#rof
  
  estd_m <- mean(estd)
  SE_E <- sqrt((nx-1) * mean((estd - estd_m)^2))
  Bias_E <- (nx-1) * (estd_m - est(idx))
  
  list(bias_j = Bias_E, se_j = SE_E)
}
