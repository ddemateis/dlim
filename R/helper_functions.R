#' @description This is a helper function in predict.dlim()
#' @param B_mod_i i^{th} row of modifier basis
#' @param B_lag lag basis
#' @param coef model coefficients
#' @return This function returns estimated point-wise effects

dlf_betas <- function(B_mod_i, B_lag, coef){
  K <- (matrix(B_mod_i,nrow=1)%x%B_lag)
  betas <- t(K%*%coef) #eq. 2.4.1
  return(betas)
}

#' @description This is a helper function in predict.dlim()
#' @param B_mod_i i^{th} row of modifier basis
#' @param B_lag lag basis
#' @param fit fitted model object
#' @param idx index of estimates to include
#' @return This function returns estimated variance of point-wise effects

dlf_covs <- function(B_mod_i, B_lag, fit, idx){
  K <- (matrix(B_mod_i,nrow=1)%x%B_lag)
  cov_mat <- K%*%vcov(fit)[idx,idx]%*%t(K) #standard covariance
  covs <- sqrt(diag(cov_mat))
  return(covs)
}
