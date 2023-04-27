#' DLIM Predictions
#' @description Predicted values based on a \code{dlim} object.
#' @seealso \link[dlim]{dlim}
#' @export
#' @import splines
#' @import dlnm
#' @param object an object of class dlim
#' @param newdata a vector of new modifier values for prediction (default is NULL)
#' @param type Type of prediction ("DLF" for the estimated distributed lag functions, "CE" for cumulative effects, "response" for fitted values, or any combination of these in a vector)
#' @param conf.level The confidence level
#' @return This function returns a list of 3 elements:
#' \item{est_dlim}{cumulative and/or point-wise estimates, standard errors, and confidence intervals (list)}
#' \item{cb}{cross-basis object ()}
#' \item{model}{model object (gam)}

predict.dlim <- function(object, newdata=NULL, type=c("DLF","CE", "response"), conf.level = 0.95){

  if(class(object)!="dlim"){
    stop("Object not of class dlim")
  }

  alpha <- 1 - conf.level

  cb <- object$cb
  fit <- object$fit
  cb_dlm <- object$cb_dlm
  model_dlm <- object$model_dlm
  if(is.null(newdata)){
    modifiers <- object$modifiers
  }else{
    modifiers <- newdata
  }

  m <- length(modifiers)

  est_dlim <- list()

  #reconstruct B_mod for given modifiers
  if(attr(object,"model_type")==4){
    if(class(cb$B_mod)[1]=="ps"){
      B_mod <- ps(modifiers, knots=attr(cb$B_mod,"knots"),intercept = T)#mxdf_m
    }else{
      B_mod <- predict(cb$B_mod,modifiers) #mxdf_m
    }
  }else if(attr(object,"model_type")==2){
    B_mod <- cbind(rep(1,m), modifiers)
  }else if(attr(object,"model_type")==3){
    B_mod <- cbind(rep(1,m), modifiers, modifiers^2)
  }


  idx <- grep("CB",names(fit$coefficients))#includes only cross-basis elements
  coef <- fit$coefficients[idx]

  #estimate cumulative effects
  if("CE" %in% type){

    cb_est <- matrix(NA, m, cb$df_l*cb$df_m) #initialize w* matrix
    for(i in 1:m){
      B_i <- matrix(rep(1,cb$L+1),ncol=1)%x%matrix(B_mod[i,],nrow=1)  #L+1 x df_m
      B_temp <- t(B_i) %*% cb$B_lag  # df_m x df_l
      cb_est[i,] <- c(t(B_temp)) #1 x df_m*df_l, vectorizes by column, so this is grouped by modifier
    }

    #Cumulative effects
    betas_cumul <- cb_est %*% coef #eq. 11
    rownames(betas_cumul) <- modifiers
    colnames(betas_cumul) <- "cumululative"

    #cov_mat <- diag(cb_est%*%vcov(model)[idx,idx]%*%t(cb_est)) #standard covariance
    cov_mat <- apply(cb_est, 1, function(w) tcrossprod(w, vcov(fit)[idx,idx]) %*% w) #eq. 12
    covs_cumul <- sqrt(cov_mat)
    cumul_UB <- betas_cumul + qnorm(1-alpha/2)*covs_cumul
    cumul_LB <- betas_cumul - qnorm(1-alpha/2)*covs_cumul

    est_dlim$betas_cumul <- betas_cumul
    est_dlim$cumul_LB <- cumul_LB
    est_dlim$cumul_UB <- cumul_UB
    est_dlim$cumul_SE <- covs_cumul
  }



  #DLFs
  if("DLF" %in% type){
    betas <- c()
    covs <- c()

    betas <- t(apply(B_mod, 1, dlf_betas, cb$B_lag, coef)) #eq. 7
    covs <- t(apply(B_mod, 1, dlf_covs, cb$B_lag, fit, idx)) #eq. 9
    UB <- betas + qnorm(1-alpha/2)*covs
    LB <- betas - qnorm(1-alpha/2)*covs

    est_dlim$betas <- betas
    est_dlim$LB <- LB
    est_dlim$UB <- UB
    est_dlim$SE <- covs
  }

  if("response" %in% type){
    print("Returning fitted values is in progress.")
  }

  est_dlim$modifiers <- modifiers
  result <- c()
  result$est_dlim <- est_dlim
  result$cb <- object$cb
  result$fit <- object$fit

  return(result)
}
