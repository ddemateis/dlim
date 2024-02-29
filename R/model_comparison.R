#' Model Comparison
#' @description Compare models to test for interaction 
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @seealso \link[dlim]{dlim}
#' @export
#' @import dlnm
#' @import mgcv
#' @importFrom stats rpois
#' @importFrom stats logLik
#' @importFrom stats quantile
#' @param fit dlim object (must be fit with REML)
#' @param null "DLM", "linear" to specify the null model
#' @param x exposure
#' @param B number of bootstrap samples
#' @param conf.level The confidence level (class "\code{numeric}")
#' @return The function returns a decision to either reject or fail to reject the null model

model_comparison <- function(fit, null = "DLM", x, B, conf.level = 0.95){
  
  #make sure family is supported
  family <- fit$fit$family$family
  if(!(family=="gaussian" | family=="poisson")){
    stop("family not supported. Only Gaussian and Poisson currently supported")
  }
  
  #make sure method is supported
  method <- fit$fit$method
  if(method != "REML"){
    stop("method not supported. Only REML is currently supported")
  }
  
  #get info from dlim object
  df_l <- fit$cb$df_l
  df_m <- fit$cb$df_m
  method <- fit$fit$method
  model_type <- attr(fit, "model_type")
  y <- fit$fit$model$y
  modifiers <- fit$fit$model$Z[,2]
  if(dim(fit$fit$model$Z)[2]>2){
    z <- fit$fit$model$Z[,3:ncol(fit$fit$model$Z)]
  }else{
    z <- NULL
  }
  
  #Step 1: fit initial null model
  if(null=="DLM"){
    cb_dlm <- crossbasis(x=x,
                         argvar=list(fun="lin"),
                         arglag = list(fun="ps",
                                       df=df_l))
    penalty <- cbPen(cb_dlm)
    covariates <- as.data.frame(cbind(z, modifiers))
    design2 <- model.matrix(~ .,
                            model.frame(~ ., covariates, na.action=na.pass))
    initial_null <- gam(y~cb_dlm+design2,
                     paraPen = list(cb_dlm = penalty), 
                     method = method,
                     family = family)
  }else{
    null_model <- dlim(y = y, 
                       x = x,
                       z = z,
                       modifiers = modifiers,
                       df_m = df_m,
                       df_l = df_l,
                       penalize = TRUE,
                       method = method,
                       model_type = null,
                       family = family)
    initial_null <- null_model$fit
  }
  
  LLR <- c()
  for(i in 1:B){
    #Step 2: generate bootstrap sample using null model
    if(family == "gaussian"){
      res_sample <- sample(initial_null$residuals, 
                           size = length(initial_null$residuals),
                           replace = TRUE)
      boot_y <- initial_null$fitted.values + res_sample
    }else if(family == "poisson"){
      boot_y <- rpois(length(y), initial_null$fitted.values)
    }
    
    #Step 3: Calculate bootstrap likelihood ratio using sample in step 2
    if(null=="DLM"){
      null_model <- gam(boot_y~cb_dlm+design2,
                        paraPen = list(cb_dlm = penalty), 
                        method = method)
    }else{
      null_model <- dlim(y = boot_y, 
                         x = x,
                         z = z,
                         modifiers = modifiers,
                         df_m = df_m,
                         df_l = df_l,
                         penalize = TRUE,
                         method = method,
                         model_type = null)$fit
    }
    
    full_model <- dlim(y = boot_y, 
                       x = x,
                       z = z,
                       modifiers = modifiers,
                       df_m = df_m,
                       df_l = df_l,
                       penalize = TRUE,
                       method = method,
                       model_type = model_type,
                       family = family)$fit
    LLR[i] <- as.numeric(logLik(full_model) - logLik(null_model))
  }
  
  
  #compute observed test statistic
  obs_LLR <- logLik(fit$fit) - logLik(initial_null)
  crit_LLR <- quantile(LLR, conf.level)
  
  return(ifelse(obs_LLR > crit_LLR, "reject", "FTR"))
}
