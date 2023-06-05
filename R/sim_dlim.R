#' Fit DLIM for simulation
#' @description Fit DLIM for simulation
#' @seealso \link[dlim]{dlim}
#' @export
#' @import mgcv
#' @import dlnm
#' @param data output from sim_data()
#' @param df_m degrees of freedom for modifiers
#' @param df_l degrees of freedom for lags
#' @param penalize True to penalize model
#' @param fit_dlm True to additionally fit dlm for comparison
#' @param model_type "linear" for a DLIM with linear interaction, "quadratic" for a DLIM with quadratic interaction, "standard" for a DLIM with splines
#' @param ... arguments to pass to model fitting function
#' @return This function returns an object of class dlim

sim_dlim <- function(data, df_m, df_l, penalize=T, fit_dlm=F, model_type="standard",...){


  model <- dlim(y = data$y,
                x = data$x,
                modifiers = data$modifiers,
                z = data$Z,
                df_m = df_m,
                df_l = df_l,
                penalize=penalize,
                model_type=model_type,
                ...)

  if(fit_dlm){
    #set up
    modifiers <- matrix(data$modifiers, ncol=1)
    if(!is.null(data$Z)){
      design1 <- data.frame(intercept=rep(1,length(data$modifiers)), data$modifiers,data$Z)
    }else{
      design1 <- data.frame(intercept=rep(1,length(data$modifiers)), data$modifiers)
    }
    y <- data$y

    #cross-basis
    if(penalize){
      cb_dlm <- crossbasis(x=data$x,argvar=list(fun="lin"),arglag = list(fun="ps",df=df_l))
    }else{
      cb_dlm <- crossbasis(x=data$x,argvar=list(fun="lin"),arglag = list(fun="ns",df=df_l))
    }

    #fit model
    Z <- model.matrix(~ 0+.,model.frame(~ ., design1, na.action=na.pass)) #handles factor covariates and missing values
    if(penalize){
      penalty <- cbPen(cb_dlm)
      model_dlm <- do.call("gam",list(formula=0+y~cb_dlm+Z,paraPen = list(cb_dlm = penalty), method = "REML"))
    }else{
      model_dlm <- do.call("gam",list(formula=y~0+cb_dlm+Z))
    }

    # if(penalize){
    #   penalty <- cbPen(cb_dlm)
    #   if(!is.null(data$Z)){
    #     model_dlm <- gam(data$y~cb_dlm+data$Z,paraPen = list(cb_dlm = penalty), method = "REML")
    #   }else{
    #     model_dlm <- gam(data$y~cb_dlm,paraPen = list(cb_dlm = penalty), method = "REML")
    #   }
    # }else{
    #   if(!is.null(data$Z)){
    #     model_dlm <- gam(data$y~cb_dlm+data$Z)
    #   }else{
    #     model_dlm <- gam(data$y~cb_dlm)
    #   }
    # }


  }

  if(fit_dlm){
    results <- list("cb" = model$cb, "fit" = model$fit, "cb_dlm" = cb_dlm, "model_dlm" = model_dlm, "true_betas" = data$betas, "modifiers" = data$modifiers, data = data)
  }else{
    results <- list("cb" = model$cb, "fit" = model$fit, "true_betas" = data$betas, "modifiers" = data$modifiers, data = data)
  }

  attr(results, "model_type") <- model_type
  class(results) <- "sim_dlim"

  return(results)

}
