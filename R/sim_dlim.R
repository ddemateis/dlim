#' Fit DLIM for simulation
#' @description Fit DLIM for simulation
#' @seealso \link[dlim]{dlim}
#' @seealso \link[dlim]{sim_data}
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @export
#' @import mgcv
#' @import dlnm 
#' @param data output from \code{sim_data}
#' @param df_m degrees of freedom for modifiers
#' @param df_l degrees of freedom for lags
#' @param penalize True to penalize model
#' @param pen_fn if penalizing, can specify "ps" for penalized B-splines or "cr" for cubic regression splines with penalties on second derivatives
#' @param mod_args a list of additional arguments for the spline function (must be named by argument)
#' @param lag_args a list of additional arguments for the spline function (must be named by argument)
#' @param fit_dlm True to additionally fit dlm for comparison
#' @param model_type "linear" for a DLIM with linear interaction (linear modifier basis), "quadratic" for a DLIM with quadratic interaction (quadratic modifier basis), "nonlinear" for a DLIM with non-linear interaction (spline modifier basis)
#' @param ... arguments to pass to model fitting function
#' @return This function returns an object of class "\code{sim_dlim}"
#' \item{cb}{DLIM cross-basis (class "\code{cross-basis}")}
#' \item{fit}{DLIM model fit (class "\code{lm}", "\code{glm}", "\code{gam}")}
#' \item{cb_dlm}{DLM cross-basis (class "\code{crossbasis}")}
#' \item{model_dlm}{DLM model fit (class "\code{lm}", "\code{glm}", "\code{gam}")}
#' \item{true_betas}{true linear effect of the exposure on the response for each individual and time point (class "\code{matrix}")}
#' \item{modifiers}{\code{modifiers} from \code{numeric}}
#' \item{data}{\code{data} (class "\code{list}")}

sim_dlim <- function(data, df_m, df_l, penalize=TRUE, pen_fn = "ps", mod_args=NULL, lag_args=NULL, fit_dlm=FALSE, model_type="nonlinear",...){
  
  if(model_type == "standard"){
    lifecycle::deprecate_warn("0.3.0", "model_type = 'standard'", "model_type = 'nonlinear'")
    model_type <- "nonlinear"
  }
  
  #fit DLIM
  model <- dlim(y = data$y,
                x = data$x,
                modifiers = data$modifiers,
                z = data$Z,
                df_m = df_m,
                df_l = df_l,
                penalize=penalize,
                model_type=model_type,
                mod_args = mod_args,
                lag_args = lag_args,
                ...)

  #fit DLM
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
      cb_dlm <- crossbasis(x=data$x,argvar=list(fun="lin"),arglag = list(fun=pen_fn,df=df_l))
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
