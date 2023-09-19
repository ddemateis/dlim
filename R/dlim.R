#' Fit DLIM
#' @description Fit distributed lag interaction model
#' @export
#' @importFrom mgcv
#' @importFrom dlnm
#' @importFrom splines
#' @param y vector of response values (class "\code{numeric}")
#' @param x matrix of exposure history (columns) for individuals (rows) (class "\code{matrix}")
#' @param modifiers vector of modifying values (class "\code{numeric}")
#' @param z matrix of covariates, not including the modifier (class "\code{matrix}")
#' @param df_m degrees of freedom for modifier basis (class "\code{numeric}")
#' @param df_l degrees of freedom for exposure time basis (class "\code{numeric}")
#' @param penalize \code{TRUE} to penalize model (class "\code{logical}")
#' @param fit_fn specify "gam" to use the \code{gam} function for data sets that are not very large, and specify "bam" to use the \code{bam} function for data sets that are very large. Default will fit using \code{gam}. (class "\code{character}")
#' @param model_type "linear" for a DLIM with linear interaction, "quadratic" for a DLIM with quadratic interaction, "standard" for a DLIM with splines (class "\code{character}")
#' @param ... Other arguments to pass to model fitting function
#' @return This function returns a list that is an object of class "\code{dlim}" with the following components
#' \item{cb}{cross-basis (class "\code{matrix}")}
#' \item{fit}{model object (class "\code{lm}", "\code{glm}", "\code{gam}")}
#' \item{modifiers}{modifying values (class "\code{numeric}")}
#' \item{call}{model call}

dlim <- function(y, x, modifiers, z=NULL, df_m, df_l, penalize=T, fit_fn="gam", model_type="standard", ...){

  #set up design matrix for covariates and/or modifiers
  modifiers <- matrix(modifiers, ncol=1)
  if(!is.null(z)){
    design1 <- data.frame(intercept=rep(1,length(modifiers)), modifiers,z)
  }else{
    design1 <- data.frame(intercept= rep(1,length(modifiers)),modifiers)
  }

  #cross-basis
  if(penalize){
    cb <- cross_basis(x=x,M=modifiers,argmod=list(fun="ps",df=df_m),arglag=list(fun="ps",df=df_l), model_type = model_type)
  }else{
    cb <- cross_basis(x=x,M=modifiers,argmod=list(fun="ns",df=df_m),arglag=list(fun="ns",df=df_l), model_type = model_type)
  }

  #fit model
  CB <- cb$cb
  Z <- model.matrix(~ 0+.,model.frame(~ ., design1, na.action=na.pass)) #handles factor covariates and missing values
  if(penalize){
    model <- do.call(fit_fn,list(formula=y~0+CB+Z, paraPen = list(CB = cb$Slist), ...))
  }else{
    model <- do.call(fit_fn,list(formula=y~0+CB+Z, ...))
  }

  results <- list("cb" = cb, "fit" = model, "modifiers" = modifiers, call = match.call)

  class(results) <- "dlim"

  attr(results, "model_type") <- model_type

  return(results)

}
