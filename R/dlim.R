#' Fit DLIM
#' @description Fit distributed lag interaction model
#' @export
#' @import mgcv
#' @import dlnm
#' @import splines
#' @param y vector of response values
#' @param x matrix of exposure history
#' @param modifiers vector of modifying values
#' @param z matrix of z
#' @param df_m degrees of freedom for modifier basis
#' @param df_l degrees of freedom for exposure time basis
#' @param penalize True to penalize model
#' @param fit_fn specify "gam" to use the \code{gam} function for data sets that are not very large, and specify "bam" to use the \code{bam} function for data sets that are very large. Default will fit using \code{gam}.
#' @param ... Other arguments to pass to model fitting function
#' @return This function returns an object of class dlim
#' \item{cb}{cross-basis (matrix)}
#' \item{fit}{model object (gam)}
#' \item{modifiers}{modifying values (vector)}
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
