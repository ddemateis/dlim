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
#' @param MaME True to include modifier as linear main effect
#' @return This function returns an object of class dlim
#' \item{cb}{cross-basis (matrix)}
#' \item{fit}{model object (gam)}
#' \item{modifiers}{modifying values (vector)}
#' \item{call}{model call}

dlim <- function(y, x, modifiers, z=NULL, df_m, df_l, penalize=T, fit_fn="gam", MaME=T, new_approach = NA, model_type=4, penalty_same=F, method="GCV.Cp", family="gaussian"){

  #set up design matrix for covariates and/or modifiers
  modifiers <- matrix(modifiers, ncol=1)
  if(!is.null(z)){
    if(MaME){
      design1 <- data.frame(intercept=rep(1,length(modifiers)), modifiers,z)
    }else{
      design1 <- data.frame(intercept=rep(1,length(modifiers)),z)
    }
  }else{
    if(MaME){
      design1 <- data.frame(intercept= rep(1,length(modifiers)),modifiers)
    }else{
      design1 <- data.frame(intercept=rep(1,length(modifiers)))
    }
  }

  #cross-basis
  if(penalize){
    cb <- cross_basis(x=x,M=modifiers,argmod=list(fun="ps",df=df_m),arglag=list(fun="ps",df=df_l), model_type = model_type, penalty_same = penalty_same)
  }else{
    cb <- cross_basis(x=x,M=modifiers,argmod=list(fun="ns",df=df_m),arglag=list(fun="ns",df=df_l), model_type = model_type, penalty_same = penalty_same)
  }

  #fit model
  CB <- cb$cb
  Z <- model.matrix(~ 0+.,model.frame(~ ., design1, na.action=na.pass)) #handles factor covariates and missing values
  if(penalize){
    if(is.na(new_approach)){
      model <- do.call(fit_fn,list(formula=y~0+CB+Z, paraPen = list(CB = cb$Slist), method=method))
    }else if(new_approach=="regress_out"){
      model1 <- lm(y~modifiers)
      resid <- model1$residuals
      model <- do.call(fit_fn,list(formula=resid~0+CB+Z[,colnames(Z)!="modifiers"], paraPen = list(CB = cb$Slist), method = "REML"))
    }else if(new_approach=="mod_spline"){
      m_basis <- bs(Z[,colnames(Z)=="modifiers"], df=3) #return this
      model <- do.call(fit_fn,list(formula=y~0+CB+ m_basis + Z[,colnames(Z)!="modifiers"],paraPen = list(CB = cb$Slist), method = "REML"))
    }else if(new_approach=="ridge_penalty"){
      modifier_design <- Z[,colnames(Z)=="modifiers"]
      covariates_design <- Z[,colnames(Z)!="modifiers"]
      model <- do.call(fit_fn,list(formula=y~0+CB+modifier_design+covariates_design, paraPen = list(CB = cb$Slist, modifier_design = list(rank = 1, diag(1))), method = "REML"))
    }else if(new_approach=="ridge_penalty_fixed"){
      modifier_design <- Z[,colnames(Z)=="modifiers"]
      covariates_design <- Z[,colnames(Z)!="modifiers"]
      model <- do.call(fit_fn,list(formula=y~0+CB+modifier_design+covariates_design, paraPen = list(CB = cb$Slist, modifier_design = list(sp=100, rank = 1, diag(1))), method = "REML"))
    }
  }else{
    model <- do.call(fit_fn,list(formula=y~0+CB+Z, family=family))
  }

  results <- list("cb" = cb, "fit" = model, "modifiers" = modifiers, call = match.call)

  class(results) <- "dlim"

  return(results)

}
