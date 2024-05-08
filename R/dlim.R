#' Fit DLIM
#' @description Fit distributed lag interaction model
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @seealso \link[dlim]{predict.dlim}
#' @seealso \link[dlim]{plot_cumulative}
#' @seealso \link[dlim]{plot_DLF}
#' @export
#' @import mgcv
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#' @importFrom stats na.pass
#' @param y vector of response values (class "\code{numeric}")
#' @param x matrix of exposure history (columns) for individuals (rows) (class "\code{matrix}")
#' @param modifiers vector of modifying values (class "\code{numeric}")
#' @param z matrix of covariates, not including the modifier (class "\code{matrix}")
#' @param df_m degrees of freedom for modifier basis. Cannot specify for linear modification (model_type = "linear") (class "\code{numeric}")
#' @param df_l degrees of freedom for exposure time basis (class "\code{numeric}")
#' @param penalize \code{TRUE} to penalize model (class "\code{logical}")
#' @param pen_fn if penalizing, can specify "ps" for penalized B-splines or "cr" for cubic regression splines with penalties on second derivatives
#' @param mod_args a list of additional arguments for the spline function (must be named by argument)
#' @param lag_args a list of additional arguments for the spline function (must be named by argument)
#' @param fit_fn specify "gam" to use the \code{gam} function for data sets that are not very large, and specify "bam" to use the \code{bam} function for data sets that are very large. Default will fit using \code{gam}. (class "\code{character}")
#' @param model_type "linear" for a DLIM with linear interaction, "quadratic" for a DLIM with quadratic interaction, "standard" for a DLIM with splines (class "\code{character}")
#' @param ID group identifier for random intercept, only supported for penalized models 
#' @param ... Other arguments to pass to model fitting function
#' @example inst/examples/ex_dlim.R
#' @return This function returns a list that is an object of class "\code{dlim}" with the following components
#' \item{cb}{cross-basis (class "\code{matrix}")}
#' \item{fit}{model object (class "\code{lm}", "\code{glm}", "\code{gam}")}
#' \item{modifiers}{modifying values (class "\code{numeric}")}
#' \item{call}{model call}

dlim <- function(y, x, modifiers, z=NULL, df_m=NULL, df_l, penalize=TRUE, pen_fn = "ps", mod_args=NULL, lag_args=NULL, fit_fn="gam", model_type="standard", ID=NULL, ...){

  if(model_type == "standard" & is.null(df_m)){
    stop("Please specify df_m argument, the number of basis functions for the modifier basis.")
  }
  
  #set up design matrix for covariates and/or modifiers
  modifiers <- matrix(modifiers, ncol=1)
  if(!is.null(z)){
    design1 <- data.frame(intercept=rep(1,length(modifiers)), modifiers,z)
  }else{
    design1 <- data.frame(intercept= rep(1,length(modifiers)),modifiers)
  }

  #cross-basis
  if(penalize){
    cb <- cross_basis(x=x,M=modifiers,argmod=list(fun=pen_fn,df=df_m,arg=mod_args),arglag=list(fun=pen_fn,df=df_l,arg=lag_args), model_type = model_type)
  }else{
    cb <- cross_basis(x=x,M=modifiers,argmod=list(fun="ns",df=df_m,arg=mod_args),arglag=list(fun="ns",df=df_l,arg=lag_args), model_type = model_type)
  }

  #fit model
  CB <- cb$cb
  Z <- model.matrix(~ 0+.,model.frame(~ ., design1, na.action=na.pass)) #handles factor covariates and missing values
  if(penalize){
    if(!is.null(ID)){
      ID <- as.factor(ID) #make sure ID is a factor
      model <- do.call(fit_fn,list(formula=y~0+CB+Z+s(ID, bs="re"), paraPen = list(CB = cb$Slist), ...))
    }else{
      model <- do.call(fit_fn,list(formula=y~0+CB+Z, paraPen = list(CB = cb$Slist), ...))
    }
  }else{
    #add stop error if ID is provided for un-penalized model
    model <- do.call(fit_fn,list(formula=y~0+CB+Z, ...))
  }

  results <- list("cb" = cb, "fit" = model, "modifiers" = modifiers, call = match.call())

  class(results) <- "dlim"

  attr(results, "model_type") <- model_type

  return(results)

}
