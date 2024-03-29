#' Print DLIM Information
#' @description prints information about an object of class \code{dlim}
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @export
#' @param x a \code{dlim} object
#' @param ... further arguments passed to or from other methods
#' @return This function returns information about an object of class \code{dlim}

print.dlim <- function(x, ...){

  #check that object is class dlim
  if(!inherits(x,"dlim")){
    stop("Object not of class dlim")
  }

  cat("Object of class dlim", "\n")

  print(x$fit$family)

  print(x$call)

  cat("Modifier basis degrees of freedom:", x$cb$df_m,
      "\n")

  cat("Exposure time basis degrees of freedom:", x$cb$df_l,
      "\n \n")

  cat("Number of exposure time points:", x$cb$L+1,
      "\n \n")

  cat("Penalization:", ifelse(is.null(x$cb$Slist),"No", "Yes"),
      "\n \n")

  cat("n =", length(x$modifiers), "\n")


}
