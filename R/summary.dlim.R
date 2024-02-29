#' Summarizing DLIM
#' @description prints summary of object of class \code{dlim}
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @export
#' @param object a \code{dlim} object
#' @param ... additional arguments affecting the summary produced
#' @return This function returns a summary for an object of class \code{dlim}

summary.dlim <- function(object, ...){

  #check that object is class dlim
  if(!inherits(object,"dlim")){
    stop("Object not of class dlim")
  }

  print("For predictions, use the predict.dlim() function to see estimates and predictions of the distributed lag function. See help(predict.dlim) for more information.")

}
