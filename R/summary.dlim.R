#' Summarizing DLIM
#' @description prints summary of object of class \code{dlim}
#' @export
#' @param x a \code{dlim} object
#' @return This function returns a summary for an object of class \code{dlim}

summary.dlim <- function(x){

  #check that object is class dlim
  if(!inherits(x,"dlim")){
    stop("Object not of class dlim")
  }

  print("For predictions, use the predict.dlim() function to see estimates and predictions of the distributed lag function. See help(predict.dlim) for more information.")

}
