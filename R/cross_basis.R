#' Build crossbasis
#' @description Creates cross-basis using natural splines for regression in DLIM
#' @seealso \link[dlim]{dlim}
#' @export
#' @import tsModel
#' @import splines
#' @import dlnm
#' @param x a numeric time series vector of length n or matrix of lagged exposures (columns) for n individuals (rows)
#' @param M vector of length n containing modifier values
#' @param L a numeric vector of length 1 containing the number of lag terms. This is required if \code{x} is vector, and is not used if \code{x} is a matrix.
#' @param argmod a list: $fun is the spline function for the modifier ("ps" or "cr" to penalize), $arg is a list of arguments for the spline function (must be named by argument), $df is the degrees of freedom, $sp is optional smoothing parameter
#' @param arglag a list: $fun is the spline function for the lag ("ps" or "cr" to penalize), $arg is a list of arguments for the spline function (must be named by argument), $df is the degrees of freedom, $sp is optional smoothing parameter
#' @param model_type "linear" for a DLIM with linear interaction, "quadratic" for a DLIM with quadratic interaction, "standard" for a DLIM with splines
#' @return This function returns a list of 5 or 6 elements:
#' \item{cb }{cross-basis (matrix)}
#' \item{B_lag}{lag basis (basis matrix)}
#' \item{B_mod}{modifier basis (basis matrix)}
#' \item{df_l}{lag degrees of freedom (numeric)}
#' \item{df_m}{modifier degrees of freedom (numeric)}
#' \item{L}{number of lags (numeric)}
#' \item{Slist}{lag and modifier penalty matrices, if penalizing (list)}
cross_basis <- function(x,M,L=NULL,argmod=list(),arglag=list(), model_type="standard"){
  #set up
  if(is.vector(x)){
    X <- Lag(x,0:L)[-c(1:L),]
    M <- M[-c(1:L)]
  }else{
    L <- ncol(x)-1
    X <- x
  }
  n <- length(M)
  df_l <- arglag$df
  df_m <- argmod$df

  #Create bases
  B_lag <- do.call(arglag$fun,list(0:L, df=df_l, intercept = T, arglag$arg)) #L+1xdf_l
  if(model_type=="linear"){
    B_mod <- cbind(rep(1,n), M)
    df_m <- 2 #number of col in B_mod
  }else if(model_type=="quadratic"){
    B_mod <- cbind(rep(1,n), M, M^2)
    df_m <- 3 #number of col in B_mod
  }else if(model_type=="standard"){
    B_mod <- do.call(argmod$fun,list(M,df=df_m,intercept = T)) #nxdf_m
  }


  #Cross the modified exposures with the lag basis.
  # create_cb <- function(i){
  #   B_i <- B_mod[i,] #eqn 2
  #   X_i <- matrix(c(as.matrix(X[i,])),ncol=1)
  #   X_tilde <- X_i%*%B_i #L+1 x df_m
  #   B_temp <- t(X_tilde) %*% B_lag  # df_m x df_l
  #   return(c(t(B_temp))) #1 x df_m*df_l, vectorizes by column, so this is grouped by modifier
  # }
  # cb <- t(sapply(1:n, create_cb))
  
  #same as above, but with a for loop not apply function
  cb <- matrix(NA, n, df_l*df_m)
  for(i in 1:n){
    B_i <- B_mod[i,] #eqn 2
    X_i <- matrix(c(as.matrix(X[i,])),ncol=1)
    X_tilde <- X_i%*%B_i #L+1 x df_m
    B_temp <- t(X_tilde) %*% B_lag  # df_m x df_l
    cb[i,] <- c(t(B_temp)) #1 x df_m*df_l, vectorizes by column, so this is grouped by modifier
  }

  #create penalty argument for model if penalizing
  if((argmod$fun=="ps" | argmod$fun=="cr") & (arglag$fun=="ps" | argmod$fun=="cr")){

    #the S attribute is t(D_2) %*% D_2, where D_2 is second order difference matrix. See Gasparini 2017

    if(model_type=="standard"){
      Slag_tmp <- diag(df_m) %x% attr(B_lag,"S") #S* for exposure time
      Slag <- Slag_tmp/eigen(Slag_tmp, symmetric = TRUE,only.values = TRUE)$values[1] #penalty for exposure time basis
      Smod_tmp <- attr(B_mod,"S") %x% diag(df_l) #S* for modifier
      Smod <- Smod_tmp/eigen(Smod_tmp, symmetric = TRUE,only.values = TRUE)$values[1] #penalty for modifier basis
    }else if(model_type=="linear" | model_type=="quadratic"){
      if(model_type=="linear"){
        Slag_tmp1 <- diag(c(1,0)) %x% attr(B_lag,"S")
        Slag1 <- Slag_tmp1/eigen(Slag_tmp1, symmetric = TRUE,only.values = TRUE)$values[1]
        Slag_tmp2 <- diag(c(0,1)) %x% attr(B_lag,"S")
        Slag2 <- Slag_tmp2/eigen(Slag_tmp2, symmetric = TRUE,only.values = TRUE)$values[1]
      }else{
        Slag_tmp1 <- diag(c(1,0,0)) %x% attr(B_lag,"S")
        Slag1 <- Slag_tmp1/eigen(Slag_tmp1, symmetric = TRUE,only.values = TRUE)$values[1]
        Slag_tmp2 <- diag(c(0,1,0)) %x% attr(B_lag,"S")
        Slag2 <- Slag_tmp2/eigen(Slag_tmp2, symmetric = TRUE,only.values = TRUE)$values[1]
        Slag_tmp3 <- diag(c(0,0,1)) %x% attr(B_lag,"S")
        Slag3 <- Slag_tmp3/eigen(Slag_tmp3, symmetric = TRUE,only.values = TRUE)$values[1]
      }
    }

    #save Slist for paraPen
    if(is.null(arglag$sp)&is.null(argmod$sp)){ #without specifying smoothing parameters
      if(model_type=="standard"){
        Slist <- list(Smod,Slag)
        names(Slist) <- c("Smod","Slag")
      }else if(model_type=="linear" | model_type=="quadratic"){
        if(model_type=="linear"){
          Slist <- list(Slag1, Slag2)
          names(Slist) <- c("Slag1", "Slag2")
        }else if(model_type=="quadratic"){
          Slist <- list(Slag1, Slag2, Slag3)
          names(Slist) <- c("Slag1", "Slag2", "Slag3")
        }
      }
    }else if(!is.null(arglag$sp)&!is.null(argmod$sp)){ #with specifying smoothing parameters
      if(model_type=="linear" | model_type=="quadratic"){
        Slist <- list(Slag,arglag$sp)
        names(Slist) <- c("Slag","sp")
        attr(Slist$sp,"names") <- "Slag"
      }else if(model_type=="standard"){
        Slist <- list(Smod,Slag,c(argmod$sp,arglag$sp))
        names(Slist) <- c("Smod","Slag","sp")
        attr(Slist$sp,"names") <- c("Svar", "Slag")
      }
    }else{
      stop("Only set up to take smoothing parameters for both mod and lag")
    }
    result <- list(cb=cb,B_lag=B_lag,B_mod=B_mod,df_l=df_l,df_m=df_m,L=L,Slist=Slist)
  }else{
    result <- list(cb=cb,B_lag=B_lag,B_mod=B_mod,df_l=df_l,df_m=df_m,L=L)
  }

  class(result) <- "cross_basis"

  return(result)
}
