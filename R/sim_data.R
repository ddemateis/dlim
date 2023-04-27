#' Simulate Data
#' @description Simulate data to use with the \pkg{dlim} package. There are different scenarios to choose for simulation.
#' @export
#' @import tsModel
#' @param x a numeric time series vector of length n or matrix of lagged exposures for n individuals
#' @param L a numeric vector of length 1 containing the number of lag terms. This is required if \code{x} is vector, and is not used if \code{x} is a matrix.
#' @param modifiers vector of length n containing modifying values
#' @param noise a vector of length 1 containing the standard deviation for a normal distribution with mean 0 used to add noise to the simulated response values.
#' @param type a vector containing the number 1, 2, or 3 for simulation modification type: none, linear, complex. Default is 2
#' @param SNR value to magnify the signal of the true distributed lag functions, default is no magnification
#' @param ncovariates number of covariates to add to the model, numeric vector of length 1.
#' @return This returns a list of 9 items:
#' \item{x}{a lagged exposure matrix. If \code{x} was a matrix, it is unchanged.}
#' \item{L}{a numeric vector o flength 1 containing the number of lag terms}
#' \item{modifiers}{the \code{modifiers} argument}
#' \item{y}{a numeric vector of length \code{nrow(x)} containing the perturbed simulated response values. obtained by adding a \code{N(0,noise)} random draw to \code{y}}
#' \item{betas}{a matrix containing true coefficients for each lag/modifier combination, with each row representing a lag and each column a modifier}
#' \item{betas_cumul}{a numeric vector of length \code{L+1} containing cumululative true coefficients for the lag terms, summed over modifiers}
#' \item{sim_type}{simulation type 1, 2, or 3 (numeric)}
#' \item{Z}{covariates (matrix)}
#' \item{gammas}{true coefficients for the covariates (numeric)}



sim_data <- function(x, L=NULL, modifiers, noise=1, type=2, SNR, ncovariates=0, MaME=F, gamma=1){

  #create lagged structure
  if(is.vector(x)){
    X <- Lag(x,0:L)[-c(1:L),]
    modifiers <- modifiers[-c(1:L)]
  }else{
    L <- ncol(x)-1
    X <- x
    modifiers <- modifiers
  }

  #Create Betas
  betas <- sim_dlf(L,modifiers,type)
  betas_cumul <- colSums(betas)

  y_mean <- colSums(t(X)*betas)

  #if SNR is provided, but noise is not, reset noise to the SD based on data
  #if SNR is not provided, noise must be provided (default to noise of 1)
  if(!missing(SNR)){
    noise <- sd(y_mean)/SNR
  }

  #Create gammas and covariates
  if(MaME){##modifier as main effect
    if(ncovariates!=0){
      gammas <- c(gamma,matrix(rnorm(ncovariates),ncol=1))
      Z <- matrix(rnorm(length(modifiers)*(ncovariates)), ncol = ncovariates)
      mod_Z <- cbind(modifiers,Z)
      y <- y_mean + mod_Z%*%gammas + rnorm(length(y_mean),0,noise)
    }else{
      Z <- NULL
      gammas <- gamma
      y <- y_mean + matrix(M,ncol=1)*gammas + rnorm(length(y_mean),0,noise)
    }
  }else{##modifier not as main effect
    if(ncovariates!=0){
      gammas <- matrix(rnorm(ncovariates),ncol=1)
      Z <- matrix(rnorm(length(modifiers)*ncovariates), ncol = ncovariates)
      y <- y_mean + Z%*%gammas + rnorm(length(y_mean),0,noise)
    }else{
      Z <- NULL
      gammas <- NULL
      y <- y_mean + rnorm(length(y_mean),0,noise)
    }
  }

  result <- list(x=X,
                 L=L,
                 modifiers=modifiers,
                 y=y,
                 betas=betas,
                 betas_cumul=betas_cumul,
                 Z = Z,
                 gammas = gammas)

  return(result)
}
