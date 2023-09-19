#' Simulate Data
#' @description Simulate data to use with the \pkg{dlim} package. There are different effect modification scenarios to choose for simulation.
#' @export
#' @importFrom tsModel
#' @param x a time series vector of length \code{n} or matrix of lagged exposures for \code{n} individuals (class "\code{numeric}", "\code{matrix}")
#' @param L a vector of length 1 containing the number of lag terms. This is required if \code{x} is vector, and is not used if \code{x} is a matrix (class "\code{numeric}")
#' @param modifiers vector of length \code{n} containing modifying values (class "\code{numeric}")
#' @param noise a vector of length 1 containing the standard deviation for a normal distribution with mean 0 used to add noise to the simulated response values. Must proivde if \code{SNR} is not provided (class "\code{numeric}")
#' @param type a vector containing the number 1, 2, 3, or 4 for simulation modification type: none, linear, non-linear shift, non-linear shift with linear scale (class "\code{numeric}")
#' @param SNR The signal-to-noise ratio. If \code{SNR} is provided, but \code{noise} is not, \code{noise} is reset to be the standard deviation of the response, before adding noise.   (class "\code{numeric}")
#' @param ncovariates number of covariates to add to the model, numeric vector of length 1.
#' @param gamma True coefficient for the main effect of the modifier (class "\code{numeric}")
#' @return This returns a list of 8 items:
#' \item{x}{a lagged exposure matrix. If \code{x} was a matrix, it is unchanged. (class "\code{matrix}")}
#' \item{L}{a numeric vector of length 1 containing the number of lag terms (class "\code{numeric}")}
#' \item{modifiers}{the \code{modifiers} argument (class "\code{numeric}")}
#' \item{y}{a numeric vector of length \code{nrow(x)} containing the perturbed simulated response values. (class "\code{numeric}")}
#' \item{betas}{a matrix containing true coefficients for each lag/modifier combination, with each row representing a lag and each column a modifier (class "\code{matrix}")}
#' \item{betas_cumul}{a numeric vector of length \code{L+1} containing cumulative true coefficients for the lag terms, summed over modifiers (class "\code{numeric}")}
#' \item{Z}{covariates (class "\code{matrix}")}
#' \item{gammas}{true coefficients for the covariates (class "\code{numeric}")}



sim_data <- function(x, L=NULL, modifiers, noise=1, type=2, SNR, ncovariates=0, gamma=1){

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
