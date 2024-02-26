\name{dlim-package}
\alias{dlim-package}
\docType{package}

\title{
Distributed Lag Interaction Model (DLIM)
}

\description{
The package \pkg{dlim} contains functions to fit, perform inference and estimation on, and visualize a distributed lag interaction model (DLIM).
}

\section{Modelling framework}{
A distributed lag interaction model (DLIM) is an extension of a distributed lag model. A DLIM regresses an outcome onto repeated measures of an exposure and allows for associations to vary by a single continuous modifier. More details on methodology are provided in the reference listed below.
}

\section{Functions and data included in the package}{
To fit a DLIM using this package, use the \code{\link{dlim}} function, which calls the \code{\link{cross_basis}} function to create the cross-basis and estimates regression coefficients using \code{gam} from \strong{mgcv} package.

The \code{\link{predict.dlim}} S3 function provides point-wise or cumulative effect estimates and uncertainty measures. 

The \code{\link{plot_DLF}} and \code{\link{plot_cumulative}} functions provide plots of the modified distributed lag functions and the cumulative effect estimate curve. 
}

\section{Additional information}{
Additonal details on the package \pkg{dlim} are available in the vignette, available by typing:

\code{vignette("dlimOverview")}

The \pkg{dlim} package is available on the Comprehensive R Archive Network (CRAN). A development website is available on GitHub (\href{https://github.com/ddemateis/dlim}{github.com/ddemateis/dlim}).

Please use \code{citation("dlim")} to cite this package.
}

\author{
Danielle Demateis, Kayleigh Keller, and Ander Wilson

Maintainer: Danielle Demateis <\email{Danielle.Demateis@colostate.edu}>
}

\references{
Demateis et al. (2024) <doi:10.1002/env.2843>, avaibable at (\href{https://arxiv.org/abs/2401.02939}{arxiv.org/abs/2401.02939}).
}

\keyword{package}

\seealso{

Type \code{'vignette(dlimOverview)'} for a detailed description.
}