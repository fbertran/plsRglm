\name{bordeaux}
\alias{bordeaux}
\docType{data}
\title{Quality of wine dataset}
\description{
Quality of Bordeaux wines (\code{Quality}) and four potentially predictive variables (\code{Temperature}, \code{Sunshine}, \code{Heat} and \code{Rain}).
}
\usage{data(bordeaux)}
\format{
  A data frame with 34 observations on the following 5 variables.
  \describe{
    \item{\code{Temperature}}{a numeric vector}
    \item{\code{Sunshine}}{a numeric vector}
    \item{\code{Heat}}{a numeric vector}
    \item{\code{Rain}}{a numeric vector}
    \item{\code{Quality}}{an ordered factor with levels \code{1} < \code{2} < \code{3}}
  }
}
\source{
P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS generalised linear regression. \emph{Computational Statistics & Data Analysis}, 48(1):17-46.
}
\references{
M. Tenenhaus. (2005). La regression logistique PLS. In J.-J. Droesbeke, M. Lejeune, and G. Saporta, editors, Modeles statistiques pour donnees qualitatives. Editions Technip, Paris.}
\examples{
data(bordeaux)
str(bordeaux)
}
\keyword{datasets}
