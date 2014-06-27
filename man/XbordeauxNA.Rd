\name{XbordeauxNA}
\alias{XbordeauxNA}
\docType{data}
\title{Incomplete dataset for the quality of wine dataset}
\description{
Quality of Bordeaux wines (\code{Quality}) and four potentially predictive variables (\code{Temperature}, \code{Sunshine}, \code{Heat} and \code{Rain}).\cr
The value of Temperature for the first observation was remove from the matrix of predictors on purpose.
}
\usage{data(XbordeauxNA)}
\format{
  A data frame with 34 observations on the following 4 variables.
  \describe{
    \item{\code{Temperature}}{a numeric vector}
    \item{\code{Sunshine}}{a numeric vector}
    \item{\code{Heat}}{a numeric vector}
    \item{\code{Rain}}{a numeric vector}
  }
}
\source{
P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS generalised linear regression. \emph{Computational Statistics & Data Analysis}, 48(1):17-46.
}
\references{
M. Tenenhaus. (2005). La regression logistique PLS. In J.-J. Droesbeke, M. Lejeune, and G. Saporta, editors, Modeles statistiques pour donnees qualitatives. Editions Technip, Paris.}
\examples{
data(XbordeauxNA)
str(XbordeauxNA)
}
\keyword{datasets}
