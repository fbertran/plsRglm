\name{Cornell}
\alias{Cornell}
\docType{data}
\title{Cornell dataset}
\description{
The famous Cornell dataset. A mixture experiment on \code{X1}, \code{X2}, \code{X3}, \code{X4}, \code{X5}, \code{X6} and \code{X7} to analyse octane degree (\code{Y}) in gazoline.
}
\usage{data(Cornell)}
\format{
  A data frame with 12 observations on the following 8 variables.
  \describe{
    \item{\code{X1}}{a numeric vector}
    \item{\code{X2}}{a numeric vector}
    \item{\code{X3}}{a numeric vector}
    \item{\code{X4}}{a numeric vector}
    \item{\code{X5}}{a numeric vector}
    \item{\code{X6}}{a numeric vector}
    \item{\code{X7}}{a numeric vector}
    \item{\code{Y}}{response value: a numeric vector}
  }
}
\source{
M. Tenenhaus. (1998). \emph{La regression PLS, Theorie et pratique}. Editions Technip, Paris.
}
\references{
N. Kettaneh-Wold. Analysis of mixture data with partial least squares. (1992). \emph{Chemometrics and Intelligent Laboratory Systems}, 14(1):57-69.
}
\examples{
data(Cornell)
str(Cornell)
}
\keyword{datasets}
