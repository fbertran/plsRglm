\name{CorMat}
\alias{CorMat}
\docType{data}
\title{Correlation matrix for simulating plsR datasets}
\description{
A correlation matrix to simulate datasets
}
\usage{data(CorMat)}
\format{
  A data frame with 17 observations on the following 17 variables.
  \describe{
    \item{\code{y}}{a numeric vector}
    \item{\code{x11}}{a numeric vector}
    \item{\code{x12}}{a numeric vector}
    \item{\code{x13}}{a numeric vector}
    \item{\code{x21}}{a numeric vector}
    \item{\code{x22}}{a numeric vector}
    \item{\code{x31}}{a numeric vector}
    \item{\code{x32}}{a numeric vector}
    \item{\code{x33}}{a numeric vector}
    \item{\code{x34}}{a numeric vector}
    \item{\code{x41}}{a numeric vector}
    \item{\code{x42}}{a numeric vector}
    \item{\code{x51}}{a numeric vector}
    \item{\code{x61}}{a numeric vector}
    \item{\code{x62}}{a numeric vector}
    \item{\code{x63}}{a numeric vector}
    \item{\code{x64}}{a numeric vector}
  }
}
\source{Handmade.}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\examples{
data(CorMat)
str(CorMat)
}
\keyword{datasets}
