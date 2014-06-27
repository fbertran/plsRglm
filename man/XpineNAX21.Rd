\name{XpineNAX21}
\alias{XpineNAX21}
\docType{data}
\title{Incomplete dataset from the pine caterpillars example}
\description{
The caterpillar dataset was extracted from a 1973 study on pine processionary caterpillars. It assesses the influence of some forest settlement characteristics on the development of caterpillar colonies. There are k=10 potentially explanatory variables defined on n=33 areas.\cr
The value of x2 for the first observation was remove from the matrix of predictors on purpose.}
\usage{data(XpineNAX21)}
\format{
  A data frame with 33 observations on the following 10 variables.
  \describe{
    \item{\code{x1}}{altitude (in meters)}
    \item{\code{x2}}{slope (en degrees)}
    \item{\code{x3}}{number of pines in the area}
    \item{\code{x4}}{height (in meters) of the tree sampled at the center of the area}
    \item{\code{x5}}{diameter (in meters) of the tree sampled at the center of the area}
    \item{\code{x6}}{index of the settlement density}
    \item{\code{x7}}{orientation of the area (from 1 if southbound to 2 otherwise)}
    \item{\code{x8}}{height (in meters) of the dominant tree}
    \item{\code{x9}}{number of vegetation strata}
    \item{\code{x10}}{mix settlement index (from 1 if not mixed to 2 if mixed)}
  }
}
\details{
These caterpillars got their names from their habit of moving over the ground in incredibly long head-to-tail processions when leaving their nest to create a new colony.\cr
The \code{XpineNAX21} is a dataset with a missing value for testing purpose.
}
\source{
Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992). \dQuote{La \enc{regression}{r\'egression}, nouveaux regards sur une ancienne \enc{methode}{m\'ethode} statistique}, INRA, \emph{\enc{Actualit?s}{Actualit\'es} Scientifiques et Agronomiques}, Masson, Paris.}
%\references{}
\examples{
data(XpineNAX21)
str(XpineNAX21)
}
\keyword{datasets}
