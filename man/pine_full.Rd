\name{pine_full}
\alias{pine_full}
\docType{data}
\title{Complete Pine dataset}
\description{
This is the complete caterpillar dataset from a 1973 study on pine_full processionary caterpillars. It assesses the influence of some forest settlement characteristics on the development of caterpillar colonies. The response variable is the logarithmic transform of the average number of nests of caterpillars per tree in an area of 500 square meters (\code{x11}). There are k=10 potentially explanatory variables defined on n=55 areas.
}
\usage{data(pine_full)}
\format{
  A data frame with 55 observations on the following 11 variables.
  \describe{
    \item{\code{x1}}{altitude (in meters)}
    \item{\code{x2}}{slope (en degrees)}
    \item{\code{x3}}{number of pine_fulls in the area}
    \item{\code{x4}}{height (in meters) of the tree sampled at the center of the area}
    \item{\code{x5}}{diameter (in meters) of the tree sampled at the center of the area}
    \item{\code{x6}}{index of the settlement density}
    \item{\code{x7}}{orientation of the area (from 1 if southbound to 2 otherwise)}
    \item{\code{x8}}{height (in meters) of the dominant tree}
    \item{\code{x9}}{number of vegetation strata}
    \item{\code{x10}}{mix settlement index (from 1 if not mixed to 2 if mixed)}
    \item{\code{x11}}{logarithmic transform of the average number of nests of caterpillars per tree}
  }
}
\details{
These caterpillars got their names from their habit of moving over the ground in incredibly long head-to-tail processions when leaving their nest to create a new colony.
}
\source{
Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992), \dQuote{La \enc{regression}{r\'egression}, nouveaux regards sur une ancienne \enc{methode}{m\'ethode} statistique}, INRA, \emph{\enc{Actualites}{Actualit\'es} Scientifiques et Agronomiques}, Masson, Paris.
}
\references{
J.-M. Marin, C. Robert. (2007). \emph{Bayesian Core: A Practical Approach to Computational Bayesian Statistics}. Springer, New-York, pages 48-49.}
\examples{
data(pine_full)
str(pine_full)
}\keyword{datasets}
