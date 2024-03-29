% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.summary.plsRglmmodel.R
\name{print.summary.plsRglmmodel}
\alias{print.summary.plsRglmmodel}
\title{Print method for summaries of plsRglm models}
\usage{
\method{print}{summary.plsRglmmodel}(x, ...)
}
\arguments{
\item{x}{an object of the class \code{"summary.plsRglmmodel"}}

\item{\dots}{not used}
}
\value{
\item{language}{call of the model}
}
\description{
This function provides a print method for the class
\code{"summary.plsRglmmodel"}
}
\examples{

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
class(modplsglm)
print(summary(modplsglm))
rm(list=c("XCornell","yCornell","modplsglm"))

}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et
Frédéric Bertrand (2010). Comparaison de la
régression PLS et de la régression
logistique PLS : application aux données
d'allélotypage. \emph{Journal de la Société Française
de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\seealso{
\code{\link{print}} and \code{\link{summary}}
}
\author{
Frédéric Bertrand\cr
\email{frederic.bertrand@utt.fr}\cr
\url{https://fbertran.github.io/homepage/}
}
\keyword{methods}
\keyword{print}
