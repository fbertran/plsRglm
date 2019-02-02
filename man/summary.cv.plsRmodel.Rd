\name{summary.cv.plsRmodel}
\alias{summary.cv.plsRmodel}
\title{Summary method for plsR models}
\description{
This function provides a summary method for the class \code{"cv.plsRmodel"}
}
\usage{
\method{summary}{cv.plsRmodel}(object, \dots)
}
\arguments{
  \item{object}{an object of the class \code{"cv.plsRmodel"}}
  \item{\dots}{further arguments to be passed to or from methods.}  
}
%\details{}
\value{An object of class \code{"summary.cv.plsRglmmodel"}.}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{summary}}}
\examples{
data(Cornell)
summary(cv.plsR(Y~.,data=Cornell,nt=10,K=6, verbose=FALSE), verbose=FALSE)
}
\keyword{methods}
\keyword{print}
