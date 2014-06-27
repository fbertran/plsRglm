\name{print.cv.plsRmodel}
\alias{print.cv.plsRmodel}
\title{Print method for plsR models}
\description{
This function provides a print method for the class \code{"cv.plsRmodel"}
}
\usage{
\method{print}{cv.plsRmodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"cv.plsRmodel"}}
  \item{\dots}{not used}
}
%\details{}
\value{
\code{NULL}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsR(dataY=yCornell,dataX=XCornell,nt=10,K=6)
print(bbb)
rm(list=c("XCornell","yCornell","bbb"))
}
\keyword{methods}
\keyword{print}
