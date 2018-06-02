\name{plot.table.summary.cv.plsRglmmodel}
\alias{plot.table.summary.cv.plsRglmmodel}
\title{Plot method for table of summary of cross validated plsRglm models}
\description{
This function provides a table method for the class \code{"summary.cv.plsRglmmodel"}
}
\usage{
\method{plot}{table.summary.cv.plsRglmmodel}(x, type=c("CVMC","CVQ2Chi2","CVPreChi2"), ...)
}
\arguments{
  \item{x}{an object of the class \code{"table.summary.cv.plsRglmmodel"}}
  \item{type}{the type of cross validation criterion to plot.}  
  \item{\dots}{further arguments to be passed to or from methods.}  
}
%\details{}
\value{
  \code{NULL}
}
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
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=10,NK=1,
modele="pls-glm-family",family=gaussian())
plot(cvtable(summary(bbb)),type="CVQ2Chi2")
rm(list=c("XCornell","yCornell","bbb"))

\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=10,NK=1,
modele="pls-glm-family",family=gaussian())
plot(cvtable(summary(bbb)),type="CVQ2Chi2")
rm(list=c("XCornell","yCornell","bbb"))
}
}
\keyword{methods}
\keyword{print}
