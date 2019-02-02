\name{plot.table.summary.cv.plsRmodel}
\alias{plot.table.summary.cv.plsRmodel}
\title{Plot method for table of summary of cross validated plsR models}
\description{
This function provides a table method for the class \code{"summary.cv.plsRmodel"}
}
\usage{
\method{plot}{table.summary.cv.plsRmodel}(x, type=c("CVMC","CVQ2","CVPress"), ...)
}
\arguments{
  \item{x}{an object of the class \code{"table.summary.cv.plsRmodel"}}
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
bbb <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=5, verbose=FALSE)
plot(cvtable(summary(bbb)),type="CVQ2")
rm(list=c("bbb"))

\donttest{
data(Cornell)
plot(cvtable(summary(cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=100, verbose=FALSE))),type="CVQ2")
}
}
\keyword{methods}
\keyword{print}
