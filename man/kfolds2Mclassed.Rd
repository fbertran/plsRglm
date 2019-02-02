\name{kfolds2Mclassed}
\alias{kfolds2Mclassed}
\title{Number of missclassified individuals for k-fold cross validated partial least squares regression models.}
\description{This function indicates the total number of missclassified individuals for k-fold cross validated partial least squares regression models.
}
\usage{
kfolds2Mclassed(pls_kfolds)
}
\arguments{
  \item{pls_kfolds}{a k-fold cross validated partial least squares regression model used on binary data}
}
%\details{}
\value{
  \item{list}{Total number of missclassified individuals vs number of components for the first group partition}
  \item{\dots}{\dots}
  \item{list}{Total number of missclassified individuals vs number of components for the last group partition}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{cv.plsR}} to create k-fold cross validated partial least squares regression models.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Pressind}} and \code{\link{kfolds2Mclassedind}} to extract and transforms results from k-fold cross validation.}
\examples{
\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
kfolds2Mclassed(cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=1,verbose=FALSE))
kfolds2Mclassed(cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=2,verbose=FALSE))
rm(list=c("Xaze_compl","yaze_compl"))
}
}\keyword{models}
\keyword{regression}
