\name{kfolds2Mclassedind}
\alias{kfolds2Mclassedind}
\title{Number of missclassified individuals per group for k-fold cross validated partial least squares regression models.}
\description{This function indicates the number of missclassified individuals per group for k-fold cross validated partial least squares regression models.
}
\usage{
kfolds2Mclassedind(pls_kfolds)
}
\arguments{
  \item{pls_kfolds}{a k-fold cross validated partial least squares regression model used on binary data}
}
%\details{}
\value{
  \item{list}{Number of missclassified individuals per group vs number of components for the first group partition}
  \item{\dots}{\dots}
  \item{list}{Number of missclassified individuals per group vs number of components for the last group partition}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{cv.plsR}} or \code{\link{cv.plsRglm}} to create k-fold cross validated partial least squares regression models or generalized linear ones.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Pressind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{
\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=1)
bbb2 <- cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=2)
kfolds2Mclassedind(bbb)
kfolds2Mclassedind(bbb2)
rm(list=c("Xaze_compl","yaze_compl","bbb","bbb2"))
}
}
\keyword{models}
\keyword{regression}
