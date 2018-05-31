\name{kfolds2Press}
\alias{kfolds2Press}
\title{Computes PRESS for k-fold cross validated partial least squares regression models.}
\description{This function computes PRESS for k-fold cross validated partial least squares regression models.}
\usage{
kfolds2Press(pls_kfolds)
}
\arguments{
  \item{pls_kfolds}{a k-fold cross validated partial least squares regression model}
}
%\details{}
\value{
  \item{list}{Press vs number of components for the first group partition}
  \item{\dots}{\dots}
  \item{list}{Press vs number of components for the last group partition}
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
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross validation.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsR(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1)
bbb2 <- cv.plsR(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=6,NK=1)
kfolds2Press(bbb)
kfolds2Press(bbb2)
rm(list=c("XCornell","yCornell","bbb","bbb2"))

\donttest{
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsR(dataY=ypine,dataX=Xpine,nt=10,NK=1)
bbb2 <- cv.plsR(dataY=ypine,dataX=Xpine,nt=10,NK=2)
kfolds2Press(bbb)
kfolds2Press(bbb2)

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- cv.plsR(dataY=ypine,dataX=XpineNAX21,nt=10,NK=1)
bbbNA2 <- cv.plsR(dataY=ypine,dataX=XpineNAX21,nt=10,NK=2)
kfolds2Press(bbbNA)
kfolds2Press(bbbNA2)
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2"))
}
}
\keyword{models}
\keyword{regression}
