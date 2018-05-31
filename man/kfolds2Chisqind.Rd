\name{kfolds2Chisqind}
\alias{kfolds2Chisqind}
\title{Computes individual Predicted Chisquare for k-fold cross validated partial least squares regression models.}
\description{This function computes individual Predicted Chisquare for k-fold cross validated partial least squares regression models.}
\usage{
kfolds2Chisqind(pls_kfolds)
}
\arguments{
  \item{pls_kfolds}{a k-fold cross validated partial least squares regression glm model}
}
%\details{}
\value{
  \item{list}{Individual PChisq vs number of components for the first group partition}
  \item{\dots}{\dots}
  \item{list}{Individual PChisq vs number of components for the last group partition}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{cv.plsRglm}} to create k-fold cross validated partial least squares regression glm models.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Chisq}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=16)
bbb2 <- cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=5)
kfolds2Chisqind(bbb)
kfolds2Chisqind(bbb2)
rm(list=c("XCornell","yCornell","bbb","bbb2"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsRglm(dataY=ypine,dataX=Xpine,nt=4,modele="pls-glm-gaussian")
bbb2 <- cv.plsRglm(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-gaussian",K=10)
kfolds2Chisqind(bbb)
kfolds2Chisqind(bbb2)
                  
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- cv.plsRglm(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls",K=10)
kfolds2Pressind(bbbNA)
kfolds2Chisqind(bbbNA)
bbbNA2 <- cv.plsRglm(dataY=ypine,dataX=XpineNAX21,nt=4,modele="pls-glm-gaussian")
bbbNA3 <- cv.plsRglm(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian",K=10)
kfolds2Chisqind(bbbNA2)
kfolds2Chisqind(bbbNA3)
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2","bbbNA3"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",
family=binomial())
bbb <- cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic")
bbb2 <- cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-family",
family=binomial(),K=10)
bbb2 <- cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-logistic",K=10)
kfolds2Chisqind(bbb)
kfolds2Chisqind(bbb2)
rm(list=c("Xaze_compl","yaze_compl","bbb","bbb2"))
}
}\keyword{models}
\keyword{regression}
