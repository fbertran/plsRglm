\name{kfolds2CVinfos_lm}
\alias{kfolds2CVinfos_lm}
\title{Extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares models}
\description{
This function extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares models for both formula or classic specifications of the model.
}
\usage{
kfolds2CVinfos_lm(pls_kfolds, MClassed = FALSE,verbose=TRUE)
}
\arguments{
  \item{pls_kfolds}{an object computed using \code{\link{PLS_lm_kfoldcv}}}
  \item{MClassed}{should number of miss classed be computed}
  \item{verbose}{should infos be displayed ?}
}
\details{
The Mclassed option should only set to \code{TRUE} if the response is binary. 
}
\value{
  \item{list}{table of fit statistics for first group partition}
  \item{\dots}{\dots}
  \item{list}{table of fit statistics for last group partition}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{summary}} and \code{\link{cv.plsR}} instead.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsR(dataY=yCornell,dataX=XCornell,nt=10,K=6)
summary(bbb)
rm(list=c("XCornell","yCornell","bbb"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb <- cv.plsR(dataY=ypine,dataX=Xpine,nt=10,NK=3)
summary(bbb)
bbb2 <- cv.plsR(dataY=ypine,dataX=XpineNAX21,nt=10,NK=3)
summary(bbb2)
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- cv.plsR(yaze_compl,Xaze_compl,nt=10,K=8,NK=3)
summary(bbb)
summary(bbb,MClassed=TRUE)
rm(list=c("Xaze_compl","yaze_compl","bbb"))
}
}\keyword{internal}
\keyword{models}
\keyword{regression}
