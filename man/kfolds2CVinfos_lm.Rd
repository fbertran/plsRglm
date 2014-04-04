\name{kfolds2CVinfos_lm}
\alias{kfolds2CVinfos_lm}
\title{Extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares models}
\description{
This function extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares models for both formula or classic specifications of the model.
}
\usage{
kfolds2CVinfos_lm(pls_kfolds, MClassed = FALSE)
}
\arguments{
  \item{pls_kfolds}{an object computed using \code{\link{PLS_lm_kfoldcv}}}
  \item{MClassed}{should number of miss classed be computed}
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
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{r?gression}{r\'egression} PLS et de la \enc{r?gression}{r\'egression} logistique PLS : application aux \enc{donn?es}{donn\'ees} \enc{d'all?lotypage}{d'all\'elotypage}. \emph{Journal de la Soci?t? Fran?aise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
%\note{Use \code{\link{kfolds2CVinfos_lm}} and \code{\link{PLS_lm_kfoldcv}} instead.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb,MClassed=TRUE)
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,keepcoeffs=TRUE)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb,MClassed=TRUE)
rm(list=c("XCornell","yCornell","bbb"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=10,K=12,NK=3,
keepfolds=FALSE,keepdataY=TRUE)
kfolds2CVinfos_lm(bbb)
bbb2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=XpineNAX21,nt=10,K=12,NK=3,
random=TRUE,keepfolds=FALSE,keepdataY=TRUE)
kfolds2CVinfos_lm(bbb2)
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- PLS_lm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=3,keepdataY=TRUE)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb,MClassed=TRUE)
rm(list=c("Xaze_compl","yaze_compl","bbb"))
}
}\keyword{internal}
\keyword{models}
\keyword{regression}
