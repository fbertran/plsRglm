\name{kfolds2CVinfos_glm}
\alias{kfolds2CVinfos_glm}
\title{Extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares glm models}
\description{
This function extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares glm models for both formula or classic specifications of the model.
}
\usage{
kfolds2CVinfos_glm(pls_kfolds, MClassed = FALSE,verbose=TRUE)
}
\arguments{
  \item{pls_kfolds}{an object computed using \code{\link{cv.plsRglm}}}
  \item{MClassed}{should number of miss classed be computed ?}
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
\note{Use \code{\link{summary}} and \code{\link{cv.plsRglm}} instead.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{
\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsRglm(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),
nt=6,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
summary(bbb,MClassed=TRUE)
rm(list=c("XCornell","yCornell","bbb"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- cv.plsR(yaze_compl,Xaze_compl,nt=10,K=8,modele="pls")
summary(bbb,MClassed=TRUE)
bbbbis <- cv.plsRglm(yaze_compl,Xaze_compl,nt=10,K=8,modele="pls")
summary(bbbbis,MClassed=TRUE)
bbba <- cv.plsRglm(yaze_compl,Xaze_compl,nt=10,K=8,
modele="pls-glm-family",family=gaussian())
summary(bbba,MClassed=TRUE)
bbb2 <- cv.plsRglm(yaze_compl,Xaze_compl,nt=10,K=8,
modele="pls-glm-logistic")
summary(bbb2,MClassed=TRUE)
bbb2a <- cv.plsRglm(yaze_compl,Xaze_compl,nt=10,K=8,
modele="pls-glm-family",family=binomial())
summary(bbb2a,MClassed=TRUE)
rm(list=c("Xaze_compl","yaze_compl","bbb","bbbbis","bbb2","bbba","bbb2a"))


if(require(chemometrics)){
data(hyptis)
hyptis
yhyptis <- factor(hyptis$Group,ordered=TRUE)
Xhyptis <- as.data.frame(hyptis[,c(1:6)])
options(contrasts = c("contr.treatment", "contr.poly"))
modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
modpls2$Coeffsmodel_vals
modpls2$InfCrit
modpls2$Coeffs
modpls2$std.coeffs

table(yhyptis,predict(modpls2$FinalModel,type="class"))

modpls3 <- PLS_glm(yhyptis[-c(1,2,3)],Xhyptis[-c(1,2,3),],3,modele="pls-glm-polr",
dataPredictY=Xhyptis[c(1,2,3),])
bbb <- cv.plsRglm(yhyptis,Xhyptis,nt=4,K=10,random=TRUE,modele="pls-glm-polr",
keepcoeffs=TRUE)
summary(bbb,MClassed=TRUE)
}
}
}
\keyword{models}
\keyword{regression}
