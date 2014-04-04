\name{kfolds2CVinfos_glm}
\alias{kfolds2CVinfos_glm}
\title{Extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares glm models}
\description{
This function extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares glm models for both formula or classic specifications of the model.
}
\usage{
kfolds2CVinfos_glm(pls_kfolds, MClassed = FALSE)
}
\arguments{
  \item{pls_kfolds}{an object computed using \code{\link{PLS_glm_kfoldcv}}}
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
%\note{Use \code{\link{kfolds2CVinfos_glm}} and \code{\link{PLS_glm_kfoldcv}} instead.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),
nt=6,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
kfolds2CVinfos_glm(bbb)
rm(list=c("XCornell","yCornell","bbb"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,
keepdataY=TRUE,modele="pls")
kfolds2CVinfos_glm(bbb,MClassed=TRUE)
bbba <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,
keepdataY=TRUE,modele="pls-glm-family",family=gaussian())
kfolds2CVinfos_glm(bbba,MClassed=TRUE)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,
keepdataY=TRUE,modele="pls-glm-logistic")
kfolds2CVinfos_glm(bbb2,MClassed=TRUE)
bbb2a <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,
keepdataY=TRUE,modele="pls-glm-family",family=binomial())
kfolds2CVinfos_glm(bbb2a,MClassed=TRUE)
rm(list=c("Xaze_compl","yaze_compl","bbb","bbb2","bbba","bbb2a"))


#install.packages("chemometrics")
library(chemometrics)
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
bbb <- PLS_glm_kfoldcv(yhyptis,Xhyptis,nt=4,K=10,random=TRUE,modele="pls-glm-polr",
keepcoeffs=T)
kfolds2CVinfos_v2(bbb,MClassed=TRUE)
}
}\keyword{models}
\keyword{regression}
