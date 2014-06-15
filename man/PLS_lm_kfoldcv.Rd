\name{PLS_lm_kfoldcv}
\alias{PLS_lm_kfoldcv}
\title{Partial least squares regression models with k-fold cross-validation}
\description{
This function implements k-fold cross-validation on complete or incomplete datasets for partial least squares regression models
}
\usage{
PLS_lm_kfoldcv(dataY, dataX, nt = 2, limQ2set = 0.0975, modele = "pls", 
K = nrow(dataX), NK = 1, grouplist = NULL, random = FALSE, scaleX = TRUE, 
scaleY = NULL, keepcoeffs = FALSE, keepfolds = FALSE, keepdataY = TRUE, 
keepMclassed=FALSE, tol_Xi = 10^(-12), weights)
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{modele}{name of the PLS model to be fitted, only (\code{"pls"} available for this fonction.}
  \item{K}{number of groups}
  \item{NK}{number of times the group division is made}
  \item{grouplist}{to specify the members of the \code{K} groups}
  \item{random}{should the \code{K} groups be made randomly}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{keepcoeffs}{shall the coefficients for each model be returned}
  \item{keepfolds}{shall the groups' composition be returned}
  \item{keepdataY}{shall the observed value of the response for each one of the predicted value be returned}
  \item{keepMclassed}{shall the number of miss classed be returned}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
}
\details{
Predicts 1 group with the \code{K-1} other groups. Leave one out cross validation is thus obtained for \code{K==nrow(dataX)}.

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 
 }
\value{
  \item{results_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    \item{\dots}{\dots}
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    }
  }
  \item{folds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{}{list of \code{K} vectors of length about \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were used as a training set}
    \item{\dots}{\dots}
    \item{}{list of \code{K} vectors of length about \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were used as a training set}
    }
  }
  \item{dataY_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * 1} with the observed values of the response}
    \item{\dots}{\dots}
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * 1} with the observed values of the response}
    }
  }
  \item{call}{the call of the function}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{r?gression}{r\'egression} PLS et de la \enc{r?gression}{r\'egression} logistique PLS : application aux \enc{donn?es}{donn\'ees} \enc{d'all?lotypage}{d'all\'elotypage}. \emph{Journal de la Soci?t? Fran?aise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Work for complete and incomplete datasets.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}}, \code{\link{kfolds2Mclassed}} and \code{\link{kfolds2CVinfos_lm}} to extract and transform results from k-fold cross-validation.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=12,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=12,keepfolds=FALSE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=6,NK=2,random=FALSE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=6,NK=2,random=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE,keepfolds=FALSE)

bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm(yCornell,XCornell,6,typeVC="standard")$InfCrit
rm(list=c("XCornell","yCornell","bbb","bbb2"))


\donttest{
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=10,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=10,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=XpineNAX21,nt=10,K=12,NK=1)
bbbNA2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=XpineNAX21,nt=10,K=6,NK=1)
kfolds2CVinfos_lm(bbbNA)
kfolds2CVinfos_lm(bbbNA2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$InfCrit
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2"))
}
}
\keyword{models}
\keyword{regression}
