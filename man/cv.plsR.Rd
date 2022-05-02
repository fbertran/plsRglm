\name{cv.plsR}
\alias{cv.plsR}
\alias{cv.plsRmodel.default}
\alias{cv.plsRmodel.formula}
\alias{PLS_lm_kfoldcv}
\alias{PLS_lm_kfoldcv_formula}
\title{Partial least squares regression models with k-fold cross-validation}
\description{
This function implements k-fold cross-validation on complete or incomplete datasets for partial least squares regression models
}
\usage{
cv.plsR(x, \dots)
\method{cv.plsRmodel}{default}(dataY,dataX,nt=2,limQ2set=.0975,modele="pls", 
K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, 
scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12), weights, verbose=TRUE)
\method{cv.plsRmodel}{formula}(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", 
K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, 
scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL, verbose=TRUE)
PLS_lm_kfoldcv(dataY, dataX, nt = 2, limQ2set = 0.0975, modele = "pls", 
K = 5, NK = 1, grouplist = NULL, random = TRUE, scaleX = TRUE, 
scaleY = NULL, keepcoeffs = FALSE, keepfolds = FALSE, keepdataY = TRUE, 
keepMclassed=FALSE, tol_Xi = 10^(-12), weights, verbose=TRUE)
PLS_lm_kfoldcv_formula(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", 
K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, 
scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL,verbose=TRUE)
}
\arguments{
  \item{x}{a formula or a response (training) dataset}
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{formula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{plsRglm} is called.}    
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{modele}{name of the PLS model to be fitted, only (\code{"pls"} available for this fonction.}
  \item{K}{number of groups. Defaults to 5.}
  \item{NK}{number of times the group division is made}
  \item{grouplist}{to specify the members of the \code{K} groups}
  \item{random}{should the \code{K} groups be made randomly. Defaults to \code{TRUE}}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{keepcoeffs}{shall the coefficients for each model be returned}
  \item{keepfolds}{shall the groups' composition be returned}
  \item{keepdataY}{shall the observed value of the response for each one of the predicted value be returned}
  \item{keepMclassed}{shall the number of miss classed be returned}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  \item{verbose}{should info messages be displayed ?}
\item{\dots}{arguments to pass to \code{cv.plsRmodel.default} or to \code{cv.plsRmodel.formula}}
}
\details{
Predicts 1 group with the \code{K-1} other groups. Leave one out cross validation is thus obtained for \code{K==nrow(dataX)}.

A typical predictor has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with any duplicates removed. 

A specification of the form first:second indicates the the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second. 

The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a terms object as the formula. 

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 
 }
\value{An object of class \code{"cv.plsRmodel"}.\cr 
  \item{results_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{list}{ of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    \item{\dots}{\dots}
    \item{list}{ of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    }
  }
  \item{folds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{list}{ of \code{K} vectors of length about \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were used as a training set}
    \item{\dots}{\dots}
    \item{list}{ of \code{K} vectors of length about \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were used as a training set}
    }
  }
  \item{dataY_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{list}{ of \code{K} matrices of size about \code{nrow(dataX)/K * 1} with the observed values of the response}
    \item{\dots}{\dots}
    \item{list}{ of \code{K} matrices of size about \code{nrow(dataX)/K * 1} with the observed values of the response}
    }
  }
  \item{call}{the call of the function}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@utt.fr}\cr
\url{https://fbertran.github.io/homepage/}
}
\note{Work for complete and incomplete datasets.}
\seealso{Summary method \code{summary.cv.plsRmodel}. \code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}}, \code{\link{kfolds2Mclassed}} and \code{\link{kfolds2CVinfos_lm}} to extract and transform results from k-fold cross-validation.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

#Leave one out CV (K=nrow(Cornell)) one time (NK=1)
bbb <- cv.plsR(dataY=yCornell,dataX=XCornell,nt=6,K=nrow(Cornell),NK=1)
bbb2 <- cv.plsR(Y~.,data=Cornell,nt=6,K=12,NK=1,verbose=FALSE)
(sum1<-summary(bbb2))

#6-fold CV (K=6) two times (NK=2)
#use random=TRUE to randomly create folds for repeated CV
bbb3 <- cv.plsR(dataY=yCornell,dataX=XCornell,nt=6,K=6,NK=2)
bbb4 <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=2,verbose=FALSE)
(sum3<-summary(bbb4))

cvtable(sum1)
cvtable(sum3)
rm(list=c("XCornell","yCornell","bbb","bbb2","bbb3","bbb4"))
}
\keyword{models}
\keyword{regression}
