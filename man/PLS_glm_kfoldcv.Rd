\name{PLS_glm_kfoldcv}
\alias{PLS_glm_kfoldcv}
\title{Partial least squares regression glm models with k-fold cross validation}
\description{
This function implements k-fold cross-validation on complete or incomplete datasets for partial least squares regression generalized linear models
}
\usage{
PLS_glm_kfoldcv(dataY, dataX, nt = 2, limQ2set = 0.0975, modele = "pls", 
family = NULL, K = nrow(dataX), NK = 1, grouplist = NULL, random = FALSE, 
scaleX = TRUE, scaleY = NULL, keepcoeffs = FALSE, keepfolds = FALSE, 
keepdataY = TRUE, keepMclassed=FALSE, tol_Xi = 10^(-12), weights, method)
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{modele}{name of the PLS glm model to be fitted (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}). Use \code{"modele=pls-glm-family"} to enable the \code{family} option.}
  \item{family}{a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link[stats]{family}} for details of family functions.) To use the family option, please set \code{modele="pls-glm-family"}. User defined families can also be defined. See details.}
  \item{K}{number of groups}
  \item{NK}{number of times the group division is made}
  \item{grouplist}{to specify the members of the \code{K} groups}
  \item{random}{should the \code{K} groups be made randomly}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{keepcoeffs}{shall the coefficients for each model be returned}
  \item{keepfolds}{shall the groups' composition be returned}
  \item{keepdataY}{shall the observed value of the response for each one of the predicted value be returned}
  \item{keepMclassed}{shall the number of miss classed be returned (unavailable)}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{method}{for \code{pls-glm-polr} models, logistic, probit, complementary log-log (cloglog) or cauchit (corresponding to a Cauchy latent variable).}
}
\details{
Predicts 1 group with the \code{K-1} other groups. Leave one out cross validation is thus obtained for \code{K==nrow(dataX)}.

There are seven different predefined models with predefined link functions available :
\describe{
\item{\code{"pls"}}{ordinary pls models}
\item{\code{"pls-glm-Gamma"}}{glm gaussian with inverse link pls models}
\item{\code{"pls-glm-gaussian"}}{glm gaussian with identity link pls models}
\item{\code{"pls-glm-inverse-gamma"}}{glm binomial with square inverse link pls models}
\item{\code{"pls-glm-logistic"}}{glm binomial with logit link pls models}
\item{\code{"pls-glm-poisson"}}{glm poisson with log link pls models}
\item{\code{"pls-glm-polr"}}{glm polr with logit link pls models}
}
Using the \code{"family="} option and setting \code{"modele=pls-glm-family"} allows changing the family and link function the same way as for the \code{\link[stats]{glm}} function. As a consequence user-specified families can also be used. 
\describe{
\item{The \code{gaussian} family}{accepts the links (as names) \code{identity}, \code{log} and \code{inverse}.}
\item{The \code{binomial} family}{accepts the links \code{logit}, \code{probit}, \code{cauchit}, (corresponding to logistic, normal and Cauchy CDFs respectively) \code{log} and \code{cloglog} (complementary log-log).} 
\item{The \code{Gamma} family}{accepts the links \code{inverse}, \code{identity} and \code{log}.}
\item{The \code{poisson} family}{accepts the links \code{log}, \code{identity}, and \code{sqrt}.}
\item{The \code{inverse.gaussian} family}{accepts the links \code{1/mu^2}, \code{inverse}, \code{identity} and \code{log}.}
\item{The \code{quasi} family}{accepts the links \code{logit}, \code{probit}, \code{cloglog}, \code{identity}, \code{inverse}, \code{log}, \code{1/mu^2} and \code{sqrt}.}
\item{The function \code{power}}{can be used to create a power link function.} 
}

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 
}
\value{
  \item{results_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    \item{\dots}{\dots}
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    }}
  \item{folds}{list of \code{NK}. Each element of the list sums up the informations for a group division:
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
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}}, \code{\link{kfolds2Mclassed}} and \code{\link{kfolds2CVinfos_glm}} to extract and transform results from k-fold cross validation.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)

PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=6,
NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds

\donttest{
#Different ways of model specifications
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=6,
NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
family=gaussian,K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
family=gaussian(),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
family=gaussian(link=log),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds

bbb2 <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,
modele="pls-glm-gaussian",keepcoeffs=TRUE)
bbb2 <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-family",family=gaussian(link=log),K=6,keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(yCornell),XCornell,10,typeVC="standard")$CVinfos
rm(list=c("XCornell","yCornell","bbb","bbb2"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-family",
family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
 bb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,
 modele="pls-glm-gaussian",K=10,keepcoeffs=TRUE,keepfolds=FALSE)
b
#For Jackknife computations
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])

kfolds2Chisqind(bbb)
kfolds2Chisq(bbb)
kfolds2CVinfos_glm(bbb)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$CVinfos

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-family"
,family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,
modele="pls-glm-gaussian",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$CVinfos
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=10,modele="pls",keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=10,modele="pls-glm-family",
family=binomial(probit),keepcoeffs=TRUE)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=10,
modele="pls-glm-logistic",keepcoeffs=TRUE)
kfolds2CVinfos_v2(bbb,MClassed=TRUE)
kfolds2CVinfos_v2(bbb2,MClassed=TRUE)
kfolds2coeff(bbb2)

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
rm(list=c("Xaze_compl","yaze_compl","bbb","bbb2"))



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=Xpine,nt=10,modele="pls-glm-family",
family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=Xpine,nt=10,
modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])

kfolds2Chisqind(bbb)
kfolds2Chisq(bbb)
kfolds2CVinfos_glm(bbb)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$CVinfos

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb2 <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=XpineNAX21,nt=10,
modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb2 <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=XpineNAX21,nt=10,
modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$CVinfos
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-family",
family=Gamma,K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-Gamma",
K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])

kfolds2Chisqind(bbb)
kfolds2Chisq(bbb)
kfolds2CVinfos_glm(bbb)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$CVinfos

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,
modele="pls-glm-family",family=Gamma(),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,
modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$CVinfos
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))



data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)

PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-inverse.gaussian",K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-family",family=inverse.gaussian,K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-inverse.gaussian",K=6,NK=2,random=TRUE,
keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-family",family=inverse.gaussian(),K=6,NK=2,random=TRUE,
keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-inverse.gaussian",K=6,NK=2,random=FALSE,
keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,
modele="pls-glm-family",family=inverse.gaussian(link = "1/mu^2"),K=6,NK=2,
random=FALSE,keepfolds=TRUE)$results_kfolds

bbb2 <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,
modele="pls-glm-inverse.gaussian",keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(yCornell),XCornell,10,typeVC="standard")$CVinfos
rm(list=c("XCornell","yCornell","bbb","bbb2"))
}
}
\keyword{models}
\keyword{regression}
