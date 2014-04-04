\name{plsRglm}
\alias{plsRglm}
\alias{plsRglmmodel.default}
\alias{plsRglmmodel.formula}
\title{Partial least squares Regression generalized linear models}
\description{
This function implements Partial least squares Regression generalized linear models complete or incomplete datasets.
}
\usage{
plsRglm(x, \dots)
\method{plsRglmmodel}{default}(dataY,dataX,nt=2,limQ2set=.0975,
dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",
EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,
alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,
sparse=FALSE,sparseStop=TRUE,naive=FALSE)
\method{plsRglmmodel}{formula}(formula,data=NULL,nt=2,limQ2set=.0975,
dataPredictY,modele="pls",family=NULL,typeVC="none",
EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,
alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,subset,
start=NULL,etastart,mustart,offset,method="glm.fit",control= list(),
contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE)
}
\arguments{
  \item{x}{a formula or a response (training) dataset}
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{formula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{plsRglm} is called.}  
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{dataPredictY}{predictor(s) (testing) dataset}
  \item{modele}{name of the PLS glm model to be fitted (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}). Use \code{"modele=pls-glm-family"} to enable the \code{family} option.}
  \item{family}{a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link[stats]{family}} for details of family functions.) To use the family option, please set \code{modele="pls-glm-family"}. User defined families can also be defined. See details.}
  \item{typeVC}{type of leave one out cross validation. For back compatibility purpose.
    \describe{
      \item{\code{none}}{no cross validation}
      }
  }
  \item{EstimXNA}{only for \code{modele="pls"}. Set whether the missing X values have to be estimated.}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{pvals.expli}{should individual p-values be reported to tune model selection ?}
  \item{alpha.pvals.expli}{level of significance for predictors when pvals.expli=TRUE}
  \item{MClassed}{number of missclassified cases, should only be used for binary responses}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{start}{starting values for the parameters in the linear predictor.}
  \item{etastart}{starting values for the linear predictor.}
  \item{mustart}{starting values for the vector of means.}
  \item{offset}{this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more \code{\link{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.}
  \item{method}{the method to be used in fitting the model. The default method \code{"glm.fit"} uses iteratively reweighted least squares (IWLS). User-supplied fitting functions can be supplied either as a function or a character string naming a function, with a function which takes the same arguments as \code{glm.fit}.}
  \item{control}{a list of parameters for controlling the fitting process. For \code{glm.fit} this is passed to \code{\link{glm.control}}.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  \item{sparse}{should the coefficients of non-significant predictors (<\code{alpha.pvals.expli}) be set to 0}
  \item{sparseStop}{should component extraction stop when no significant predictors (<\code{alpha.pvals.expli}) are found}
  \item{naive}{Use the naive estimates for the Degrees of Freedom in plsR? Default is \code{FALSE}.}
  \item{\dots}{arguments to pass to \code{plsRmodel.default} or to \code{plsRmodel.formula}}
}
\details{
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

A typical predictor has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with any duplicates removed. 

A specification of the form first:second indicates the the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second. 

The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a terms object as the formula. 

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 

The default estimator for Degrees of Freedom is the Kramer and Sugiyama's one which only works for classical plsR models. For these models, Information criteria are computed accordingly to these estimations. Naive Degrees of Freedom and Information Criteria are also provided for comparison purposes. For more details, see Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, http://arxiv.org/abs/1002.4112. 
}
\value{
Depends on the model that was used to fit the model.
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{
Use \code{plsRglm} instead.
}
\seealso{\code{\link{plsR}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
plsRglm(yCornell,XCornell,3)$uscores
plsRglm(yCornell,XCornell,3)$pp
plsRglm(yCornell,XCornell,3)$Coeffs

\donttest{
plsRglm(yCornell,XCornell,10)$InfCrit
plsRglm(yCornell,XCornell,10,modele="pls-glm-gaussian")$InfCrit
rm(list=c("XCornell","yCornell"))

data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
plsRglm(log(ypine),Xpine,1)$Std.Coeffs
plsRglm(log(ypine),Xpine,1)$Coeffs
plsRglm(log(ypine),Xpine,4)$Std.Coeffs
plsRglm(log(ypine),Xpine,4)$Coeffs
plsRglm(log(ypine),Xpine,4)$PredictY[1,]
plsRglm(log(ypine),Xpine,4,dataPredictY=Xpine[1,])$PredictY[1,]

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
str(plsRglm(log(ypine),XpineNAX21,2))
plsRglm(log(ypine),XpineNAX21,4)$Std.Coeffs
plsRglm(log(ypine),XpineNAX21,4)$YChapeau[1,]
plsRglm(log(ypine),Xpine,4)$YChapeau[1,]
plsRglm(log(ypine),XpineNAX21,4)$CoeffC
plsRglm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeau
plsRglm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA

# compare pls-glm-gaussian with classic plsR
cbind(plsRglm(log(ypine),Xpine,4,modele="pls")$Std.Coeffs,plsRglm(log(ypine),
Xpine,4,modele="pls-glm-gaussian")$Std.Coeffs)

# without missing data
cbind(log(ypine),plsRglm(log(ypine),Xpine,4,modele="pls")$YChapeau,
plsRglm(log(ypine),Xpine,4,modele="pls-glm-gaussian")$YChapeau)
cbind(log(ypine),plsRglm(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,
plsRglm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau)

# with missing data
cbind((log(ypine)),plsRglm(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,
plsRglm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau)
cbind((log(ypine)),plsRglm(log(ypine),XpineNAX21,4,modele="pls")$ValsPredictY,
plsRglm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$ValsPredictY)
rm(list=c("Xpine","ypine"))

data(fowlkes)
Xfowlkes <- fowlkes[,2:13]
yfowlkes <- fowlkes[,1]
modpls <- plsRglm(yfowlkes,Xfowlkes,4,modele="pls-glm-logistic",pvals.expli=TRUE)
modpls$pvalstep
rm(list=c("Xfowlkes","yfowlkes","modpls"))
}

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",
MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
modpls$valpvalstep
modpls$Coeffsmodel_vals

plot(plsRglm(yaze_compl,Xaze_compl,4,modele="pls-glm-logistic")$FinalModel)
plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,
modele="pls-glm-logistic",pvals.expli=TRUE)$pvalstep
plot(plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,
modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)
rm(list=c("Xaze_compl","yaze_compl","modpls"))


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
modpls <- plsRglm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit

XbordeauxNA<-Xbordeaux
XbordeauxNA[1,1] <- NA
modplsNA <- plsRglm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit
rm(list=c("Xbordeaux","XbordeauxNA","ybordeaux","modplsNA"))


\donttest{
#install.packages(chemometrics)
library(chemometrics)
data(hyptis)
yhyptis <- factor(hyptis$Group,ordered=TRUE)
Xhyptis <- as.data.frame(hyptis[,c(1:6)])
options(contrasts = c("contr.treatment", "contr.poly"))
modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
modpls2$Coeffsmodel_vals
modpls2$InfCrit
modpls2$Coeffs
modpls2$std.coeffs

table(yhyptis,predict(modpls2$FinalModel,type="class"))
rm(list=c("yhyptis","Xhyptis","modpls2"))


dimX <- 6
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 2
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar2)[,1]
Xsimbin1 <- dicho(dataAstar2)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar2","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 3
dataAstar3 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar3)[,1]
Xsimbin1 <- dicho(dataAstar3)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar3","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 5
dataAstar5 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar5)[,1]
Xsimbin1 <- dicho(dataAstar5)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar5","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 6
dataAstar6 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar6)[,1]
Xsimbin1 <- dicho(dataAstar6)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar6","ysimbin1","Xsimbin1","modplsglm"))

}
}\keyword{models}
\keyword{regression}
