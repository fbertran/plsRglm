\name{PLS_glm_formula}
\alias{PLS_glm_formula}
\title{Partial least squares Regression generalized linear models}
\description{
This function implements Partial least squares Regression generalized linear models complete or incomplete datasets (formula specification of the model).
}
\usage{
PLS_glm_formula(formula,data=NULL,nt=2,limQ2set=.0975,dataPredictY=dataX,
modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,
scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,
tol_Xi=10^(-12),weights,subset,start=NULL,etastart,mustart,offset,method,
control= list(),contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE)
}
\arguments{
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
  \item{scaleY}{scale the response : Yes/No. Ignored since not always possible for glm responses.}
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
  \item{method}{\describe{
  \item{for fitting glms with glm (\code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-poisson"}, \code{"modele=pls-glm-family"})}{the method to be used in fitting the model. The default method \code{"glm.fit"} uses iteratively reweighted least squares (IWLS). User-supplied fitting functions can be supplied either as a function or a character string naming a function, with a function which takes the same arguments as \code{glm.fit}. If \code{"model.frame"}, the model frame is returned.}
  \item{\code{pls-glm-polr}}{logistic, probit, complementary log-log (cloglog) or cauchit (corresponding to a Cauchy latent variable).}}}
  \item{control}{a list of parameters for controlling the fitting process. For \code{glm.fit} this is passed to \code{\link{glm.control}}.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  \item{sparse}{should the coefficients of non-significant predictors (<\code{alpha.pvals.expli}) be set to 0}
  \item{sparseStop}{should component extraction stop when no significant predictors (<\code{alpha.pvals.expli}) are found}
  \item{naive}{Use the naive estimates for the Degrees of Freedom in plsR? Default is \code{FALSE}.}
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
\seealso{\code{\link[plsRglm]{PLS_glm_wvc}} and \code{\link[plsRglm]{PLS_glm_kfoldcv_formula}}}
\examples{
data(Cornell)
PLS_glm_formula(Y~.,data=Cornell,3)$uscores
PLS_glm_formula(Y~.,data=Cornell,3)$pp
PLS_glm_formula(Y~.,data=Cornell,3)$Coeffs
PLS_glm_formula(Y~.,data=Cornell,10)$InfCrit
PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-gaussian")$InfCrit
data.frame(pls=PLS_glm_formula(Y~.,data=Cornell,3)$Coeffs,PLS_glm_formula_gaussian=
PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-gaussian")$Coeffs,
PLS_glm_formula_family_gaussian=PLS_glm_formula(Y~.,data=Cornell,10,
modele="pls-glm-family",family=gaussian())$Coeffs)

mod <- PLS_glm_formula(Y~.,data=Cornell,10,pvals.expli =TRUE)
mod2 <- PLS_glm_formula(Y~.,data=Cornell,10,sparse=TRUE)
mod3 <- PLS_glm_formula(Y~.,data=Cornell,10,sparse=TRUE,sparseStop=FALSE)


\donttest{
## User specified links can be used.
## Example of user-specified link, a logit model for p^days
## See Shaffer, T.  2004. Auk 121(2): 526-540 and ?family.
logexp <- function(days = 1)
{
    linkfun <- function(mu) qlogis(mu^(1/days))
    linkinv <- function(eta) plogis(eta)^days
    mu.eta <- function(eta) days * plogis(eta)^(days-1) *
      .Call("logit_mu_eta", eta, PACKAGE = "stats")
    valideta <- function(eta) TRUE
    link <- paste("logexp(", days, ")", sep="")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
binomial(logexp(3))

data(aze_compl)
modpls <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",
family=binomial(link=logexp(3)),MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
rm(list=c("modpls","logexp","aze_compl"))


data(pine)
PLS_glm_formula(log(x11)~.,data=pine,1)$Std.Coeffs
PLS_glm_formula(log(x11)~.,data=pine,1)$Coeffs
PLS_glm_formula(log(x11)~.,data=pine,4)$Std.Coeffs
PLS_glm_formula(log(x11)~.,data=pine,4)$Coeffs
PLS_glm_formula(log(x11)~.,data=pine,4)$PredictY[1,]
PLS_glm_formula(log(x11)~.,data=pine,4,dataPredictY=pine[1,-11])$PredictY[1,]
PLS_glm_formula(log(x11)~.,data=pine,4,dataPredictY=pine[1,-11])$ValsPredictY[1]

pineNAX21 <- pine
pineNAX21[1,2] <- NA
str(PLS_glm_formula(log(x11)~.,data=pineNAX21,2))
PLS_glm_formula(log(x11)~.,data=pineNAX21,4)$Std.Coeffs
PLS_glm_formula(log(x11)~.,data=pineNAX21,4)$YChapeau[1,]
PLS_glm_formula(log(x11)~.,data=pine,4)$YChapeau[1,]
PLS_glm_formula(log(x11)~.,data=pineNAX21,4)$CoeffC
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,EstimXNA=TRUE)$XChapeau
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,EstimXNA=TRUE)$XChapeauNA

# compare pls-glm-gaussian with classic plsR
cbind(PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls")$Std.Coeffs,
PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-gaussian")$Std.Coeffs)

# without missing data
cbind(log(pine$x11),PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-gaussian")$YChapeau)
cbind(log(pine$x11),PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-gaussian")$YChapeau)

# with missing data
cbind((log(pine$x11)),PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-gaussian")$YChapeau)
cbind((log(pine$x11)),PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls")$ValsPredictY,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-gaussian")$ValsPredictY)


# compare pls-glm-gaussian with log link with classic plsR on the log
cbind(PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls")$Std.Coeffs,
PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-gaussian")$Std.Coeffs,
PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-family",
family=gaussian(link="identity"))$Std.Coeffs,PLS_glm_formula(x11~.,data=pine,4,
modele="pls-glm-family",family=gaussian(link=log))$Std.Coeffs)

# without missing data
cbind(log(pine$x11),PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-gaussian")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-family",
family=gaussian(link="identity"))$YChapeau,log(PLS_glm_formula(x11~.,data=pineNAX21,
4,modele="pls-glm-family",family=gaussian(link=log))$YChapeau))

# with missing data
cbind(log(pine$x11),PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-gaussian")$YChapeau,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-family",
family=gaussian(link="identity"))$YChapeau,log(PLS_glm_formula(x11~.,data=pineNAX21,4,
modele="pls-glm-family",family=gaussian(link=log))$YChapeau))
cbind(log(pine$x11),PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls")$ValsPredictY,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-gaussian")$ValsPredictY,
PLS_glm_formula(log(x11)~.,data=pineNAX21,4,modele="pls-glm-family",
family=gaussian(link="identity"))$ValsPredictY,log(PLS_glm_formula(x11~.,data=pineNAX21,4,
modele="pls-glm-family",family=gaussian(link=log))$ValsPredictY))


#other links
data.frame(pls=PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls")$Std.Coeffs,
identity=PLS_glm_formula(log(x11)~.,data=pine,4,modele="pls-glm-family",
family=gaussian(link="identity"))$Std.Coeffs,log=PLS_glm_formula(x11~.,data=pine,4,
modele="pls-glm-family",family=gaussian(link=log))$Std.Coeffs,
inverse=PLS_glm_formula(x11~.,data=pine,4,modele="pls-glm-family",
family=gaussian(link=inverse))$Std.Coeffs)



data(fowlkes)
modpls <- PLS_glm_formula(Y~.,data=fowlkes,4,modele="pls-glm-logistic",
pvals.expli=TRUE)
modpls$pvalstep
rm(list=c("modpls"))


data(aze_compl)
PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-logistic",
MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
modpls$valpvalstep
modpls$Coeffsmodel_vals
modpls2 <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",
family=binomial(link=logit),MClassed=TRUE,pvals.expli=TRUE)
modpls3 <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",
family=binomial(link=probit),MClassed=TRUE,pvals.expli=TRUE)
modpls4 <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",
family=binomial(link=cauchit),MClassed=TRUE,pvals.expli=TRUE)
#fails modpls5 <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",
family=binomial(link=log),MClassed=TRUE,pvals.expli=TRUE)
modpls6 <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",
family=binomial(link=cloglog),MClassed=TRUE,pvals.expli=TRUE)
data.frame(logit=modpls2$Std.Coeffs,probit=modpls3$Std.Coeffs,
cauchit=modpls4$Std.Coeffs,log=NA,cloglog=modpls6$Std.Coeffs)

plot(PLS_glm_formula(y~.,data=aze_compl,4,modele="pls-glm-logistic")$FinalModel)
PLS_glm_formula(y~.,data=aze_compl[-c(99,72),],4,modele="pls-glm-logistic",
pvals.expli=TRUE)$pvalstep
plot(PLS_glm_formula(y~.,data=aze_compl[-c(99,72),],4,
modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)

modpls7 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",
MClassed=TRUE,sparse=TRUE)
modpls8 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",
MClassed=TRUE,sparse=TRUE,sparseStop=FALSE)
modpls7$InfCrit
modpls7$valpvalstep
modpls7$Coeffs

modpls$InfCrit
modpls7$InfCrit
colSums(modpls$pvalstep)
colSums(modpls7$pvalstep)
rm(list=c("modpls"))


data(bordeaux)
modpls <- PLS_glm_formula(Quality~.,data=bordeaux,10,modele="pls-glm-polr")
modpls <- PLS_glm_formula(Quality~Temperature+Sunshine+Heat+Rain,
data=bordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit

bordeauxNA<-bordeaux
bordeauxNA[1,1] <- NA
modplsNA <- PLS_glm_formula(Quality~Temperature+Sunshine+Heat+Rain,
data=bordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit
rm(list=c("bordeauxNA"))

modpls2 <- PLS_glm_formula(Quality~.,data=bordeaux,10,modele="pls-glm-polr",
method="logistic")
modpls3 <- PLS_glm_formula(Quality~.,data=bordeaux,10,modele="pls-glm-polr",
method="probit")
modpls4 <- PLS_glm_formula(Quality~.,data=bordeaux,10,modele="pls-glm-polr",
method="cloglog")
modpls5 <- PLS_glm_formula(Quality~.,data=bordeaux,2,modele="pls-glm-polr",
method="cauchit")

cbind(modpls2$Coeffs,modpls3$Coeffs,modpls4$Coeffs,modpls5$Coeffs)
rbind(modpls2$InfCrit[,"Chi2_Pearson_Y"],modpls3$InfCrit[,"Chi2_Pearson_Y"],
modpls4$InfCrit[,"Chi2_Pearson_Y"],c(modpls5$InfCrit[,"Chi2_Pearson_Y"],NA,NA))

modpls6 <- PLS_glm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr",pvals.expli=TRUE)
modpls6$Coeffsmodel_vals
modpls6$InfCrit
modpls6$valpvalstep
modpls6$Coeffs

modpls7 <- PLS_glm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr",
alpha.pvals.expli=.15,sparse=TRUE)
modpls8 <- PLS_glm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr",
alpha.pvals.expli=.15,sparse=TRUE,sparseStop=FALSE)
modpls7$Coeffsmodel_vals
modpls7$InfCrit
modpls7$valpvalstep
modpls7$Coeffs

modpls6$InfCrit
modpls7$InfCrit
colSums(modpls6$pvalstep)
colSums(modpls7$pvalstep)
rm(list=c("modpls","modplsNA","modpls2","modpls3","modpls4","modpls5",
"modpls6","modpls7"))


# Test of other families and links on the same datasets.
data(pine)
modpls <- PLS_glm_formula(x11~.,data=pine,10,modele="pls")
modpls$computed_nt
modpls$InfCrit
modpls2 <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-Gamma")
modpls2$InfCrit
modpls2a <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-family",
family=Gamma(link=inverse))
modpls2a$InfCrit
modpls2b <- PLS_glm_formula(x11+1~.,data=pine,10,modele="pls-glm-family",
family=Gamma(link=identity))
modpls2b$InfCrit
modpls2c <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-family",
family=Gamma(link=log))
modpls2c$InfCrit


modpls3 <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-gaussian")
modpls3$InfCrit
modpls3a <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-family",
family=gaussian(link=identity))
modpls3a$InfCrit
modpls3b <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-family",
family=gaussian(link=log))
modpls3b$InfCrit
modpls3c <- PLS_glm_formula(x11~.,data=pine,10,modele="pls-glm-family",
family=gaussian(link=inverse))
modpls3c$InfCrit


modpls4 <- PLS_glm_formula(round(x11)~.,data=pine,10,modele="pls-glm-poisson")
modpls4$InfCrit
modpls4a <- PLS_glm_formula(round(x11)~.,data=pine,10,modele="pls-glm-family",
family=poisson(link=log))
modpls4a$InfCrit
modpls4b <- PLS_glm_formula(round(x11)+1~.,data=pine,10,modele="pls-glm-family",
family=poisson(link=identity))
modpls4b$InfCrit
modpls4c <- PLS_glm_formula(round(x11)+1~.,data=pine,10,modele="pls-glm-family",
family=poisson(link=sqrt))
modpls4c$InfCrit
rm(list=c("pine","modpls","modpls2","modpls3","modpls4","modpls2a","modpls3a",
"modpls4a","modpls2b","modpls3b","modpls4b","modpls2c","modpls3c","modpls4c"))

data(Cornell)
modpls <- PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-inverse.gaussian")
modpls$InfCrit
modplsa <- PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-family",
family=inverse.gaussian(link=1/mu^2))
modplsa$InfCrit
modplsb <- PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-family",
family=inverse.gaussian(link=inverse))
modplsb$InfCrit
modplsc <- PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-family",
family=inverse.gaussian(link=identity))
modplsc$InfCrit
modplsd <- PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-family",
family=inverse.gaussian(link=log))
modplsd$InfCrit
rm(list=c("modpls","modplsa","modplsb","modplsc","modplsd"))


dimX <- 6
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- PLS_glm_formula(ysimbin1~Xsimbin1,nt=10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 2
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar2)[,1]
Xsimbin1 <- dicho(dataAstar2)[,2:(dimX+1)]
modplsglm <- PLS_glm_formula(ysimbin1~Xsimbin1,nt=10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar2","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 3
dataAstar3 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar3)[,1]
Xsimbin1 <- dicho(dataAstar3)[,2:(dimX+1)]
modplsglm <- PLS_glm_formula(ysimbin1~Xsimbin1,nt=10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar3","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- PLS_glm_formula(ysimbin1~Xsimbin1,nt=10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 5
dataAstar5 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar5)[,1]
Xsimbin1 <- dicho(dataAstar5)[,2:(dimX+1)]
modplsglm <- PLS_glm_formula(ysimbin1~Xsimbin1,nt=10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar5","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 6
dataAstar6 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar6)[,1]
Xsimbin1 <- dicho(dataAstar6)[,2:(dimX+1)]
modplsglm <- PLS_glm_formula(ysimbin1~Xsimbin1,nt=10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar6","ysimbin1","Xsimbin1","modplsglm"))
}
}
\keyword{models}
\keyword{regression}
