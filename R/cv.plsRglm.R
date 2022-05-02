#' Partial least squares regression glm models with k-fold cross validation
#' 
#' This function implements k-fold cross-validation on complete or incomplete
#' datasets for partial least squares regression generalized linear models
#' 
#' Predicts 1 group with the \code{K-1} other groups. Leave one out cross
#' validation is thus obtained for \code{K==nrow(dataX)}.
#' 
#' There are seven different predefined models with predefined link functions
#' available : \describe{ \item{list("\"pls\"")}{ordinary pls models}
#' \item{list("\"pls-glm-Gamma\"")}{glm gaussian with inverse link pls models}
#' \item{list("\"pls-glm-gaussian\"")}{glm gaussian with identity link pls
#' models} \item{list("\"pls-glm-inverse-gamma\"")}{glm binomial with square
#' inverse link pls models} \item{list("\"pls-glm-logistic\"")}{glm binomial
#' with logit link pls models} \item{list("\"pls-glm-poisson\"")}{glm poisson
#' with log link pls models} \item{list("\"pls-glm-polr\"")}{glm polr with
#' logit link pls models} } Using the \code{"family="} option and setting
#' \code{"modele=pls-glm-family"} allows changing the family and link function
#' the same way as for the \code{\link[stats]{glm}} function. As a consequence
#' user-specified families can also be used.  \describe{ \item{The }{accepts
#' the links (as names) \code{identity}, \code{log} and
#' \code{inverse}.}\item{list("gaussian")}{accepts the links (as names)
#' \code{identity}, \code{log} and \code{inverse}.}\item{ family}{accepts the
#' links (as names) \code{identity}, \code{log} and \code{inverse}.} \item{The
#' }{accepts the links \code{logit}, \code{probit}, \code{cauchit},
#' (corresponding to logistic, normal and Cauchy CDFs respectively) \code{log}
#' and \code{cloglog} (complementary log-log).}\item{list("binomial")}{accepts
#' the links \code{logit}, \code{probit}, \code{cauchit}, (corresponding to
#' logistic, normal and Cauchy CDFs respectively) \code{log} and \code{cloglog}
#' (complementary log-log).}\item{ family}{accepts the links \code{logit},
#' \code{probit}, \code{cauchit}, (corresponding to logistic, normal and Cauchy
#' CDFs respectively) \code{log} and \code{cloglog} (complementary log-log).}
#' \item{The }{accepts the links \code{inverse}, \code{identity} and
#' \code{log}.}\item{list("Gamma")}{accepts the links \code{inverse},
#' \code{identity} and \code{log}.}\item{ family}{accepts the links
#' \code{inverse}, \code{identity} and \code{log}.} \item{The }{accepts the
#' links \code{log}, \code{identity}, and
#' \code{sqrt}.}\item{list("poisson")}{accepts the links \code{log},
#' \code{identity}, and \code{sqrt}.}\item{ family}{accepts the links
#' \code{log}, \code{identity}, and \code{sqrt}.} \item{The }{accepts the links
#' \code{1/mu^2}, \code{inverse}, \code{identity} and
#' \code{log}.}\item{list("inverse.gaussian")}{accepts the links \code{1/mu^2},
#' \code{inverse}, \code{identity} and \code{log}.}\item{ family}{accepts the
#' links \code{1/mu^2}, \code{inverse}, \code{identity} and \code{log}.}
#' \item{The }{accepts the links \code{logit}, \code{probit}, \code{cloglog},
#' \code{identity}, \code{inverse}, \code{log}, \code{1/mu^2} and
#' \code{sqrt}.}\item{list("quasi")}{accepts the links \code{logit},
#' \code{probit}, \code{cloglog}, \code{identity}, \code{inverse}, \code{log},
#' \code{1/mu^2} and \code{sqrt}.}\item{ family}{accepts the links
#' \code{logit}, \code{probit}, \code{cloglog}, \code{identity},
#' \code{inverse}, \code{log}, \code{1/mu^2} and \code{sqrt}.} \item{The
#' function }{can be used to create a power link
#' function.}\item{list("power")}{can be used to create a power link function.}
#' \item{list()}{arguments to pass to \code{cv.plsRglmmodel.default} or to
#' \code{cv.plsRglmmodel.formula}} }
#' 
#' A typical predictor has the form response ~ terms where response is the
#' (numeric) response vector and terms is a series of terms which specifies a
#' linear predictor for response. A terms specification of the form first +
#' second indicates all the terms in first together with all the terms in
#' second with any duplicates removed.
#' 
#' A specification of the form first:second indicates the the set of terms
#' obtained by taking the interactions of all terms in first with all terms in
#' second. The specification first*second indicates the cross of first and
#' second. This is the same as first + second + first:second.
#' 
#' The terms in the formula will be re-ordered so that main effects come first,
#' followed by the interactions, all second-order, all third-order and so on:
#' to avoid this pass a terms object as the formula.
#' 
#' Non-NULL weights can be used to indicate that different observations have
#' different dispersions (with the values in weights being inversely
#' proportional to the dispersions); or equivalently, when the elements of
#' weights are positive integers w_i, that each response y_i is the mean of w_i
#' unit-weight observations.
#' 
#' @aliases cv.plsRglm cv.plsRglmmodel.default cv.plsRglmmodel.formula
#' PLS_glm_kfoldcv PLS_glm_kfoldcv_formula
#' @param x a formula or a response (training) dataset
#' @param dataY response (training) dataset
#' @param dataX predictor(s) (training) dataset
#' @param formula an object of class "\code{\link{formula}}" (or one that can
#' be coerced to that class): a symbolic description of the model to be fitted.
#' The details of model specification are given under 'Details'.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in
#' the model. If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{plsRglm} is called.
#' @param nt number of components to be extracted
#' @param limQ2set limit value for the Q2
#' @param modele name of the PLS glm model to be fitted (\code{"pls"},
#' \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"},
#' \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"},
#' \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}). Use
#' \code{"modele=pls-glm-family"} to enable the \code{family} option.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link[stats]{family}} for details of family functions.) To use
#' the family option, please set \code{modele="pls-glm-family"}. User defined
#' families can also be defined. See details.
#' @param K number of groups. Defaults to 5.
#' @param NK number of times the group division is made
#' @param grouplist to specify the members of the \code{K} groups
#' @param random should the \code{K} groups be made randomly. Defaults to
#' \code{TRUE}
#' @param scaleX scale the predictor(s) : must be set to TRUE for
#' \code{modele="pls"} and should be for glms pls.
#' @param scaleY scale the response : Yes/No. Ignored since non always possible
#' for glm responses.
#' @param keepcoeffs shall the coefficients for each model be returned
#' @param keepfolds shall the groups' composition be returned
#' @param keepdataY shall the observed value of the response for each one of
#' the predicted value be returned
#' @param keepMclassed shall the number of miss classed be returned
#' (unavailable)
#' @param tol_Xi minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times
#' pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It
#' defaults to \eqn{10^{-12}}{10^{-12}}
#' @param weights an optional vector of 'prior weights' to be used in the
#' fitting process. Should be \code{NULL} or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param start starting values for the parameters in the linear predictor.
#' @param etastart starting values for the linear predictor.
#' @param mustart starting values for the vector of means.
#' @param offset this can be used to specify an \emph{a priori} known component
#' to be included in the linear predictor during fitting. This should be
#' \code{NULL} or a numeric vector of length equal to the number of cases. One
#' or more \code{\link{offset}} terms can be included in the formula instead or
#' as well, and if more than one is specified their sum is used. See
#' \code{\link{model.offset}}.
#' @param method \describe{ \item{for fitting glms with glm (}{the method to be
#' used in fitting the model. The default method \code{"glm.fit"} uses
#' iteratively reweighted least squares (IWLS). User-supplied fitting functions
#' can be supplied either as a function or a character string naming a
#' function, with a function which takes the same arguments as \code{glm.fit}.
#' If "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-Gamma\"")}{the method to be used in fitting
#' the model. The default method \code{"glm.fit"} uses iteratively reweighted
#' least squares (IWLS). User-supplied fitting functions can be supplied either
#' as a function or a character string naming a function, with a function which
#' takes the same arguments as \code{glm.fit}. If "model.frame", the model
#' frame is returned.}\item{, }{the method to be used in fitting the model. The
#' default method \code{"glm.fit"} uses iteratively reweighted least squares
#' (IWLS). User-supplied fitting functions can be supplied either as a function
#' or a character string naming a function, with a function which takes the
#' same arguments as \code{glm.fit}. If "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-gaussian\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-inverse.gaussian\"")}{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-logistic\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-poisson\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"modele=pls-glm-family\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{)}{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}
#' \item{list("pls-glm-polr")}{logistic, probit, complementary log-log or
#' cauchit (corresponding to a Cauchy latent variable).}}
#' @param control a list of parameters for controlling the fitting process. For
#' \code{glm.fit} this is passed to \code{\link{glm.control}}.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{model.matrix.default}.
#' @param verbose should info messages be displayed ?
#' @param \dots arguments to pass to \code{cv.plsRglmmodel.default} or to
#' \code{cv.plsRglmmodel.formula}
#' @return An object of class \code{"cv.plsRglmmodel"}.\cr
#' \item{results_kfolds}{list of \code{NK}. Each element of the list sums up
#' the results for a group division: \describe{ \item{list}{ of \code{K}
#' matrices of size about \code{nrow(dataX)/K * nt} with the predicted values
#' for a growing number of components} \item{list()}{\dots{}} \item{list}{ of
#' \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted
#' values for a growing number of components} }} \item{folds}{list of
#' \code{NK}. Each element of the list sums up the informations for a group
#' division: \describe{ \item{list}{ of \code{K} vectors of length about
#' \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were
#' used as a training set} \item{list()}{\dots{}} \item{list}{ of \code{K}
#' vectors of length about \code{nrow(dataX)} with the numbers of the rows of
#' \code{dataX} that were used as a training set} } } \item{dataY_kfolds}{list
#' of \code{NK}. Each element of the list sums up the results for a group
#' division: \describe{ \item{list}{ of \code{K} matrices of size about
#' \code{nrow(dataX)/K * 1} with the observed values of the response}
#' \item{list()}{\dots{}} \item{list}{ of \code{K} matrices of size about
#' \code{nrow(dataX)/K * 1} with the observed values of the response} } }
#' \item{call}{the call of the function}
#' @note Work for complete and incomplete datasets.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso Summary method \code{summary.cv.plsRglmmodel}.
#' \code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}},
#' \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}},
#' \code{\link{kfolds2Mclassed}} and \code{\link{summary}} to extract and
#' transform results from k-fold cross validation.
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' @keywords models regression
#' @examples
#' 
#' data(Cornell)
#' bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10)
#' (sum1<-summary(bbb))
#' cvtable(sum1)
#' 
#' bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=3,
#' modele="pls-glm-family",family=gaussian(),K=12,verbose=FALSE)
#' (sum2<-summary(bbb2))
#' cvtable(sum2)
#' 
#' \donttest{
#' #random=TRUE is the default to randomly create folds for repeated CV
#' bbb3 <- cv.plsRglm(Y~.,data=Cornell,nt=3,
#' modele="pls-glm-family",family=gaussian(),K=6,NK=10, verbose=FALSE)
#' (sum3<-summary(bbb3))
#' plot(cvtable(sum3))
#' 
#' data(aze_compl)
#' bbb <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,modele="pls",keepcoeffs=TRUE, verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,modele="pls-glm-family",
#' family=binomial(probit),keepcoeffs=TRUE, verbose=FALSE)
#' bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,
#' modele="pls-glm-logistic",keepcoeffs=TRUE, verbose=FALSE)
#' summary(bbb,MClassed=TRUE)
#' summary(bbb2,MClassed=TRUE)
#' kfolds2coeff(bbb2)
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' 
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,modele="pls-glm-family",
#' family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,
#' modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' boxplot(kfolds2coeff(bbb)[,1])
#' 
#' kfolds2Chisqind(bbb)
#' kfolds2Chisq(bbb)
#' summary(bbb)
#' PLS_lm(ypine,Xpine,10,typeVC="standard")$InfCrit
#' 
#' data(pineNAX21)
#' bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
#' modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
#' modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' 
#' data(XpineNAX21)
#' PLS_lm(ypine,XpineNAX21,10,typeVC="standard")$InfCrit
#' rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))
#' 
#' 
#' 
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",
#' family=Gamma,K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-Gamma",
#' K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' boxplot(kfolds2coeff(bbb)[,1])
#' 
#' kfolds2Chisqind(bbb)
#' kfolds2Chisq(bbb)
#' summary(bbb)
#' PLS_lm(ypine,Xpine,10,typeVC="standard")$InfCrit
#' 
#' data(pineNAX21)
#' bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
#' modele="pls-glm-family",family=Gamma(),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
#' modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' XpineNAX21 <- Xpine
#' XpineNAX21[1,2] <- NA
#' PLS_lm(ypine,XpineNAX21,10,typeVC="standard")$InfCrit
#' rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))
#' 
#' 
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,modele="pls",verbose=FALSE)
#' summary(bbb)
#' 
#' cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=12,verbose=FALSE)
#' cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
#' family=inverse.gaussian,K=12,verbose=FALSE)
#' cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
#' NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=inverse.gaussian(),
#' K=6,NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
#' NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
#' family=inverse.gaussian(link = "1/mu^2"),K=6,NK=2,verbose=FALSE)$results_kfolds
#' 
#' bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=10,
#' modele="pls-glm-inverse.gaussian",keepcoeffs=TRUE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' PLS_lm(yCornell,XCornell,10,typeVC="standard")$InfCrit
#' rm(list=c("XCornell","yCornell","bbb","bbb2"))
#' }
#' data(Cornell)
#' bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,modele="pls")
#' summary(bbb)
#' 
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=12)
#' 
#' \donttest{
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,
#' NK=2,random=TRUE,keepfolds=TRUE,verbose=FALSE)$results_kfolds
#' 
#' #Different ways of model specifications
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,
#' NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian,
#' K=6,NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),
#' K=6,NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),
#' K=6,NK=2,verbose=FALSE)$results_kfolds
#' 
#' bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=10,
#' modele="pls-glm-gaussian",keepcoeffs=TRUE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",
#' family=gaussian(link=log),K=6,keepcoeffs=TRUE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' PLS_lm_formula(Y~.,data=Cornell,10,typeVC="standard")$InfCrit
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' data(pine)
#' bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",
#' family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",family=gaussian(),
#' K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' boxplot(kfolds2coeff(bbb)[,1])
#' 
#' kfolds2Chisqind(bbb)
#' kfolds2Chisq(bbb)
#' summary(bbb)
#' PLS_lm_formula(x11~.,data=pine,nt=10,typeVC="standard")$InfCrit
#' 
#' data(pineNAX21)
#' bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
#' modele="pls-glm-family",family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
#' modele="pls-glm-gaussian",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' PLS_lm_formula(x11~.,data=pineNAX21,nt=10,typeVC="standard")$InfCrit
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' data(aze_compl)
#' bbb <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,modele="pls",
#' keepcoeffs=TRUE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=3,K=10,
#' modele="pls-glm-family",family=binomial(probit),keepcoeffs=TRUE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=3,K=10,
#' modele="pls-glm-logistic",keepcoeffs=TRUE,verbose=FALSE)
#' summary(bbb,MClassed=TRUE)
#' summary(bbb2,MClassed=TRUE)
#' kfolds2coeff(bbb2)
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' 
#' data(pine)
#' bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,
#' modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,
#' modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' boxplot(kfolds2coeff(bbb)[,1])
#' 
#' kfolds2Chisqind(bbb)
#' kfolds2Chisq(bbb)
#' summary(bbb)
#' PLS_lm_formula(x11~.,data=pine,10,typeVC="standard")$InfCrit
#' 
#' data(pineNAX21)
#' bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
#' modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
#' modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' PLS_lm_formula(x11~.,data=pineNAX21,10,typeVC="standard")$InfCrit
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' 
#' data(pine)
#' bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",
#' family=Gamma,K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-Gamma",
#' K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb)
#' boxplot(kfolds2coeff(bbb)[,1])
#' 
#' kfolds2Chisqind(bbb)
#' kfolds2Chisq(bbb)
#' summary(bbb)
#' PLS_lm_formula(x11~.,data=pine,10,typeVC="standard")$InfCrit
#' 
#' data(pineNAX21)
#' bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
#' modele="pls-glm-family",family=Gamma(),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
#' modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' PLS_lm_formula(x11~.,data=pineNAX21,10,typeVC="standard")$InfCrit
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' 
#' data(Cornell)
#' summary(cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,modele="pls",verbose=FALSE))
#' 
#' cv.plsRglm(Y~.,data=Cornell,nt=3,
#' modele="pls-glm-inverse.gaussian",K=12,verbose=FALSE)
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=inverse.gaussian,K=12,verbose=FALSE)
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
#' NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",
#' family=inverse.gaussian(),K=6,NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
#' NK=2,verbose=FALSE)$results_kfolds
#' cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",
#' family=inverse.gaussian(link = "1/mu^2"),K=6,NK=2,verbose=FALSE)$results_kfolds
#' 
#' bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=10,
#' modele="pls-glm-inverse.gaussian",keepcoeffs=TRUE,verbose=FALSE)
#' 
#' #For Jackknife computations
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' 
#' kfolds2Chisqind(bbb2)
#' kfolds2Chisq(bbb2)
#' summary(bbb2)
#' PLS_lm_formula(Y~.,data=Cornell,10,typeVC="standard")$InfCrit
#' rm(list=c("bbb","bbb2"))
#' 
#' 
#' data(bordeaux)
#' summary(cv.plsRglm(Quality~.,data=bordeaux,10,
#' modele="pls-glm-polr",K=7))
#' 
#' data(bordeauxNA)
#' summary(cv.plsRglm(Quality~.,data=bordeauxNA,
#' 10,modele="pls-glm-polr",K=10,verbose=FALSE))
#' 
#' summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
#' modele="pls-glm-polr",method="logistic",verbose=FALSE))
#' summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
#' modele="pls-glm-polr",method="probit",verbose=FALSE))
#' summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
#' modele="pls-glm-polr",method="cloglog",verbose=FALSE))
#' suppressWarnings(summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
#' modele="pls-glm-polr",method="cauchit",verbose=FALSE)))
#' }
#' 
#' @export cv.plsRglm
cv.plsRglm <- function (x, ...) {UseMethod("cv.plsRglmmodel")}

#' @rdname cv.plsRglm
#' @aliases cv.plsRglm
#' @export cv.plsRglmmodel
cv.plsRglmmodel <- cv.plsRglm
