#' Partial least squares Regression generalized linear models
#' 
#' This function implements Partial least squares Regression generalized linear
#' models complete or incomplete datasets.
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
#' }
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
#' The default estimator for Degrees of Freedom is the Kramer and Sugiyama's
#' one which only works for classical plsR models. For these models,
#' Information criteria are computed accordingly to these estimations. Naive
#' Degrees of Freedom and Information Criteria are also provided for comparison
#' purposes. For more details, see N. Kraemer and M. Sugiyama. (2011). The
#' Degrees of Freedom of Partial Least Squares Regression. \emph{Journal of the
#' American Statistical Association}, 106(494), 697-705, 2011.
#' 
#' @aliases plsRglm plsRglmmodel.default plsRglmmodel.formula PLS_glm
#' PLS_glm_formula
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
#' @param dataPredictY predictor(s) (testing) dataset
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
#' @param typeVC type of leave one out cross validation. For back compatibility
#' purpose.  \describe{ \item{list("none")}{no cross validation} }
#' @param EstimXNA only for \code{modele="pls"}. Set whether the missing X
#' values have to be estimated.
#' @param scaleX scale the predictor(s) : must be set to TRUE for
#' \code{modele="pls"} and should be for glms pls.
#' @param scaleY scale the response : Yes/No. Ignored since non always possible
#' for glm responses.
#' @param pvals.expli should individual p-values be reported to tune model
#' selection ?
#' @param alpha.pvals.expli level of significance for predictors when
#' pvals.expli=TRUE
#' @param MClassed number of missclassified cases, should only be used for
#' binary responses
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
#' @param method For a glm model (\code{modele="pls-glm-family"}), the method
#' to be used in fitting the model. The default method \code{"glm.fit"} uses
#' iteratively reweighted least squares (IWLS). User-supplied fitting functions
#' can be supplied either as a function or a character string naming a
#' function, with a function which takes the same arguments as \code{glm.fit}.
#' For a polr model (\code{modele="pls-glm-polr"}), \code{logistic} or
#' \code{probit} or (complementary) log-log (\code{loglog} or \code{cloglog})
#' or \code{cauchit} (corresponding to a Cauchy latent variable).
#' @param control a list of parameters for controlling the fitting process. For
#' \code{glm.fit} this is passed to \code{\link{glm.control}}.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{model.matrix.default}.
#' @param sparse should the coefficients of non-significant predictors
#' (<\code{alpha.pvals.expli}) be set to 0
#' @param sparseStop should component extraction stop when no significant
#' predictors (<\code{alpha.pvals.expli}) are found
#' @param naive Use the naive estimates for the Degrees of Freedom in plsR?
#' Default is \code{FALSE}.
#' @param verbose Should details be displayed ?
#' @param \dots arguments to pass to \code{plsRmodel.default} or to
#' \code{plsRmodel.formula}
#' @return Depends on the model that was used to fit the model. You can
#' generally at least find these items.\cr \item{nr}{Number of observations}
#' \item{nc}{Number of predictors} \item{nt}{Number of requested components}
#' \item{ww}{raw weights (before L2-normalization)} \item{wwnorm}{L2 normed
#' weights (to be used with deflated matrices of predictor variables)}
#' \item{wwetoile}{modified weights (to be used with original matrix of
#' predictor variables)} \item{tt}{PLS components} \item{pp}{loadings of the
#' predictor variables} \item{CoeffC}{coefficients of the PLS components}
#' \item{uscores}{scores of the response variable} \item{YChapeau}{predicted
#' response values for the dataX set} \item{residYChapeau}{residuals of the
#' deflated response on the standardized scale} \item{RepY}{scaled response
#' vector} \item{na.miss.Y}{is there any NA value in the response vector}
#' \item{YNA}{indicatrix vector of missing values in RepY}
#' \item{residY}{deflated scaled response vector} \item{ExpliX}{scaled matrix
#' of predictors} \item{na.miss.X}{is there any NA value in the predictor
#' matrix} \item{XXNA}{indicator of non-NA values in the predictor matrix}
#' \item{residXX}{deflated predictor matrix} \item{PredictY}{response values
#' with NA replaced with 0} \item{RSS}{residual sum of squares (original
#' scale)} \item{RSSresidY}{residual sum of squares (scaled scale)}
#' \item{R2residY}{R2 coefficient value on the standardized scale} \item{R2}{R2
#' coefficient value on the original scale} \item{press.ind}{individual PRESS
#' value for each observation (scaled scale)} \item{press.tot}{total PRESS
#' value for all observations (scaled scale)} \item{Q2cum}{cumulated Q2
#' (standardized scale)} \item{family}{glm family used to fit PLSGLR model}
#' \item{ttPredictY}{PLS components for the dataset on which prediction was
#' requested} \item{typeVC}{type of leave one out cross-validation used}
#' \item{dataX}{predictor values} \item{dataY}{response values}
#' \item{weights}{weights of the observations} \item{computed_nt}{number of
#' components that were computed} \item{AIC}{AIC vs number of components}
#' \item{BIC}{BIC vs number of components} \item{Coeffsmodel_vals}{}
#' \item{ChisqPearson}{} \item{CoeffCFull}{matrix of the coefficients of the
#' predictors} \item{CoeffConstante}{value of the intercept (scaled scale)}
#' \item{Std.Coeffs}{Vector of standardized regression coefficients}
#' \item{Coeffs}{Vector of regression coefficients (used with the original data
#' scale)} \item{Yresidus}{residuals of the PLS model}
#' \item{residusY}{residuals of the deflated response on the standardized
#' scale} \item{InfCrit}{table of Information Criteria:\cr \describe{
#' \item{list("AIC")}{AIC vs number of components} \item{list("BIC")}{BIC vs
#' number of components} \item{list("MissClassed")}{Number of miss classed
#' results} \item{list("Chi2_Pearson_Y")}{Q2 value (standardized scale)}
#' \item{list("RSS")}{residual sum of squares (original scale)}
#' \item{list("R2")}{R2 coefficient value on the original scale}
#' \item{list("R2residY")}{R2 coefficient value on the standardized scale}
#' \item{list("RSSresidY")}{residual sum of squares (scaled scale)} } }
#' \item{Std.ValsPredictY}{predicted response values for supplementary dataset
#' (standardized scale)} \item{ValsPredictY}{predicted response values for
#' supplementary dataset (original scale)} \item{Std.XChapeau}{estimated values
#' for missing values in the predictor matrix (standardized scale)}
#' \item{FinalModel}{final GLR model on the PLS components}
#' \item{XXwotNA}{predictor matrix with missing values replaced with 0}
#' \item{call}{call}
#' 
#' \item{AIC.std}{AIC.std vs number of components (AIC computed for the
#' standardized model}
#' @note Use \code{\link{cv.plsRglm}} to cross-validate the plsRglm models and
#' \code{\link{bootplsglm}} to bootstrap them.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{plsR}}.
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparaison de la
#' régression PLS et de la régression
#' logistique PLS : application aux données
#' d'allélotypage. \emph{Journal de la Societe Francaise
#' de Statistique}, 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords models regression
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' 
#' modplsglm <- plsRglm(yCornell,XCornell,10,modele="pls-glm-gaussian")
#' 
#' #To retrieve the final GLR model on the PLS components
#' finalmod <- modplsglm$FinalModel
#' #It is a glm object.
#' plot(finalmod)
#' 
#' \donttest{
#' #Cross validation
#' cv.modplsglm<-cv.plsRglm(Y~.,data=Cornell,6,NK=100,modele="pls-glm-gaussian", verbose=FALSE)
#' res.cv.modplsglm<-cvtable(summary(cv.modplsglm))
#' plot(res.cv.modplsglm)
#' 
#' #If no model specified, classic PLSR model
#' modpls <- plsRglm(Y~.,data=Cornell,6)
#' modpls
#' modpls$tt
#' modpls$uscores
#' modpls$pp
#' modpls$Coeffs
#' 
#' #rm(list=c("XCornell","yCornell",modpls,cv.modplsglm,res.cv.modplsglm))
#' }
#' 
#' data(aze_compl)
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE, verbose=FALSE)$InfCrit
#' modpls <- plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",
#' MClassed=TRUE,pvals.expli=TRUE, verbose=FALSE)
#' modpls
#' colSums(modpls$pvalstep)
#' modpls$Coeffsmodel_vals
#' 
#' plot(plsRglm(yaze_compl,Xaze_compl,4,modele="pls-glm-logistic")$FinalModel)
#' plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,
#' modele="pls-glm-logistic",pvals.expli=TRUE)$pvalstep
#' plot(plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,
#' modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)
#' rm(list=c("Xaze_compl","yaze_compl","modpls"))
#' 
#' 
#' data(bordeaux)
#' Xbordeaux<-bordeaux[,1:4]
#' ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
#' modpls <- plsRglm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr",pvals.expli=TRUE)
#' modpls
#' colSums(modpls$pvalstep)
#' 
#' 
#' XbordeauxNA<-Xbordeaux
#' XbordeauxNA[1,1] <- NA
#' modplsNA <- plsRglm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr",pvals.expli=TRUE)
#' modpls
#' colSums(modpls$pvalstep)
#' rm(list=c("Xbordeaux","XbordeauxNA","ybordeaux","modplsNA"))
#' 
#' \donttest{
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' modpls1 <- plsRglm(ypine,Xpine,1)
#' modpls1$Std.Coeffs
#' modpls1$Coeffs
#' modpls4 <- plsRglm(ypine,Xpine,4)
#' modpls4$Std.Coeffs
#' modpls4$Coeffs
#' modpls4$PredictY[1,]
#' plsRglm(ypine,Xpine,4,dataPredictY=Xpine[1,])$PredictY[1,]
#' 
#' XpineNAX21 <- Xpine
#' XpineNAX21[1,2] <- NA
#' modpls4NA <- plsRglm(ypine,XpineNAX21,4)
#' modpls4NA$Std.Coeffs
#' modpls4NA$YChapeau[1,]
#' modpls4$YChapeau[1,]
#' modpls4NA$CoeffC
#' plsRglm(ypine,XpineNAX21,4,EstimXNA=TRUE)$XChapeau
#' plsRglm(ypine,XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA
#' 
#' # compare pls-glm-gaussian with classic plsR
#' modplsglm4 <- plsRglm(ypine,Xpine,4,modele="pls-glm-gaussian")
#' cbind(modpls4$Std.Coeffs,modplsglm4$Std.Coeffs)
#' 
#' # without missing data
#' cbind(ypine,modpls4$ValsPredictY,modplsglm4$ValsPredictY)
#' 
#' # with missing data
#' modplsglm4NA <- plsRglm(ypine,XpineNAX21,4,modele="pls-glm-gaussian")
#' cbind((ypine),modpls4NA$ValsPredictY,modplsglm4NA$ValsPredictY)
#' rm(list=c("Xpine","ypine","modpls4","modpls4NA","modplsglm4","modplsglm4NA"))
#' 
#' data(fowlkes)
#' Xfowlkes <- fowlkes[,2:13]
#' yfowlkes <- fowlkes[,1]
#' modpls <- plsRglm(yfowlkes,Xfowlkes,4,modele="pls-glm-logistic",pvals.expli=TRUE)
#' modpls
#' colSums(modpls$pvalstep)
#' rm(list=c("Xfowlkes","yfowlkes","modpls"))
#' 
#' 
#' if(require(chemometrics)){
#' data(hyptis)
#' yhyptis <- factor(hyptis$Group,ordered=TRUE)
#' Xhyptis <- as.data.frame(hyptis[,c(1:6)])
#' options(contrasts = c("contr.treatment", "contr.poly"))
#' modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
#' modpls2$Coeffsmodel_vals
#' modpls2$InfCrit
#' modpls2$Coeffs
#' modpls2$Std.Coeffs
#' 
#' table(yhyptis,predict(modpls2$FinalModel,type="class"))
#' rm(list=c("yhyptis","Xhyptis","modpls2"))
#' }
#' 
#' dimX <- 24
#' Astar <- 6
#' dataAstar6 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
#' ysimbin1 <- dicho(dataAstar6)[,1]
#' Xsimbin1 <- dicho(dataAstar6)[,2:(dimX+1)]
#' modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
#' modplsglm
#' 
#' simbin=data.frame(dicho(dataAstar6))
#' cv.modplsglm <- suppressWarnings(cv.plsRglm(Y~.,data=simbin,nt=10,
#' modele="pls-glm-logistic",NK=100, verbose=FALSE))
#' res.cv.modplsglm <- cvtable(summary(cv.modplsglm,MClassed=TRUE,
#' verbose=FALSE))
#' plot(res.cv.modplsglm) #defaults to type="CVMC"
#' 
#' rm(list=c("dimX","Astar","dataAstar6","ysimbin1","Xsimbin1","modplsglm","cv.modplsglm",
#' "res.cv.modplsglm"))
#' }
#' 
#' @export plsRglm
plsRglm <- function(x, ...) UseMethod("plsRglmmodel")

#' @rdname plsRglm
#' @aliases plsRglm
#' @export plsRglmmodel
plsRglmmodel <- plsRglm
