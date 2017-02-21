\name{plsR}
\alias{plsR}
\alias{plsRmodel.default}
\alias{plsRmodel.formula}
\alias{PLS_lm}
\alias{PLS_lm_formula}
\title{Partial least squares Regression models with leave one out cross validation}
\description{
This function implements Partial least squares Regression models with leave one out cross validation for complete or incomplete datasets.
}
\usage{
plsR(x, \dots)
\method{plsRmodel}{default}(dataY, dataX, nt = 2, limQ2set = 0.0975, 
dataPredictY = dataX, modele = "pls", family = NULL, typeVC = "none", 
EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12), weights,
sparse = FALSE, sparseStop = TRUE, naive = FALSE,verbose=TRUE)
\method{plsRmodel}{formula}(formula, data, nt = 2, limQ2set = 0.0975,
dataPredictY, modele = "pls", family = NULL, typeVC = "none",
EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12), weights,
subset, contrasts = NULL, sparse = FALSE, sparseStop = TRUE, naive = FALSE,
verbose=TRUE)
PLS_lm(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, 
modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, 
scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12),
weights,sparse=FALSE,sparseStop=FALSE,naive=FALSE,verbose=TRUE)
PLS_lm_formula(formula,data=NULL,nt=2,limQ2set=.0975,dataPredictY=dataX,
modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,
scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,
tol_Xi=10^(-12),weights,subset,contrasts=NULL,sparse=FALSE,
sparseStop=FALSE,naive=FALSE,verbose=TRUE)
}
\arguments{
  \item{x}{a formula or a response (training) dataset}
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{formula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{plsR} is called.}  
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{dataPredictY}{predictor(s) (testing) dataset}
  \item{modele}{name of the PLS model to be fitted, only (\code{"pls"} available for this fonction.}
  \item{family}{for the present moment the family argument is ignored and set thanks to the value of modele.}
  \item{typeVC}{type of leave one out cross validation. Several procedures are available. If cross validation is required, one needs to selects the way of predicting the response for left out observations. For complete rows, without any missing value, there are two different ways of computing these predictions. As a consequence, for mixed datasets, with complete and incomplete rows, there are two ways of computing prediction : either predicts any row as if there were missing values in it (\code{missingdata}) or selects the prediction method accordingly to the completeness of the row (\code{adaptative}).
    \describe{
      \item{\code{none}}{no cross validation}
      \item{\code{standard}}{as in SIMCA for datasets without any missing value. For datasets with any missing value, it is the as using \code{missingdata}}
      \item{\code{missingdata}}{all values predicted as those with missing values for datasets with any missing values}
      \item{\code{adaptative}}{predict a response value for an x with any missing value as those with missing values and for an x without any missing value as those without missing values.}
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
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  \item{sparse}{should the coefficients of non-significant predictors (<\code{alpha.pvals.expli}) be set to 0}
  \item{sparseStop}{should component extraction stop when no significant predictors (<\code{alpha.pvals.expli}) are found}
  \item{naive}{Use the naive estimates for the Degrees of Freedom in plsR? Default is \code{FALSE}.}
  \item{verbose}{should info messages be displayed ?}
  \item{\dots}{arguments to pass to \code{plsRmodel.default} or to \code{plsRmodel.formula}}
}
\details{
There are several ways to deal with missing values that leads to different computations of leave one out cross validation criteria.

A typical predictor has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with any duplicates removed. 

A specification of the form first:second indicates the the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second. 

The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a terms object as the formula. 

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 

The default estimator for Degrees of Freedom is the Kramer and Sugiyama's one. Information criteria are computed accordingly to these estimations. Naive Degrees of Freedom and Information Criteria are also provided for comparison purposes. For more details, see N. Kraemer and M. Sugiyama. (2011). The Degrees of Freedom of Partial Least Squares Regression. \emph{Journal of the American Statistical Association}, 106(494), 697-705, 2011. 
}
\value{
  \item{nr}{Number of observations}
  \item{nc}{Number of predictors}
  \item{nt}{Number of requested components}
  \item{ww}{raw weights (before L2-normalization)}
  \item{wwnorm}{L2 normed weights (to be used with deflated matrices of predictor variables)}
  \item{wwetoile}{modified weights (to be used with original matrix of predictor variables)}
  \item{tt}{PLS components}
  \item{pp}{loadings of the predictor variables}
  \item{CoeffC}{coefficients of the PLS components}
  \item{uscores}{scores of the response variable}
  \item{YChapeau}{predicted response values for the dataX set}
  \item{residYChapeau}{residuals of the deflated response on the standardized scale}
  \item{RepY}{scaled response vector}
  \item{na.miss.Y}{is there any NA value in the response vector}
  \item{YNA}{indicatrix vector of missing values in RepY}
  \item{residY}{deflated scaled response vector}
  \item{ExpliX}{scaled matrix of predictors}
  \item{na.miss.X}{is there any NA value in the predictor matrix}
  \item{XXNA}{indicator of non-NA values in the predictor matrix}
  \item{residXX}{deflated predictor matrix}
  \item{PredictY}{response values with NA replaced with 0}
  \item{press.ind}{individual PRESS value for each observation (scaled scale)}
  \item{press.tot}{total PRESS value for all observations (scaled scale)}
  \item{family}{glm family used to fit PLSGLR model}
  \item{ttPredictY}{PLS components for the dataset on which prediction was requested}
  \item{typeVC}{type of leave one out cross-validation used}
  \item{dataX}{predictor values}
  \item{dataY}{response values}
  \item{computed_nt}{number of components that were computed}
  \item{CoeffCFull}{matrix of the coefficients of the predictors}
  \item{CoeffConstante}{value of the intercept (scaled scale)}
  \item{Std.Coeffs}{Vector of standardized regression coefficients}
  \item{press.ind2}{individual PRESS value for each observation (original scale)}
  \item{RSSresidY}{residual sum of squares (scaled scale)}
  \item{Coeffs}{Vector of regression coefficients (used with the original data scale)}
  \item{Yresidus}{residuals of the PLS model}
  \item{RSS}{residual sum of squares (original scale)}
  \item{residusY}{residuals of the deflated response on the standardized scale}
  \item{AIC.std}{AIC.std vs number of components (AIC computed for the standardized model}
  \item{AIC}{AIC vs number of components}
  \item{optional}{If the response is assumed to be binary:\cr
        i.e. \code{MClassed=TRUE}.    
  \describe{
  \item{\code{MissClassed}}{Number of miss classed results}
  \item{\code{Probs}}{"Probability" predicted by the model. These are not true probabilities since they may lay outside of [0,1]}
  \item{\code{Probs.trc}}{Probability predicted by the model and constrained to belong to [0,1]}
  }
  }
  \item{ttPredictFittedMissingY}{Description of 'comp2'}
  \item{optional}{If cross validation was requested:\cr
        i.e. \code{typeVC="standard"}, \code{typeVC="missingdata"} or \code{typeVC="adaptative"}.
  \describe{
  \item{\code{R2residY}}{R2 coefficient value on the standardized scale}
  \item{\code{R2}}{R2 coefficient value on the original scale}
  \item{\code{press.tot2}}{total PRESS value for all observations (original scale)}
  \item{\code{Q2}}{Q2 value (standardized scale)}
  \item{\code{limQ2}}{limit of the Q2 value}
  \item{\code{Q2_2}}{Q2 value (original scale)}
  \item{\code{Q2cum}}{cumulated Q2 (standardized scale)}
  \item{\code{Q2cum_2}}{cumulated Q2 (original scale)}
  }
  }
  \item{InfCrit}{table of Information Criteria}
  \item{Std.ValsPredictY}{predicted response values for supplementary dataset (standardized scale)}
  \item{ValsPredictY}{predicted response values for supplementary dataset (original scale)}
  \item{Std.XChapeau}{estimated values for missing values in the predictor matrix (standardized scale)}
  \item{XXwotNA}{predictor matrix with missing values replaced with 0}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{cv.plsR}} to cross-validate the plsRglm models and \code{\link{bootpls}} to bootstrap them.}
\seealso{See also \code{\link{plsRglm}} to fit PLSGLR models.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

#maximum 6 components could be extracted from this dataset
#trying 10 to trigger automatic stopping criterion
modpls10<-plsR(yCornell,XCornell,10)
modpls10

#With iterative leave one out CV PRESS
modpls6cv<-plsR(yCornell,XCornell,6,typeVC="standard")
modpls6cv
cv.modpls<-cv.plsR(yCornell,XCornell,6,NK=100)
res.cv.modpls<-cvtable(summary(cv.modpls))
plot(res.cv.modpls)

rm(list=c("XCornell","yCornell","modpls10","modpls6cv"))

\donttest{
#A binary response example
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
modpls.aze <- plsR(yaze_compl,Xaze_compl,10,MClassed=TRUE,typeVC="standard")
modpls.aze

#Direct access to not cross validated values
modpls.aze$AIC
modpls.aze$AIC.std
modpls.aze$MissClassed

#Raw predicted values (not really probabily since not constrained in [0,1]
modpls.aze$Probs
#Truncated to [0;1] predicted values (true probabilities)
modpls.aze$Probs.trc
modpls.aze$Probs-modpls.aze$Probs.trc

#Repeated cross validation of the model (NK=100 times)
cv.modpls.aze<-cv.plsR(yaze_compl,Xaze_compl,10,NK=100)
res.cv.modpls.aze<-cvtable(summary(cv.modpls.aze,MClassed=TRUE))
#High discrepancy in the number of component choice using repeated cross validation
#and missclassed criterion
plot(res.cv.modpls.aze)

rm(list=c("Xaze_compl","yaze_compl","modpls.aze","cv.modpls.aze","res.cv.modpls.aze"))

#24 predictors
dimX <- 24
#2 components
Astar <- 2
simul_data_UniYX(dimX,Astar)
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar2 <- dataAstar2[,1]
XdataAstar2 <- dataAstar2[,2:(dimX+1)]
modpls.A2<- plsR(ydataAstar2,XdataAstar2,10,typeVC="standard")
modpls.A2
cv.modpls.A2<-cv.plsR(ydataAstar2,XdataAstar2,10,NK=100)
res.cv.modpls.A2<-cvtable(summary(cv.modpls.A2))
#Perfect choice for the Q2 criterion in PLSR
plot(res.cv.modpls.A2)

#Binarized response
ysimbin1 <- dicho(ydataAstar2)
#Binarized predictors
Xsimbin1 <- dicho(XdataAstar2)
modpls.B2 <- plsR(ysimbin1,Xsimbin1,10,typeVC="standard",MClassed=TRUE)
modpls.B2
modpls.B2$Probs
modpls.B2$Probs.trc
modpls.B2$MissClassed
plsR(ysimbin1,XdataAstar2,10,typeVC="standard",MClassed=TRUE)$InfCrit
cv.modpls.B2<-cv.plsR(ysimbin1,Xsimbin1,2,NK=100)
res.cv.modpls.B2<-cvtable(summary(cv.modpls.B2,MClassed=TRUE))
#Only one component found by repeated CV missclassed criterion
plot(res.cv.modpls.B2)

rm(list=c("dimX","Astar","dataAstar2","XdataAstar2","ydataAstar2","modpls.A2","cv.modpls.A2",
"res.cv.modpls.A2","Xsimbin1","ysimbin1","modpls.B2","cv.modpls.B2","res.cv.modpls.B2"))
}
}
\keyword{models}
\keyword{regression}
