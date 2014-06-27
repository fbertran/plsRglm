\name{predict.plsRglmmodel}
\alias{predict.plsRglmmodel}
\title{Print method for plsRcox models}
\description{
This function provides a predict method for the class \code{"plsRglmmodel"}
}
\usage{
\method{predict}{plsRglmmodel}(object,newdata,comps=object$computed_nt,
type=c("link", "response", "terms", "scores", "class", "probs"),
se.fit=FALSE,weights, dispersion = NULL,methodNA="adaptative",...)
}
\arguments{
  \item{object}{An object of the class \code{"plsRmodel"}.}
  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{comps}{A value with a single value of component to use for prediction.}
  \item{type}{Type of predicted value. Available choices are the glms ones ("\code{link}", "\code{response}", "\code{terms}"), the polr ones ("\code{class}", "\code{probs}") or the scores ("\code{scores}").}
  \item{se.fit}{If TRUE, pointwise standard errors are produced for the predictions using the Cox model.}
  \item{weights}{Vector of case weights. If \code{weights} is a vector of integers, then the estimated coefficients are equivalent to estimating the model from data with the individual \code{cases} replicated as many times as indicated by \code{weights}.}
  \item{dispersion}{the dispersion of the GLM fit to be assumed in computing the standard errors. If omitted, that returned by summary applied to the object is used.}
  \item{methodNA}{Selects the way of predicting the response or the scores of the new data. For complete rows, without any missing value, there are two different ways of computing the prediction. As a consequence, for mixed datasets, with complete and incomplete rows, there are two ways of computing prediction : either predicts any row as if there were missing values in it (\code{missingdata}) or selects the prediction method accordingly to the completeness of the row (\code{adaptative}).}
  \item{\dots}{Arguments to be passed on to \code{stats::glm} and \code{plsRglm::plsRglm}.}
}
\value{
When type is "\code{response}", a matrix of predicted response values is returned.\cr
When type is "\code{scores}", a score matrix is returned. 
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{See Also \code{\link[stats]{predict.glm}}}
\examples{
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
data(pine_sup)
Xpine_sup<-pine_sup[,1:10]
Xpine_supNA<-Xpine_sup
Xpine_supNA[1,1]<-NA

modpls=plsRglm(dataY=ypine,dataX=Xpine,nt=6,modele="pls-glm-family",family="gaussian")
modplsform=plsRglm(x11~.,data=pine,nt=6,modele="pls-glm-family",family="gaussian")
modpls2=plsRglm(dataY=ypine,dataX=Xpine,nt=6,modele="pls-glm-family",
dataPredictY=Xpine_sup,family="gaussian")
modpls2NA=plsRglm(dataY=ypine,dataX=Xpine,nt=6,modele="pls-glm-family",
dataPredictY=Xpine_supNA,family="gaussian")

#Identical to predict(modpls,type="link") or modpls$Std.ValsPredictY
cbind(modpls$Std.ValsPredictY,modplsform$Std.ValsPredictY,
predict(modpls),predict(modplsform))

#Identical to predict(modpls,type="response") or modpls$ValsPredictY
cbind(modpls$ValsPredictY,modplsform$ValsPredictY,
predict(modpls,type="response"),predict(modplsform,type="response"))

#Identical to modpls$ttPredictY
predict(modpls,type="scores")
predict(modplsform,type="scores")

\donttest{
#Identical to modpls2$ValsPredictY
cbind(predict(modpls,newdata=Xpine_sup,type="response"),
predict(modplsform,newdata=Xpine_sup,type="response"))

#Select the number of components to use to derive the prediction
predict(modpls,newdata=Xpine_sup,type="response",comps=1)    
predict(modpls,newdata=Xpine_sup,type="response",comps=3)    
predict(modpls,newdata=Xpine_sup,type="response",comps=6)    
try(predict(modpls,newdata=Xpine_sup,type="response",comps=8))

#Identical to modpls2$ttValsPredictY
predict(modpls,newdata=Xpine_sup,type="scores")    

#Select the number of components in the scores matrix
predict(modpls,newdata=Xpine_sup,type="scores",comps=1)    
predict(modpls,newdata=Xpine_sup,type="scores",comps=3)    
predict(modpls,newdata=Xpine_sup,type="scores",comps=6)    
try(predict(modpls,newdata=Xpine_sup,type="scores",comps=8))

#Identical to modpls2NA$ValsPredictY
predict(modpls,newdata=Xpine_supNA,type="response",methodNA="missingdata")    

cbind(predict(modpls,newdata=Xpine_supNA,type="response"),
predict(modplsform,newdata=Xpine_supNA,type="response"))

predict(modpls,newdata=Xpine_supNA,type="response",comps=1)    
predict(modpls,newdata=Xpine_supNA,type="response",comps=3)    
predict(modpls,newdata=Xpine_supNA,type="response",comps=6)    
try(predict(modpls,newdata=Xpine_supNA,type="response",comps=8))

#Identical to modpls2NA$ttPredictY
predict(modpls,newdata=Xpine_supNA,type="scores",methodNA="missingdata")
predict(modplsform,newdata=Xpine_supNA,type="scores",methodNA="missingdata")

predict(modpls,newdata=Xpine_supNA,type="scores")    
predict(modplsform,newdata=Xpine_supNA,type="scores")    
predict(modpls,newdata=Xpine_supNA,type="scores",comps=1)    
predict(modpls,newdata=Xpine_supNA,type="scores",comps=3)    
predict(modpls,newdata=Xpine_supNA,type="scores",comps=6)    
try(predict(modpls,newdata=Xpine_supNA,type="scores",comps=8))
}    
}
\keyword{methods}
\keyword{predict}

