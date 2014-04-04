\name{PLS_lm}
\alias{PLS_lm}
\title{Partial least squares Regression models with leave one out cross validation}
\description{
This function implements Partial least squares Regression models with leave one out cross validation for complete or incomplete datasets.
}
\usage{
PLS_lm(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, 
modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, 
scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12),
weights,sparse=FALSE,sparseStop=TRUE,naive=FALSE)
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
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
  \item{scaleY}{scale the response : Yes/No. Ignored since not always possible for glm responses.}
  \item{pvals.expli}{should individual p-values be reported to tune model selection ?}
  \item{alpha.pvals.expli}{level of significance for predictors when pvals.expli=TRUE}
  \item{MClassed}{number of missclassified cases, should only be used for binary responses}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{sparse}{should the coefficients of non-significant predictors (<\code{alpha.pvals.expli}) be set to 0}
  \item{sparseStop}{should component extraction stop when no significant predictors (<\code{alpha.pvals.expli}) are found}
  \item{naive}{Use the naive estimates for the Degrees of Freedom in plsR? Default is \code{FALSE}.}
}
\details{
There are several ways to deal with missing values that leads to different computations of leave one out cross validation criteria.\cr

The default estimator for Degrees of Freedom is the Kramer and Sugiyama's one. Information criteria are computed accordingly to these estimations. Naive Degrees of Freedom and Information Criteria are also provided for comparison purposes. For more details, see Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, http://arxiv.org/abs/1002.4112. 
}
\value{
  \item{nr}{Description of 'comp1'}
  \item{nc}{Description of 'comp2'}
  \item{nt}{Description of 'comp2'}
  \item{ww}{Description of 'comp2'}
  \item{wwnorm}{Description of 'comp2'}
  \item{wwetoile}{Description of 'comp2'}
  \item{tt}{Description of 'comp2'}
  \item{pp}{Description of 'comp2'}
  \item{CoeffC}{Description of 'comp2'}
  \item{uscores}{Description of 'comp2'}
  \item{YChapeau}{Description of 'comp2'}
  \item{residYChapeau}{Description of 'comp2'}
  \item{RepY}{Description of 'comp2'}
  \item{na.miss.Y}{Description of 'comp2'}
  \item{YNA}{Description of 'comp2'}
  \item{residY}{Description of 'comp2'}
  \item{ExpliX}{Description of 'comp2'}
  \item{na.miss.X}{Description of 'comp2'}
  \item{XXNA}{Description of 'comp2'}
  \item{residXX}{Description of 'comp2'}
  \item{PredictY}{Description of 'comp2'}
  \item{press.ind}{Description of 'comp2'}
  \item{press.tot}{Description of 'comp2'}
  \item{familly}{Description of 'comp2'}
  \item{ttPredictY}{Description of 'comp2'}
  \item{typeVC}{Description of 'comp2'}
  \item{dataY}{Description of 'comp2'}
  \item{computed_nt}{Description of 'comp2'}
  \item{CoeffCFull}{Description of 'comp2'}
  \item{CoeffConstante}{Description of 'comp2'}
  \item{Std.Coeffs}{Description of 'comp2'}
  \item{press.ind2}{Description of 'comp2'}
  \item{RSSresidY}{Description of 'comp2'}
  \item{Coeffs}{Description of 'comp2'}
  \item{Yresidus}{Description of 'comp2'}
  \item{RSS}{Description of 'comp2'}
  \item{residusY}{Description of 'comp2'}
  \item{AIC.std}{AIC.std vs number of components (AIC computed for the standardized model}
  \item{AIC}{AIC vs nber of components}
  \item{optional}{If the response is assumed to be binary:\cr
        i.e. \code{MClassed=TRUE}.    
  \describe{
  \item{\code{MissClassed}}{Number of missclassified results}
  \item{\code{Probs}}{Probability predicted by the model. These atre not true probabilities since they may lay outside of [0,1]}
  \item{\code{Probs.trc}}{Probability predicted by the model constrined in [0,1]}
  }
  }
  \item{ttPredictFittedMissingY}{Description of 'comp2'}
  \item{optional}{If cross validation was requested:\cr
        i.e. \code{typeVC="standard"}, \code{typeVC="missingdata"} or \code{typeVC="adaptative"}.
  \describe{
  \item{\code{R2residY}}{Description of 'comp2'}
  \item{\code{R2}}{Description of 'comp2'}
  \item{\code{press.tot2}}{Description of 'comp2'}
  \item{\code{Q2}}{Description of 'comp2'}
  \item{\code{limQ2}}{Description of 'comp2'}
  \item{\code{Q2_2}}{Description of 'comp2'}
  \item{\code{Q2cum}}{Description of 'comp2'}
  \item{\code{Q2cum_2}}{Description of 'comp2'}
  \item{\code{CVinfos}}{Description of 'comp2'}
  }
  }
  \item{optional}{If no cross validation was requested :\cr
          i.e. \code{typeVC="none"}.
  \describe{
  \item{\code{InfCrit}}{Several information criteria of fit statistics}
  }
  }
  \item{Std.ValsPredictY}{Description of 'comp2'}
  \item{ValsPredictY}{Description of 'comp2'}
  \item{Std.XChapeau}{Description of 'comp2'}
  \item{XXwotNA}{Description of 'comp2'}
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
Use \code{plsR} instead.
}
\seealso{\code{\link{PLS_lm_wvc}} and \code{\link{PLS_lm_kfoldcv}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm(yCornell,XCornell,10)$InfCrit
PLS_lm(yCornell,XCornell,10,typeVC="standard")$CVinfos
PLS_lm(yCornell,XCornell,6)$AIC 
PLS_lm(yCornell,XCornell,6)$AIC.std   

modpls <- PLS_lm(yCornell,XCornell,6,pvals.expli =TRUE)
modpls2 <- PLS_lm(yCornell,XCornell,6,sparse=TRUE)
modpls3 <- PLS_lm(yCornell,XCornell,6,sparse=TRUE,sparseStop=FALSE)
rm(list=c("XCornell","yCornell"))


\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

modpls <- PLS_lm(yaze_compl,Xaze_compl,10,MClassed=TRUE)

modpls$AIC
modpls$AIC.std
modpls$MissClassed
modpls$Probs
modpls$Probs.trc
modpls$Probs-modpls$Probs.trc

modpls$InfCrit
PLS_lm(yaze_compl,Xaze_compl,10)$InfCrit

PLS_lm(yaze_compl,Xaze_compl,10,typeVC="standard")$CVinfos
PLS_lm(yaze_compl,Xaze_compl,10,typeVC="standard",MClassed=TRUE)$CVinfos
rm(list=c("Xaze_compl","yaze_compl","modpls"))


dimX <- 24
Astar <- 2
simul_data_UniYX(dimX,Astar)
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar2 <- dataAstar2[,1]
XdataAstar2 <- dataAstar2[,2:(dimX+1)]
ysimbin1 <- dicho(ydataAstar2)
Xsimbin1 <- dicho(XdataAstar2)
PLS_lm(ysimbin1,Xsimbin1,10,MClassed=TRUE)$Probs
PLS_lm(ysimbin1,Xsimbin1,10,MClassed=TRUE)$Probs.trc
PLS_lm(ysimbin1,Xsimbin1,10,MClassed=TRUE)$MissClassed
PLS_lm(ysimbin1,Xsimbin1,10,typeVC="standard",MClassed=TRUE)$CVinfos
PLS_lm(ysimbin1,XdataAstar2,10,typeVC="standard",MClassed=TRUE)$CVinfos
PLS_lm(ydataAstar2,XdataAstar2,10,typeVC="standard")$CVinfos
rm(list=c("dimX","Astar","dataAstar2","ysimbin1","Xsimbin1","ydataAstar2","XdataAstar2"))


dimX <- 6
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar4 <- dataAstar4[,1]
XdataAstar4 <- dataAstar4[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar4,XdataAstar4,10,typeVC="standard")
modpls$computed_nt
modpls$CVinfos
str(modpls)
rm(list=c("dimX","Astar","dataAstar4","modpls","ydataAstar4","XdataAstar4"))


dimX <- 24
Astar <- 2
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar2 <- dataAstar2[,1]
XdataAstar2 <- dataAstar2[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar2,XdataAstar2,10,typeVC="standard")
modpls$computed_nt
modpls$CVinfos
rm(list=c("dimX","Astar","dataAstar2","modpls","ydataAstar2","XdataAstar2"))


# Comparing the results with the plspm package and SIMCA results in Tenenhaus's book
library(plspm)

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
plsreg1(x=XCornell,y=as.vector(yCornell),nc=3)$coeffs
PLS_lm(yCornell,XCornell,3)$uscores
PLS_lm(yCornell,XCornell,3)$pp
PLS_lm(yCornell,XCornell,3)$Coeffs

plsreg1(x=XCornell,y=as.vector(yCornell),nc=4,cv=TRUE)
PLS_lm(yCornell,XCornell,4,typeVC="standard")$press.ind
PLS_lm(yCornell,XCornell,4,typeVC="standard")$press.tot
PLS_lm(yCornell,XCornell,4,typeVC="standard")$CVinfos
plsreg1(x=XCornell,y=as.vector(yCornell),nc=4,cv=TRUE)$Q2


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_lm(log(ypine),Xpine,4)$Std.Coeffs
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=4)$std.coef
PLS_lm(log(ypine),Xpine,4)$Coeffs
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=4)$coeffs
PLS_lm(log(ypine),Xpine,1)$Std.Coeffs
PLS_lm(log(ypine),Xpine,1)$Coeffs
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=10)$Q2
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$CVinfos
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=4,cv=TRUE)$Q2


data(pine_full)
Xpine_full<-pine_full[,1:10]
ypine_full<-pine_full[,11]
PLS_lm(log(ypine_full),Xpine_full,1)$Std.Coeffs
PLS_lm(log(ypine_full),Xpine_full,1)$Coeffs
plsreg1(x=Xpine_full,y=log(as.vector(ypine_full)),nc=4)$coeffs
cor(cbind(Xpine,log(ypine)))


XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
PLS_lm(log(ypine),XpineNAX21,4)$Std.Coeffs

PLS_lm(log(ypine),XpineNAX21,4)$YChapeau[1,]
PLS_lm(log(ypine),Xpine,4)$YChapeau[1,]

PLS_lm(log(ypine),XpineNAX21,4)$CoeffC

plsreg1(x=XpineNAX21,y=as.vector(log(ypine)),nc=4,cv=TRUE)
PLS_lm(log(ypine),XpineNAX21,2,dataPredictY=XpineNAX21[1,])$ValsPredictY

PLS_lm(log(ypine),Xpine,10,typeVC="none")$InfCrit
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$CVinfos
PLS_lm(log(ypine),Xpine,10,typeVC="adaptative")$CVinfos
PLS_lm(log(ypine),Xpine,10,typeVC="missingdata")$CVinfos
PLS_lm(log(ypine),XpineNAX21,10,typeVC="none")$InfCrit
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$CVinfos
PLS_lm(log(ypine),XpineNAX21,10,typeVC="adaptative")$CVinfos
PLS_lm(log(ypine),XpineNAX21,10,typeVC="missingdata")$CVinfos

PLS_lm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeau
PLS_lm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA
rm(list=c("XCornell","yCornell","Xpine","ypine","Xpine_full","ypine_full","XpineNAX21"))



dimX <- 24
Astar <- 3
dataAstar3 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar3 <- dataAstar3[,1]
XdataAstar3 <- dataAstar3[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar3,XdataAstar3,10,typeVC="standard")
modpls$computed_nt
modpls$CVinfos
rm(list=c("dimX","Astar","dataAstar3","modpls","ydataAstar3","XdataAstar3"))


dimX <- 24
Astar <- 4
dataAstar4 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar4 <- dataAstar4[,1]
XdataAstar4 <- dataAstar4[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar4,XdataAstar4,10,typeVC="standard")
modpls$computed_nt
modpls$CVinfos
rm(list=c("dimX","Astar","dataAstar4","modpls","ydataAstar4","XdataAstar4"))


dimX <- 24
Astar <- 5
dataAstar5 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar5 <- dataAstar5[,1]
XdataAstar5 <- dataAstar5[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar5,XdataAstar5,10,typeVC="standard")
modpls$computed_nt
modpls$CVinfos
rm(list=c("dimX","Astar","dataAstar5","modpls","ydataAstar5","XdataAstar5"))


dimX <- 24
Astar <- 6
dataAstar6 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar6 <- dataAstar6[,1]
XdataAstar6 <- dataAstar6[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar6,XdataAstar6,10,typeVC="standard")
modpls$computed_nt
modpls$CVinfos
rm(list=c("dimX","Astar","dataAstar6","modpls","ydataAstar6","XdataAstar6"))
}
}
\keyword{models}
\keyword{regression}
