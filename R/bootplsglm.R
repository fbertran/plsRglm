#' Non-parametric Bootstrap for PLS generalized linear models
#' 
#' Provides a wrapper for the bootstrap function \code{boot} from the
#' \code{boot} R package.\cr Implements non-parametric bootstraps for PLS
#' Generalized Linear Regression models by either (Y,X) or (Y,T) resampling.
#' 
#' More details on bootstrap techniques are available in the help of the
#' \code{\link[boot:boot]{boot}} function.
#' 
#' @param object An object of class \code{plsRglmmodel} to bootstrap
#' @param typeboot The type of bootstrap. Either (Y,X) boostrap
#' (\code{typeboot="plsmodel"}) or (Y,T) bootstrap
#' (\code{typeboot="fmodel_np"}). Defaults to (Y,T) resampling.
#' @param R The number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case \code{R}
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' @param statistic A function which when applied to data returns a vector
#' containing the statistic(s) of interest. \code{statistic} must take at least
#' two arguments. The first argument passed will always be the original data.
#' The second will be a vector of indices, frequencies or weights which define
#' the bootstrap sample. Further, if predictions are required, then a third
#' argument is required which would be a vector of the random indices used to
#' generate the bootstrap predictions. Any further arguments can be passed to
#' statistic through the \code{...} argument.
#' @param sim A character string indicating the type of simulation required.
#' Possible values are \code{"ordinary"} (the default), \code{"balanced"},
#' \code{"permutation"}, or \code{"antithetic"}.
#' @param stype A character string indicating what the second argument of
#' \code{statistic} represents. Possible values of stype are \code{"i"}
#' (indices - the default), \code{"f"} (frequencies), or \code{"w"} (weights).
#' @param stabvalue A value to hard threshold bootstrap estimates computed from
#' atypical resamplings. Especially useful for Generalized Linear Models.
#' @param verbose should info messages be displayed ?
#' @param \dots Other named arguments for \code{statistic} which are passed
#' unchanged each time it is called. Any such arguments to \code{statistic}
#' should follow the arguments which \code{statistic} is required to have for
#' the simulation. Beware of partial matching to arguments of \code{boot}
#' listed above.
#' @return An object of class \code{"boot"}. See the Value part of the help of
#' the function \code{\link[boot:boot]{boot}}.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[boot:boot]{boot}}
#' @references A. Lazraq, R. Cleroux, and J.-P. Gauchi. (2003). Selecting both
#' latent and explanatory variables in the PLS1 regression model.
#' \emph{Chemometrics and Intelligent Laboratory Systems}, 66(2):117-126.\cr P.
#' Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS generalised linear
#' regression. \emph{Computational Statistics & Data Analysis}, 48(1):17-46.\cr
#' A. C. Davison and D. V. Hinkley. (1997). \emph{Bootstrap Methods and Their
#' Applications}. Cambridge University Press, Cambridge.
#' @keywords models
#' @examples
#' 
#' #Imputed aze dataset
#' data(aze_compl)
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' 
#' dataset <- cbind(y=yaze_compl,Xaze_compl)
#' modplsglm <- plsRglm(y~.,data=dataset,3,modele="pls-glm-logistic")
#' 
#' library(boot)
#' # Bastien (Y,T) PLS bootstrap
#' aze_compl.bootYT <- bootplsglm(modplsglm, R=250, verbose=FALSE)
#' boxplots.bootpls(aze_compl.bootYT)
#' confints.bootpls(aze_compl.bootYT)
#' plots.confints.bootpls(confints.bootpls(aze_compl.bootYT))
#' 
#' \donttest{
#' plot(aze_compl.bootYT,index=2)
#' jack.after.boot(aze_compl.bootYT, index=2, useJ=TRUE, nt=3)
#' plot(aze_compl.bootYT, index=2,jack=TRUE)
#' aze_compl.tilt.boot <- tilt.bootplsglm(modplsglm, statistic=coefs.plsRglm, 
#' R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#' 
#' # PLS bootstrap balanced
#' aze_compl.bootYT <- bootplsglm(modplsglm, sim="balanced", R=250, verbose=FALSE)
#' boxplots.bootpls(aze_compl.bootYT)
#' confints.bootpls(aze_compl.bootYT)
#' plots.confints.bootpls(confints.bootpls(aze_compl.bootYT))
#' 
#' 
#' plot(aze_compl.bootYT)
#' jack.after.boot(aze_compl.bootYT, index=1, useJ=TRUE, nt=3)
#' plot(aze_compl.bootYT,jack=TRUE)
#' aze_compl.tilt.boot <- tilt.bootplsglm(modplsglm, statistic=coefs.plsR,
#' R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)
#' 
#' 
#' # PLS permutation bootstrap
#' 
#' aze_compl.bootYT <- bootplsglm(modplsglm, sim="permutation", R=250, verbose=FALSE)
#' boxplots.bootpls(aze_compl.bootYT)
#' plot(aze_compl.bootYT)
#' 
#' 
#' #Original aze dataset with missing values
#' data(aze)
#' Xaze<-aze[,2:34]
#' yaze<-aze$y
#' 
#' library(boot)
#' modplsglm2 <- plsRglm(yaze,Xaze,3,modele="pls-glm-logistic")
#' aze.bootYT <- bootplsglm(modplsglm2, R=250, verbose=FALSE)
#' boxplots.bootpls(aze.bootYT)
#' confints.bootpls(aze.bootYT)
#' plots.confints.bootpls(confints.bootpls(aze.bootYT))
#' 
#' 
#' 
#' 
#' #Ordinal logistic regression
#' data(bordeaux)
#' Xbordeaux<-bordeaux[,1:4]
#' ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
#' dataset <- cbind(y=ybordeaux,Xbordeaux)
#' options(contrasts = c("contr.treatment", "contr.poly"))
#' modplsglm3 <- plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr")
#' bordeaux.bootYT<- bootplsglm(modplsglm3, sim="permutation", R=250, verbose=FALSE)
#' boxplots.bootpls(bordeaux.bootYT)
#' boxplots.bootpls(bordeaux.bootYT,ranget0=TRUE)
#' 
#' bordeaux.bootYT2<- bootplsglm(modplsglm3, sim="permutation", R=250, 
#' strata=unclass(ybordeaux), verbose=FALSE)
#' boxplots.bootpls(bordeaux.bootYT2,ranget0=TRUE)
#' 
#' 
#' if(require(chemometrics)){
#' data(hyptis)
#' hyptis
#' yhyptis <- factor(hyptis$Group,ordered=TRUE)
#' Xhyptis <- as.data.frame(hyptis[,c(1:6)])
#' dataset <- cbind(y=yhyptis,Xhyptis)
#' options(contrasts = c("contr.treatment", "contr.poly"))
#' modplsglm4 <- plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr")
#' hyptis.bootYT3<- bootplsglm(modplsglm4, sim="permutation", R=250, verbose=FALSE)
#' rownames(hyptis.bootYT3$t0)<-c("Sabi\nnene","Pin\nene",
#' "Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
#' boxplots.bootpls(hyptis.bootYT3)
#' boxplots.bootpls(hyptis.bootYT3,xaxisticks=FALSE)
#' boxplots.bootpls(hyptis.bootYT3,ranget0=TRUE)
#' boxplots.bootpls(hyptis.bootYT3,ranget0=TRUE,xaxisticks=FALSE)
#' }
#' }
#' 
#' @export bootplsglm
bootplsglm <- function(object, typeboot="fmodel_np", R=250, statistic=coefs.plsRglmnp, sim="ordinary", stype="i", stabvalue=1e6,verbose=TRUE,...){
callplsRglm <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRglm$nt)
ifbootfail <- as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset)))))
#dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(!is.null(callplsRglm$family)){family <- eval(callplsRglm$family)} else {family <- NULL}

if(typeboot=="plsmodel"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}

if(typeboot=="fmodel_np"){
dataRepYtt <- cbind(y = object$RepY,object$tt)
wwetoile <- object$wwetoile
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataRepYtt, statistic=coefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail, ...)} else {
boot(data=dataRepYtt, statistic=permcoefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}

if(typeboot=="fmodel_par"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}
}
