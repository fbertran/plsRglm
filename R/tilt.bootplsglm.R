#' Non-parametric tilted bootstrap for PLS generalized linear regression models
#' 
#' Provides a wrapper for the bootstrap function \code{tilt.boot} from the
#' \code{boot} R package.\cr Implements non-parametric tilted bootstrap for PLS
#' generalized linear regression models by case resampling : the
#' \code{tilt.boot} function will run an initial bootstrap with equal
#' resampling probabilities (if required) and will use the output of the
#' initial run to find resampling probabilities which put the value of the
#' statistic at required values. It then runs an importance resampling
#' bootstrap using the calculated probabilities as the resampling distribution.
#' 
#' 
#' @param object An object of class \code{plsRbetamodel} to bootstrap
#' @param typeboot The type of bootstrap. Either (Y,X) boostrap
#' (\code{typeboot="plsmodel"}) or (Y,T) bootstrap
#' (\code{typeboot="fmodel_np"}). Defaults to (Y,T) resampling.
#' @param statistic A function which when applied to data returns a vector
#' containing the statistic(s) of interest. \code{statistic} must take at least
#' two arguments. The first argument passed will always be the original data.
#' The second will be a vector of indices, frequencies or weights which define
#' the bootstrap sample. Further, if predictions are required, then a third
#' argument is required which would be a vector of the random indices used to
#' generate the bootstrap predictions. Any further arguments can be passed to
#' statistic through the \code{...} argument.
#' @param R The number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case \code{R}
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' @param alpha The alpha level to which tilting is required. This parameter is
#' ignored if \code{R[1]} is 0 or if \code{theta} is supplied, otherwise it is
#' used to find the values of \code{theta} as quantiles of the initial uniform
#' bootstrap. In this case \code{R[1]} should be large enough that
#' \code{min(c(alpha, 1-alpha))*R[1] > 5}, if this is not the case then a
#' warning is generated to the effect that the \code{theta} are extreme values
#' and so the tilted output may be unreliable.
#' @param sim A character string indicating the type of simulation required.
#' Possible values are \code{"ordinary"} (the default), \code{"balanced"},
#' \code{"permutation"}, or \code{"antithetic"}.
#' @param stype A character string indicating what the second argument of
#' \code{statistic} represents. Possible values of stype are \code{"i"}
#' (indices - the default), \code{"f"} (frequencies), or \code{"w"} (weights).
#' @param index The index of the statistic of interest in the output from
#' \code{statistic}. By default the first element of the output of
#' \code{statistic} is used.
#' @param stabvalue Upper bound for the absolute value of the coefficients.
#' @param \dots ny further arguments can be passed to \code{statistic}.
#' @return An object of class "boot".
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[boot:boot]{tilt.boot}}
#' @keywords models
#' @examples
#' 
#' \donttest{
#' data(aze_compl)
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' 
#' dataset <- cbind(y=yaze_compl,Xaze_compl)
#' 
#' # Lazraq-Cleroux PLS bootstrap Classic
#' 
#' aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
#' modele="pls-glm-logistic", family=NULL), statistic=coefs.plsRglm, R=c(499, 100, 100), 
#' alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#' boxplots.bootpls(aze_compl.tilt.boot,1:2)
#' 
#' aze_compl.tilt.boot2 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
#' modele="pls-glm-logistic"), statistic=coefs.plsRglm, R=c(499, 100, 100), 
#' alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#' boxplots.bootpls(aze_compl.tilt.boot2,1:2)
#' 
#' aze_compl.tilt.boot3 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
#' modele="pls-glm-family", family=binomial), statistic=coefs.plsRglm, R=c(499, 100, 100), 
#' alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#' boxplots.bootpls(aze_compl.tilt.boot3,1:2)
#' 
#' 
#' # PLS bootstrap balanced
#' 
#' aze_compl.tilt.boot4 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
#' modele="pls-glm-logistic"), statistic=coefs.plsRglm, R=c(499, 100, 100), 
#' alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)
#' boxplots.bootpls(aze_compl.tilt.boot4,1:2)
#' }
#' 
#' @export tilt.bootplsglm
tilt.bootplsglm <- function(object, typeboot="fmodel_np", statistic=coefs.plsRglm, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, stabvalue=1e6,...){
callplsRglm <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRglm$nt)
ifbootfail <- as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset)))))
#dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(!is.null(callplsRglm$family)){family <- eval(callplsRglm$family)} else {family <- NULL}

if(typeboot=="plsmodel"){
#return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family))
  temp.tilt.bootplsRglm <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail,...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail)}
  indices.temp.tilt.bootplsRglm <- !is.na(temp.tilt.bootplsRglm$t[,1])
  temp.tilt.bootplsRglm$t=temp.tilt.bootplsRglm$t[indices.temp.tilt.bootplsRglm,]
  temp.tilt.bootplsRglm$R=sum(indices.temp.tilt.bootplsRglm)
  temp.tilt.bootplsRglm$call$R<-sum(indices.temp.tilt.bootplsRglm)
  return(temp.tilt.bootplsRglm)
}
if(typeboot=="fmodel_np"){
#return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglmnp} else {permcoefs.plsRglmnp}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family))
  dataRepYtt <- cbind(y = object$RepY,object$tt)
  wwetoile <- object$wwetoile
  temp.tilt.bootplsRglm <- if(!(sim=="permutation")){tilt.boot(data=dataRepYtt, statistic=coefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail, ...)} else {
    tilt.boot(data=dataRepYtt, statistic=permcoefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail)}
  indices.temp.tilt.bootplsRglm <- !is.na(temp.tilt.bootplsRglm$t[,1])
  temp.tilt.bootplsRglm$t=temp.tilt.bootplsRglm$t[indices.temp.tilt.bootplsRglm,]
  temp.tilt.bootplsRglm$R=sum(indices.temp.tilt.bootplsRglm)
  temp.tilt.bootplsRglm$call$R<-sum(indices.temp.tilt.bootplsRglm)
  return(temp.tilt.bootplsRglm)
}
if(typeboot=="fmodel_par"){
#return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family))
  temp.tilt.bootplsRglm <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail,...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail)}
  indices.temp.tilt.bootplsRglm <- !is.na(temp.tilt.bootplsRglm$t[,1])
  temp.tilt.bootplsRglm$t=temp.tilt.bootplsRglm$t[indices.temp.tilt.bootplsRglm,]
  temp.tilt.bootplsRglm$R=sum(indices.temp.tilt.bootplsRglm)
  temp.tilt.bootplsRglm$call$R<-sum(indices.temp.tilt.bootplsRglm)
  return(temp.tilt.bootplsRglm)
}
}
