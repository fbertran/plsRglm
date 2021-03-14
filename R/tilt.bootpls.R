#' Non-parametric tilted bootstrap for PLS regression models
#' 
#' Provides a wrapper for the bootstrap function \code{tilt.boot} from the
#' \code{boot} R package.\cr Implements non-parametric tilted bootstrap for PLS
#' regression models by case resampling : the \code{tilt.boot} function will
#' run an initial bootstrap with equal resampling probabilities (if required)
#' and will use the output of the initial run to find resampling probabilities
#' which put the value of the statistic at required values. It then runs an
#' importance resampling bootstrap using the calculated probabilities as the
#' resampling distribution.
#' 
#' 
#' @param object An object of class \code{plsRmodel} to bootstrap
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
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[boot:boot]{tilt.boot}}
#' @keywords models
#' @examples
#' 
#' \dontrun{
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' 
#' set.seed(1385)
#' Cornell.tilt.boot <- tilt.bootpls(plsR(yCornell,XCornell,1), statistic=coefs.plsR, 
#' typeboot="fmodel_np", R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", 
#' stype="i", index=1)
#' Cornell.tilt.boot
#' str(Cornell.tilt.boot)
#' 
#' boxplots.bootpls(Cornell.tilt.boot,indices=2:7)
#' 
#' rm(Cornell.tilt.boot)
#' }
#' 
#' @export tilt.bootpls
tilt.bootpls <- function(object, typeboot="plsmodel", statistic=coefs.plsR, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, stabvalue=1e6,...){
callplsR <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = eval(callplsR$dataY),eval(callplsR$dataX))
nt <- eval(callplsR$nt)
ifbootfail <- as.matrix(as.numeric(rep(NA, ncol(dataset))))

if(typeboot=="plsmodel"){
  temp.tilt.bootplsR <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, ...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail)}
indices.temp.tilt.bootplsR <- !is.na(temp.tilt.bootplsR$t[,1])
temp.tilt.bootplsR$t=temp.tilt.bootplsR$t[indices.temp.tilt.bootplsR,]
temp.tilt.bootplsR$R=sum(indices.temp.tilt.bootplsR)
temp.tilt.bootplsR$call$R<-sum(indices.temp.tilt.bootplsR)
return(temp.tilt.bootplsR)
}
if(typeboot=="fmodel_np"){
  dataRepYtt <- cbind(y = object$RepY,object$tt)
  wwetoile <- object$wwetoile
  #return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
  temp.tilt.bootplsR <- if(!(sim=="permutation")){tilt.boot(data=dataRepYtt, statistic=coefs.plsRnp, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, wwetoile = wwetoile, ifbootfail = ifbootfail, ...)} else {
    tilt.boot(data=dataRepYtt, statistic=permcoefs.plsRnp, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, wwetoile = wwetoile, ifbootfail = ifbootfail)}
  indices.temp.tilt.bootplsR <- !is.na(temp.tilt.bootplsR$t[,1])
  temp.tilt.bootplsR$t=temp.tilt.bootplsR$t[indices.temp.tilt.bootplsR,]
  temp.tilt.bootplsR$R=sum(indices.temp.tilt.bootplsR)
  temp.tilt.bootplsR$call$R<-sum(indices.temp.tilt.bootplsR)
return(temp.tilt.bootplsR)
}
if(typeboot=="fmodel_par"){
  temp.tilt.bootplsR <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, ...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail)}
  indices.temp.tilt.bootplsR <- !is.na(temp.tilt.bootplsR$t[,1])
  temp.tilt.bootplsR$t=temp.tilt.bootplsR$t[indices.temp.tilt.bootplsR,]
  temp.tilt.bootplsR$R=sum(indices.temp.tilt.bootplsR)
  temp.tilt.bootplsR$call$R<-sum(indices.temp.tilt.bootplsR)
  return(temp.tilt.bootplsR)
}
}
