#' Coefficients for bootstrap computations of PLSR models
#' 
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataRepYtt components' coordinates to bootstrap
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param wwetoile values of the Wstar matrix in the original fit
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootpls}}
#' @keywords models
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' 
#' # Lazraq-Cleroux PLS (Y,X) bootstrap
#' # statistic=coefs.plsR is the default for (Y,X) resampling of PLSR models.
#' set.seed(250)
#' modpls <- plsR(yCornell,XCornell,1)
#' Cornell.bootYT <- bootpls(modpls, R=250, typeboot="fmodel_np",
#' statistic=coefs.plsRnp, verbose=FALSE)
#' 
#' @export coefs.plsRnp
coefs.plsRnp <- function(dataRepYtt,ind,nt,modele,maxcoefvalues,wwetoile, ifbootfail){
dataRepYb=dataRepYtt[ind,1]
Tb=dataRepYtt[ind,-1]
tempCb=try(solve(t(Tb)%*%Tb)%*%t(Tb)%*%dataRepYb,silent=TRUE)
tempcoefs <- rbind(Intercept=0,wwetoile%*%tempCb)
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
