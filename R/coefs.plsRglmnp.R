#' Coefficients for bootstrap computations of PLSGLR models
#' 
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataRepYtt components' coordinates to bootstrap
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsRglm}
#' @param family glm family to use, see \link{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param wwetoile values of the Wstar matrix in the original fit
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @note ~~some notes~~
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootplsglm}}
#' @keywords models
#' @examples
#' 
#' data(Cornell)
#' 
#' # (Y,X) bootstrap of a PLSGLR model
#' # statistic=coefs.plsRglm is the default for (Y,X) bootstrap of a PLSGLR models.
#' set.seed(250)
#' modplsglm <- plsRglm(Y~.,data=Cornell,1,modele="pls-glm-family",family=gaussian)
#' Cornell.bootYT <- bootplsglm(modplsglm, R=250, statistic=coefs.plsRglmnp, verbose=FALSE)
#' 
#' @export coefs.plsRglmnp
coefs.plsRglmnp <- function(dataRepYtt, ind, nt, modele, family = NULL, maxcoefvalues, wwetoile, ifbootfail) 
{
dataRepYb=dataRepYtt[ind,1]
Tb=dataRepYtt[ind,-1]
tempCb=try(solve(t(Tb)%*%Tb)%*%t(Tb)%*%dataRepYb,silent=TRUE)
tempcoefs <- wwetoile%*%tempCb
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
