#' Coefficients for permutation bootstrap computations of PLSR models
#' 
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataset dataset to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsR}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @param verbose should info messages be displayed ?
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootpls}}.
#' @keywords models
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' 
#' # Lazraq-Cleroux PLS (Y,X) bootstrap
#' # statistic=permcoefs.plsR is the default for (Y,X) permutation resampling of PLSR models.
#' set.seed(250)
#' modpls <- plsR(yCornell,XCornell,1)
#' Cornell.bootYX <- bootpls(modpls, sim="permutation", R=250, statistic=permcoefs.plsR, verbose=FALSE)
#' 
#' @export permcoefs.plsR
permcoefs.plsR <- function(dataset,ind,nt,modele,maxcoefvalues,ifbootfail,verbose){
tempcoefs <- try(PLS_lm_wvc(dataY =dataset[ind,1], dataX=dataset[,-1], nt=nt, keepstd.coeffs=TRUE, verbose=verbose)$std.coeffs,silent=TRUE)
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
