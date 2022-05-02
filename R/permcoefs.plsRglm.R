#' Coefficients for permutation bootstrap computations of PLSGLR models
#' 
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataset dataset to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsRglm}
#' @param family glm family to use, see \link{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @param verbose should info messages be displayed ?
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootplsglm}}.
#' @keywords models
#' @examples
#' 
#' data(Cornell)
#' 
#' # (Y,X) bootstrap of a PLSGLR model
#' # statistic=coefs.plsRglm is the default for (Y,X) bootstrap of a PLSGLR models.
#' set.seed(250)
#' modplsglm <- plsRglm(Y~.,data=Cornell,1,modele="pls-glm-family",family=gaussian)
#' Cornell.bootYX <- bootplsglm(modplsglm, R=250, typeboot="plsmodel", 
#' sim="permutation", statistic=permcoefs.plsRglm, verbose=FALSE)
#' 
#' @export permcoefs.plsRglm
permcoefs.plsRglm <- function(dataset, ind, nt, modele, family = NULL, maxcoefvalues,ifbootfail,verbose){
    tempcoefs <- try(PLS_glm_wvc(dataY = dataset[ind, 1], dataX = dataset[, 
        -1], nt = nt, modele = modele, family=family, keepstd.coeffs = TRUE, verbose=verbose)$std.coeffs, silent=TRUE)
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
