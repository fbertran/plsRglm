#' Bootstrap confidence intervals
#' 
#' This function is a wrapper for \code{\link[boot:boot.ci]{boot.ci}} to derive
#' bootstrap-based confidence intervals from a \code{"boot"} object.
#' 
#' 
#' @param bootobject an object of class \code{"boot"}
#' @param indices the indices of the predictor for which CIs should be
#' calculated. Defaults to \code{NULL}: all the predictors will be used.
#' @param typeBCa shall BCa bootstrap based CI derived ? Defaults to
#' \code{TRUE}. This is a safety option since sometimes computing BCa bootstrap
#' based CI fails whereas the other types of CI can still be derived.
#' @return Matrix with the limits of bootstrap based CI for all (defaults) or
#' only the selected predictors (\code{indices} option). The limits are given
#' in that order: Normal Lower then Upper Limit, Basic Lower then Upper Limit,
#' Percentile Lower then Upper Limit, BCa Lower then Upper Limit.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootpls}} and \code{\link{bootplsglm}}.
#' @keywords regression models
#' @examples
#' 
#' \donttest{
#' data(Cornell)
#' 
#' #Lazraq-Cleroux PLS (Y,X) bootstrap
#' set.seed(250)
#' modpls <- plsR(Y~.,data=Cornell,3)
#' Cornell.bootYX <- bootpls(modpls, R=250, verbose=FALSE)
#' confints.bootpls(Cornell.bootYX,2:8)
#' confints.bootpls(Cornell.bootYX,2:8,typeBCa=FALSE)
#' }
#' 
#' @export confints.bootpls
confints.bootpls <- function(bootobject,indices=NULL,typeBCa=TRUE){
nr <- length(bootobject$t0)
if(is.null(indices)){indices <- 1:nr}
ii <- indices[1]
if(typeBCa){
temptemp.ci <- c(boot::boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot::boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot::boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)],boot::boot.ci(bootobject, conf = 0.95, type = c("bca"), index=ii)$bca[,-c(1,2,3)])
for(ii in indices[-1]){
temptemp.ci <- rbind(temptemp.ci,c(boot::boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot::boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot::boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)],boot::boot.ci(bootobject, conf = 0.95, type = c("bca"), index=ii)$bca[,-c(1,2,3)]))
}
} else {
temptemp.ci <- c(boot::boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot::boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot::boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)])
for(ii in indices[-1]){
temptemp.ci <- rbind(temptemp.ci,c(boot::boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot::boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot::boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)]))
}
}
attr(temptemp.ci, "typeBCa") <- typeBCa
rownames(temptemp.ci) <- rownames(bootobject$t0)[indices]
return(temptemp.ci)
}
