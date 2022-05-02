#' coef method for plsR models
#' 
#' This function provides a coef method for the class \code{"plsRmodel"}
#' 
#' 
#' @param object an object of the class \code{"plsRmodel"}
#' @param type if \code{scaled}, the coefficients of the predictors are given
#' for the scaled predictors, if \code{original} the coefficients are to be
#' used with the predictors on their original scale.
#' @param \dots not used
#' @return An object of class \code{coef.plsRmodel}.\cr
#' \item{CoeffC}{Coefficients of the components.}
#' \item{Std.Coeffs}{Coefficients of the scaled predictors.}
#' \item{Coeffs}{Coefficients of the untransformed predictors (on their
#' original scale).}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{coef}}
#' @keywords methods coef
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
#' class(modpls)
#' coef(modpls)
#' coef(modpls,type="scaled")
#' rm(list=c("XCornell","yCornell","modpls"))
#' 
#' @export
coef.plsRmodel <- function(object,type=c("scaled","original"),...)
{
  if(missing(type)){type="original"}
  if(type=="original"){
    res<-list(CoeffC=object$CoeffC,Coeffs=object$Coeffs)
  }
  if(type=="scaled"){
    res<-list(CoeffC=object$CoeffC,Std.Coeffs=object$Std.Coeffs)
  }
  class(res) <- "coef.plsRmodel"
  res
}
