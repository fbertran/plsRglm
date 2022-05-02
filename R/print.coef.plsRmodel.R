#' Print method for plsR models
#' 
#' This function provides a print method for the class \code{"coef.plsRmodel"}
#' 
#' 
#' @param x an object of the class \code{"coef.plsRmodel"}
#' @param \dots not used
#' @return \code{NULL}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{print}}
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords methods print
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
#' class(modpls)
#' print(coef(modpls))
#' rm(list=c("XCornell","yCornell","modpls"))
#' 
#' @export
print.coef.plsRmodel <- function(x,...)
{
  if(!is.null(x$Coeffs)){
    cat("Coefficients of the components\n")
    print(x$CoeffC)
    cat("Coefficients of the predictors (original scale)\n")
    print(x$Coeffs)
  }
  if(!is.null(x$Std.Coeffs)){
    cat("Coefficients of the components\n")
    print(x$CoeffC)
    cat("Coefficients of the predictors (scaled scale)\n")
    print(x$Std.Coeffs)
  }
}
