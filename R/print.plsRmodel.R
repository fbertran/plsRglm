#' Print method for plsR models
#' 
#' This function provides a print method for the class \code{"plsRmodel"}
#' 
#' 
#' @param x an object of the class \code{"plsRmodel"}
#' @param \dots not used
#' @return \code{NULL}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{print}}
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparaison de la
#' régression PLS et de la régression
#' logistique PLS : application aux données
#' d'allélotypage. \emph{Journal de la Société Française
#' de Statistique}, 151(2), pages 1-18.
#' \url{https://www.numdam.org/item/JSFS_2010__151_2_1_0/}
#' @keywords methods print
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
#' class(modpls)
#' print(modpls)
#' rm(list=c("XCornell","yCornell","modpls"))
#' 
#' @export
print.plsRmodel <- function(x,...)
{
  cat("Number of required components:\n")
  print(x$nt)
  cat("Number of successfully computed components:\n")
  print(x$computed_nt)
  cat("Coefficients:\n")
  print(x$Coeffs)
  if (x$typeVC=="none")
  {
    cat("Information criteria and Fit statistics:\n")
    print(x$InfCrit)
  }
  if (x$typeVC %in% c("standard","missingdata","adaptative"))
  {
    cat("Leave one out cross validated PRESS, Information criteria and Fit statistics:\n")
    print(x$InfCrit)
  }
}
