#' Print method for plsRglm models
#' 
#' This function provides a print method for the class \code{"plsRglmmodel"}
#' 
#' 
#' @param x an object of the class \code{"plsRglmmodel"}
#' @param \dots not used
#' @return \code{NULL}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{print}}
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparaison de la
#' régression PLS et de la régression
#' logistique PLS : application aux données
#' d'allélotypage. \emph{Journal de la Société Française
#' de Statistique}, 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords methods print
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
#' class(modplsglm)
#' print(modplsglm)
#' rm(list=c("XCornell","yCornell","modplsglm"))
#' 
#' @export
print.plsRglmmodel <- function(x,...)
{
  cat("Number of required components:\n")
  print(x$nt)
  cat("Number of successfully computed components:\n")
  print(x$computed_nt)
  cat("Coefficients:\n")
  print(x$Coeffs)
  cat("Information criteria and Fit statistics:\n")
  print(x$InfCrit)
  if (!is.null(x$family))
  {
    cat("Model with all the required components:\n")
    print(x$FinalModel)
  }
}
