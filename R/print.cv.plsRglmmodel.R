#' Print method for plsRglm models
#' 
#' This function provides a print method for the class \code{"cv.plsRglmmodel"}
#' 
#' 
#' @param x an object of the class \code{"cv.plsRglmmodel"}
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
#' print(cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=10,NK=1,
#' modele="pls-glm-family",family=gaussian(), verbose=FALSE))
#' rm(list=c("XCornell","yCornell","bbb"))
#' 
#' @export
print.cv.plsRglmmodel <- function(x,...)
{
  cat("Number of repeated crossvalidations:\n")
  print(length(x$results_kfolds))
  cat("Number of folds for each crossvalidation:\n")
  print(length(x$results_kfolds[[1]]))
}
