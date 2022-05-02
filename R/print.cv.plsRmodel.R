#' Print method for plsR models
#' 
#' This function provides a print method for the class \code{"cv.plsRmodel"}
#' 
#' 
#' @param x an object of the class \code{"cv.plsRmodel"}
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
#' print(cv.plsR(dataY=yCornell,dataX=XCornell,nt=10,K=6, verbose=FALSE))
#' rm(list=c("XCornell","yCornell","bbb"))
#' 
#' @export
print.cv.plsRmodel <- function(x,...)
{
  cat("Number of repeated crossvalidations:\n")
  print(length(x$results_kfolds))
  cat("Number of folds for each crossvalidation:\n")
  print(length(x$results_kfolds[[1]]))
}
