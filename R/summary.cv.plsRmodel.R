#' Summary method for plsR models
#' 
#' This function provides a summary method for the class \code{"cv.plsRmodel"}
#' 
#' 
#' @param object an object of the class \code{"cv.plsRmodel"}
#' @param \dots further arguments to be passed to or from methods.
#' @return An object of class \code{"summary.cv.plsRglmmodel"}.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{summary}}
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
#' summary(cv.plsR(Y~.,data=Cornell,nt=10,K=6, verbose=FALSE), verbose=FALSE)
#' 
#' @export
summary.cv.plsRmodel <- function(object, ...)
{
  res <- kfolds2CVinfos_lm(object, ...)
  class(res) <- "summary.cv.plsRmodel"
  return(res)
}
