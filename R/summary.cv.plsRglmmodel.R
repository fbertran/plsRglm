#' Summary method for plsRglm models
#' 
#' This function provides a summary method for the class
#' \code{"cv.plsRglmmodel"}
#' 
#' 
#' @param object an object of the class \code{"cv.plsRglmmodel"}
#' @param \dots further arguments to be passed to or from methods.
#' @return An object of class \code{"summary.cv.plsRmodel"} if \code{model} is
#' missing or \code{model="pls"}. Otherwise an object of class
#' \code{"summary.cv.plsRglmmodel"}.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
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
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' summary(cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,
#' modele="pls-glm-family",family=gaussian(), verbose=FALSE))
#' rm(list=c("XCornell","yCornell","bbb"))
#' 
#' @export
summary.cv.plsRglmmodel <- function(object, ...)
{
  res <- kfolds2CVinfos_glm(object, ...)
  if(is.null(object$call$model)){
    class(res) <- "summary.cv.plsRmodel"
  } else {
  if(object$call$model=="pls"){
    class(res) <- "summary.cv.plsRmodel"
  } else {
    class(res) <- "summary.cv.plsRglmmodel"
  }
  }
  return(res)
}
