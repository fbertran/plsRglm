#' Table method for summary of cross validated PLSR and PLSGLR models
#' 
#' The function \code{cvtable} is wrapper of \code{cvtable.plsR} and
#' \code{cvtable.plsRglm} that provides a table summary for the classes
#' \code{"summary.cv.plsRmodel"} and \code{"summary.cv.plsRglmmodel"}
#' 
#' 
#' @aliases cvtable cvtable.plsR cvtable.plsRglm
#' @param x an object of the class \code{"summary.cv.plsRmodel"}
#' @param verbose should results be displayed ?
#' @param \dots further arguments to be passed to or from methods.
#' @return \code{list}List of Information Criteria computed for each fold.
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
#' cv.modpls <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=5)
#' res.cv.modpls <- cvtable(summary(cv.modpls))
#' plot(res.cv.modpls) #defaults to type="CVQ2"
#' rm(list=c("cv.modpls","res.cv.modpls"))
#' 
#' \donttest{
#' data(Cornell)
#' cv.modpls <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=25,verbose=FALSE)
#' res.cv.modpls <- cvtable(summary(cv.modpls))
#' plot(res.cv.modpls) #defaults to type="CVQ2"
#' rm(list=c("cv.modpls","res.cv.modpls"))
#' 	
#' data(Cornell)
#' cv.modpls <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=100,verbose=FALSE)
#' res.cv.modpls <- cvtable(summary(cv.modpls))
#' plot(res.cv.modpls) #defaults to type="CVQ2"
#' rm(list=c("cv.modpls","res.cv.modpls"))
#' 
#' data(Cornell)
#' cv.modplsglm <- cv.plsRglm(Y~.,data=Cornell,nt=6,K=6,
#' modele="pls-glm-gaussian",NK=100,verbose=FALSE)
#' res.cv.modplsglm <- cvtable(summary(cv.modplsglm))
#' plot(res.cv.modplsglm) #defaults to type="CVQ2Chi2"
#' rm(list=c("res.cv.modplsglm"))
#' }
#' 
#' @export cvtable
cvtable <- function(x,verbose=TRUE,...)
{
  if(class(x) %in% c("summary.cv.plsRmodel","summary.cv.plsRglmmodel")){
  if(class(x)=="summary.cv.plsRmodel"){return(cvtable.plsR(x,verbose=verbose))} 
  if(class(x)=="summary.cv.plsRglmmodel"){return(cvtable.plsRglm(x,verbose=verbose))}} else {
  stop("cvtable must be applied to a summary.cv.plsRmodel or summary.cv.plsRglmmodel object")
  }
}
