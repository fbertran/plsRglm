#' Plot method for table of summary of cross validated plsRglm models
#' 
#' This function provides a table method for the class
#' \code{"summary.cv.plsRglmmodel"}
#' 
#' 
#' @param x an object of the class \code{"table.summary.cv.plsRglmmodel"}
#' @param type the type of cross validation criterion to plot.
#' @param \dots further arguments to be passed to or from methods.
#' @return \code{NULL}
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
#' bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,
#' modele="pls-glm-family",family=gaussian(), verbose=FALSE)
#' plot(cvtable(summary(bbb,verbose=FALSE)),type="CVQ2Chi2")
#' rm(list=c("bbb"))
#' 
#' \donttest{
#' data(Cornell)
#' plot(cvtable(summary(cv.plsRglm(Y~.,data=Cornell,nt=10,NK=100,
#' modele="pls-glm-family",family=gaussian(), verbose=FALSE),
#' verbose=FALSE)),type="CVQ2Chi2")
#' }
#' 
#' @export
plot.table.summary.cv.plsRglmmodel <- function(x,type=c("CVMC","CVQ2Chi2","CVPreChi2"), ...)
{
  if(missing(type)){
    if("CVMC" %in% names(x)){type="CVMC"} else {type="CVQ2Chi2"}
  }
resCV=x[[type]]
mp<-barplot(resCV,col="lightblue")
text(mp, pmax(resCV/2,0.5), format(resCV/(sum(resCV)),digits = 2,nsmall=2), xpd = TRUE, col = "red")
}


