#' Plot method for table of summary of cross validated plsR models
#' 
#' This function provides a table method for the class
#' \code{"summary.cv.plsRmodel"}
#' 
#' 
#' @param x an object of the class \code{"table.summary.cv.plsRmodel"}
#' @param type the type of cross validation criterion to plot.
#' @param \dots further arguments to be passed to or from methods.
#' @return \code{NULL}
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
#' bbb <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=5, verbose=FALSE)
#' plot(cvtable(summary(bbb)),type="CVQ2")
#' rm(list=c("bbb"))
#' 
#' \donttest{
#' data(Cornell)
#' plot(cvtable(summary(cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=100, verbose=FALSE))),type="CVQ2")
#' }
#' 
#' @export
plot.table.summary.cv.plsRmodel <- function(x,type=c("CVMC","CVQ2","CVPress"), ...)
{
if(missing(type)){
  if("CVMC" %in% names(x)){type="CVMC"} else {type="CVQ2"}
}
resCV=x[[type]]
mp<-barplot(resCV,col="lightblue")
text(mp, pmax(resCV/2,0.5), format(resCV/(sum(resCV)),digits = 2,nsmall=2), xpd = TRUE, col = "red")
}


