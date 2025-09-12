#' Summary method for plsR models
#' 
#' This function provides a summary method for the class \code{"plsRmodel"}
#' 
#' 
#' @param object an object of the class \code{"plsRmodel"}
#' @param \dots further arguments to be passed to or from methods.
#' @return \item{call }{function call of plsRmodel}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{summary}}
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{https://www.numdam.org/item/JSFS_2010__151_2_1_0/}
#' @keywords methods print
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modpls <- plsR(yCornell,XCornell,3,modele="pls")
#' class(modpls)
#' summary(modpls)
#' rm(list=c("XCornell","yCornell","modpls"))
#' 
#' @export
summary.plsRmodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRmodel"
res
}
