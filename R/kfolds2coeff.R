#' Extracts coefficients from k-fold cross validated partial least squares
#' regression models
#' 
#' This fonction extracts coefficients from k-fold cross validated partial
#' least squares regression models
#' 
#' This fonctions works for plsR and plsRglm models.
#' 
#' @param pls_kfolds an object that is a k-fold cross validated partial least
#' squares regression models either lm or glm
#' @return \item{coef.all}{matrix with the values of the coefficients for each
#' leave one out step or \code{NULL} if another type of cross validation was
#' used.}
#' @note Only for \code{NK=1} and leave one out CV
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}},
#' \code{\link{kfolds2Mclassedind}}, \code{\link{kfolds2Mclassed}} and
#' \code{\link[plsRglm:summary.plsRmodel]{summary}} to extract and transform
#' results from k-fold cross validation.
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords models regression
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=nrow(XCornell),keepcoeffs=TRUE,
#' verbose=FALSE)
#' kfolds2coeff(bbb)
#' boxplot(kfolds2coeff(bbb)[,2])
#' rm(list=c("XCornell","yCornell","bbb"))
#' 
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' bbb2 <- cv.plsR(object=ypine,dataX=Xpine,nt=4,K=nrow(Xpine),keepcoeffs=TRUE,verbose=FALSE)
#' kfolds2coeff(bbb2)
#' boxplot(kfolds2coeff(bbb2)[,1])
#' rm(list=c("Xpine","ypine","bbb2"))
#' 
#' @export kfolds2coeff
kfolds2coeff <- function(pls_kfolds) {
if (is.null(pls_kfolds$coeffs_kfolds)) {cat("No coefficients found\n"); return(NULL)}
if (length(pls_kfolds$coeffs_kfolds) != 1) {cat("Only if NK=1 and Jackknife (Loo) computations\n"); return(NULL)}
coef.all <- matrix(unlist(pls_kfolds$coeffs_kfolds[[1]]),nrow=length(pls_kfolds$coeffs_kfolds[[1]]),byrow=T)
if (is.null(pls_kfolds$folds)) {attr(coef.all,"folds") <- pls_kfolds$folds}
return(coef.all)}
