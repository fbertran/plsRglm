#' Computes PRESS for k-fold cross validated partial least squares regression
#' models.
#' 
#' This function computes PRESS for k-fold cross validated partial least
#' squares regression models.
#' 
#' 
#' @param pls_kfolds a k-fold cross validated partial least squares regression
#' model
#' @return \item{list}{Press vs number of components for the first group
#' partition} \item{list()}{\dots{}} \item{list}{Press vs number of components
#' for the last group partition}
#' @note Use \code{\link{cv.plsR}} to create k-fold cross validated partial
#' least squares regression models.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}},
#' \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to
#' extract and transforms results from k-fold cross validation.
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
#' kfolds2Press(cv.plsR(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),
#' nt=6,K=12,NK=1,verbose=FALSE))
#' kfolds2Press(cv.plsR(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),
#' nt=6,K=6,NK=1,verbose=FALSE))
#' rm(list=c("XCornell","yCornell"))
#' 
#' \donttest{
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' kfolds2Press(cv.plsR(dataY=ypine,dataX=Xpine,nt=10,NK=1,verbose=FALSE))
#' kfolds2Press(cv.plsR(dataY=ypine,dataX=Xpine,nt=10,NK=2,verbose=FALSE))
#' 
#' XpineNAX21 <- Xpine
#' XpineNAX21[1,2] <- NA
#' kfolds2Press(cv.plsR(dataY=ypine,dataX=XpineNAX21,nt=10,NK=1,verbose=FALSE))
#' kfolds2Press(cv.plsR(dataY=ypine,dataX=XpineNAX21,nt=10,NK=2,verbose=FALSE))
#' rm(list=c("Xpine","XpineNAX21","ypine"))
#' }
#' 
#' @export kfolds2Press
kfolds2Press <- function(pls_kfolds) {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
      max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],ncol)))
      press_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      press_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],ncol)))
          press_kfolds[[jj]] <- rep(0,max_nt[jj])
        }
      rm(jj)
      }
    }
    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]]-pls_kfolds$dataY_kfolds[[nnkk]][[ii]])^2
                } else {
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]]-pls_kfolds$dataY_kfolds[[nnkk]][[ii]])^2
                }
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2) 
                } else {
                press_kfolds[[nnkk]] <- press_kfolds[[nnkk]]+colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2)           
                }
            }
        }
    }
rm(ii)
rm(nnkk)
return(press_kfolds)
}
