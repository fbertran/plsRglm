#' Number of missclassified individuals per group for k-fold cross validated
#' partial least squares regression models.
#' 
#' This function indicates the number of missclassified individuals per group
#' for k-fold cross validated partial least squares regression models.
#' 
#' 
#' @param pls_kfolds a k-fold cross validated partial least squares regression
#' model used on binary data
#' @return \item{list}{Number of missclassified individuals per group vs number
#' of components for the first group partition} \item{list()}{\dots{}}
#' \item{list}{Number of missclassified individuals per group vs number of
#' components for the last group partition}
#' @note Use \code{\link{cv.plsR}} or \code{\link{cv.plsRglm}} to create k-fold
#' cross validated partial least squares regression models or generalized
#' linear ones.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{kfolds2coeff}}, \code{\link{kfolds2Press}},
#' \code{\link{kfolds2Pressind}} and \code{\link{kfolds2Mclassed}} to extract
#' and transforms results from k-fold cross-validation.
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords models regression
#' @examples
#' 
#' \donttest{
#' data(aze_compl)
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' kfolds2Mclassedind(cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=1,verbose=FALSE))
#' kfolds2Mclassedind(cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=2,verbose=FALSE))
#' rm(list=c("Xaze_compl","yaze_compl"))
#' }
#' 
#' @export kfolds2Mclassedind
kfolds2Mclassedind <- function(pls_kfolds) {
  if(is.null(pls_kfolds$call$modele)){pls_kfolds$call$modele<-"pls"}

  if (is.null(pls_kfolds$call$modele) & !is.null(family)) {pls_kfolds$call$modele<-"pls-glm-family"}
  if (is.null(pls_kfolds$call$family)) {
    if (pls_kfolds$call$modele=="pls") {pls_kfolds$call$family<-NULL}
    if (pls_kfolds$call$modele=="pls-glm-Gamma") {pls_kfolds$call$family<-Gamma(link = "inverse")}
    if (pls_kfolds$call$modele=="pls-glm-gaussian") {pls_kfolds$call$family<-gaussian(link = "identity")}
    if (pls_kfolds$call$modele=="pls-glm-inverse.gaussian") {pls_kfolds$call$family<-inverse.gaussian(link = "1/mu^2")}
    if (pls_kfolds$call$modele=="pls-glm-logistic") {pls_kfolds$call$family<-binomial(link = "logit")}
    if (pls_kfolds$call$modele=="pls-glm-poisson") {pls_kfolds$call$family<-poisson(link = "log")}
    if (pls_kfolds$call$modele=="pls-glm-polr") {pls_kfolds$call$family<-NULL}
  }
  
  if (!is.null(pls_kfolds$call$family)) {
        if (is.character(pls_kfolds$call$family)) {pls_kfolds$call$family <- get(pls_kfolds$call$family, mode = "function", envir = parent.frame())}
        if (is.function(pls_kfolds$call$family)) {pls_kfolds$call$family <- pls_kfolds$call$family()}
        if (is.language(pls_kfolds$call$family)) {pls_kfolds$call$family <- eval(pls_kfolds$call$family)}
        fam_var <- pls_kfolds$call$family$variance
        fam_name <- paste(pls_kfolds$call$family$family,pls_kfolds$call$family$link)
    } else {
        if (pls_kfolds$call$modele=="pls") {
            fam_var <- function(vals) {return(1)}
            fam_name <- "pls"
        }
        if (pls_kfolds$call$modele=="pls-glm-polr") {
            fam_name <- "pls-glm-polr"
            Varyy <- function(piVaryy) {
            diag(piVaryy[-length(piVaryy)])-piVaryy[-length(piVaryy)]%*%t(piVaryy[-length(piVaryy)])
            }
            Chisqcomp <- function(yichisq,pichisq) {
            t(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])%*%MASS::ginv(Varyy(pichisq))%*%(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])
            }
            Chiscompmatrix <- function(rowspi,rowsyi) {
            sum(mapply(Chisqcomp,rowsyi,rowspi))
            }
        }
    }


if (!(pls_kfolds$call$modele=="pls-glm-polr")) {
    if (length(pls_kfolds$results_kfolds)==1) {Mclassedind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassedind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          Mclassedind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
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
                Mclassedind_kfolds[[nnkk]][[ii]] <- as.numeric(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]] < 0.5, 0, 1) != unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
            } else {
                Mclassedind_kfolds[[nnkk]][[ii]] <- attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*as.numeric(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]] < 0.5, 0, 1) != unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
            }
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                Mclassedind_kfolds[[nnkk]][[ii]] <- colSums(apply(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]] < 0.5, 0, 1),2,'!=',unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])))
            } else {
                Mclassedind_kfolds[[nnkk]][[ii]] <- colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*apply(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]] < 0.5, 0, 1),2,'!=',unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])))
                   }
            }
        }
    }
}
if (pls_kfolds$call$modele=="pls-glm-polr") {
    if (length(pls_kfolds$results_kfolds)==1) {Mclassedind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassedind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          Mclassedind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }
    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
            Mclassedind_kfolds[[nnkk]][[ii]] <- unlist(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]],lapply,apply,1,which.max)[[ii]],"!=",unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])),sum))
            }
            else {
            Mclassedind_kfolds[[nnkk]][[ii]] <- unlist(lapply(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]],lapply,apply,1,which.max)[[ii]],"!=",unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])),"*",attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")),sum))     
            }
        }
    }
rm(ii)
rm(nnkk)
}
return(Mclassedind_kfolds)
}
