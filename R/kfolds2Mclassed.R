#' Number of missclassified individuals for k-fold cross validated partial
#' least squares regression models.
#' 
#' This function indicates the total number of missclassified individuals for
#' k-fold cross validated partial least squares regression models.
#' 
#' 
#' @param pls_kfolds a k-fold cross validated partial least squares regression
#' model used on binary data
#' @return \item{list}{Total number of missclassified individuals vs number of
#' components for the first group partition} \item{list()}{\dots{}}
#' \item{list}{Total number of missclassified individuals vs number of
#' components for the last group partition}
#' @note Use \code{\link{cv.plsR}} to create k-fold cross validated partial
#' least squares regression models.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{kfolds2coeff}}, \code{\link{kfolds2Press}},
#' \code{\link{kfolds2Pressind}} and \code{\link{kfolds2Mclassedind}} to
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
#' \donttest{
#' data(aze_compl)
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' kfolds2Mclassed(cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=1,verbose=FALSE))
#' kfolds2Mclassed(cv.plsR(dataY=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=2,verbose=FALSE))
#' rm(list=c("Xaze_compl","yaze_compl"))
#' }
#' 
#' @export kfolds2Mclassed
kfolds2Mclassed <- function(pls_kfolds) {
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
        }
    }


if (!(pls_kfolds$call$modele=="pls-glm-polr")) {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
      max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],ncol)))
      Mclassed_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassed_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],ncol)))
          Mclassed_kfolds[[jj]] <- rep(0,max_nt[jj])
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
                Mclassed_kfolds[[nnkk]] <- Mclassed_kfolds[[nnkk]]+as.numeric(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]] < 0.5, 0, 1) != unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
                } else {
                Mclassed_kfolds[[nnkk]] <- Mclassed_kfolds[[nnkk]]+attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*as.numeric(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]] < 0.5, 0, 1) != unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
            }
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                Mclassed_kfolds[[nnkk]] <- Mclassed_kfolds[[nnkk]]+colSums(apply(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE] < 0.5, 0, 1),2,'!=',unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))) 
                } else {
                Mclassed_kfolds[[nnkk]] <- Mclassed_kfolds[[nnkk]]+colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*apply(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]] < 0.5, 0, 1),2,'!=',unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))) 
                }
            }
        }
    }
rm(ii)
rm(nnkk)
}

if (pls_kfolds$call$modele=="pls-glm-polr") {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
        max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],length)))
        Mclassed_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassed_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],length)))
          Mclassed_kfolds[[jj]] <- rep(0,max_nt[jj])
        }
      rm(jj)
      }
    }

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
            Mclassedind_kfolds[[nnkk]][[ii]] <- unlist(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]],lapply,apply,1,which.max)[[ii]],"!=",unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])),sum))[1:max_nt[nnkk]]
            }
            else {
            Mclassedind_kfolds[[nnkk]][[ii]] <- unlist(lapply(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]],lapply,apply,1,which.max)[[ii]],"!=",unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])),"*",attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")),sum))[1:max_nt[nnkk]]
          }
        }
    }

    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
            Mclassed_kfolds[[nnkk]] <- colSums(matrix(unlist(Mclassedind_kfolds[[nnkk]]),ncol=max_nt[nnkk],byrow=TRUE))
    }

rm(ii)
rm(nnkk)
}
return(Mclassed_kfolds)
}
