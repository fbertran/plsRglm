#' Computes individual Predicted Chisquare for k-fold cross validated partial
#' least squares regression models.
#' 
#' This function computes individual Predicted Chisquare for k-fold cross
#' validated partial least squares regression models.
#' 
#' 
#' @param pls_kfolds a k-fold cross validated partial least squares regression
#' glm model
#' @return \item{list}{Individual PChisq vs number of components for the first
#' group partition} \item{list()}{\dots{}} \item{list}{Individual PChisq vs
#' number of components for the last group partition}
#' @note Use \code{\link{cv.plsRglm}} to create k-fold cross validated partial
#' least squares regression glm models.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{kfolds2coeff}}, \code{\link{kfolds2Press}},
#' \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Chisq}},
#' \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to
#' extract and transforms results from k-fold cross-validation.
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords models regression
#' @examples
#' \donttest{
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' bbb <- cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=16,verbose=FALSE)
#' bbb2 <- cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=5,verbose=FALSE)
#' kfolds2Chisqind(bbb)
#' kfolds2Chisqind(bbb2)
#' rm(list=c("XCornell","yCornell","bbb","bbb2"))
#' 
#' 
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' bbb <- cv.plsRglm(dataY=ypine,dataX=Xpine,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
#' bbb2 <- cv.plsRglm(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-gaussian",K=10,verbose=FALSE)
#' kfolds2Chisqind(bbb)
#' kfolds2Chisqind(bbb2)
#'                   
#' XpineNAX21 <- Xpine
#' XpineNAX21[1,2] <- NA
#' bbbNA <- cv.plsRglm(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls",K=10,verbose=FALSE)
#' kfolds2Pressind(bbbNA)
#' kfolds2Chisqind(bbbNA)
#' bbbNA2 <- cv.plsRglm(dataY=ypine,dataX=XpineNAX21,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
#' bbbNA3 <- cv.plsRglm(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian",
#' K=10,verbose=FALSE)
#' kfolds2Chisqind(bbbNA2)
#' kfolds2Chisqind(bbbNA3)
#' rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2","bbbNA3"))
#' 
#' 
#' data(aze_compl)
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' kfolds2Chisqind(cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",
#' family=binomial(),verbose=FALSE))
#' kfolds2Chisqind(cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic",
#' verbose=FALSE))
#' kfolds2Chisqind(cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-family",
#' family=binomial(),K=10,verbose=FALSE))
#' kfolds2Chisqind(cv.plsRglm(dataY=yaze_compl,dataX=Xaze_compl,nt=10,
#' modele="pls-glm-logistic",K=10,verbose=FALSE))
#' rm(list=c("Xaze_compl","yaze_compl"))
#' }
#' 
#' @export kfolds2Chisqind
kfolds2Chisqind <- function(pls_kfolds) {
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
            Chiscompmatrixweight <- function(rowspi,rowsyi) {
            (mapply(Chisqcomp,rowsyi,rowspi))
            }
        }
        if (pls_kfolds$call$modele=="pls-beta") {            
            fam_beta <- function(vals,phis) {return(vals*(1-vals)/(1+phis))}
            fam_name <- "pls-beta"
        }
    }


    if (length(pls_kfolds$results_kfolds)==1) {preChisqind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      preChisqind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          preChisqind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }

    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (pls_kfolds$call$modele=="pls-beta") {
                    if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
                    {
                    if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                        preChisqind_kfolds[[nnkk]][[ii]] <- (pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]])^2/(fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]]))
                    } else {
                        preChisqind_kfolds[[nnkk]][[ii]] <- attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]])^2/(fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]]))
                    }
                    }                
                    else {
                        if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                        preChisqind_kfolds[[nnkk]][[ii]] <- colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]]))) 
                        } else {
                        preChisqind_kfolds[[nnkk]][[ii]] <- colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]]))) 
                        }
                    }
            } else {
            if (pls_kfolds$call$modele=="pls-glm-polr") {
                    fff <- ~pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-1
                    m <- model.frame(fff, pls_kfolds$dataY_kfolds[[nnkk]][[ii]])
                    mat <- model.matrix(fff, model.frame(fff, pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
                    if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                    preChisqind_kfolds[[nnkk]][[ii]] <- (unlist(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]][[ii]],function(xxx) {as.list(as.data.frame(t(xxx)))}),
                    Chiscompmatrix,as.list(as.data.frame(t(mat))))))
                    } else {
                    preChisqind_kfolds[[nnkk]][[ii]] <- (unlist(lapply(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]][[ii]],function(xxx) {as.list(as.data.frame(t(xxx)))}),
                    Chiscompmatrixweight,as.list(as.data.frame(t(mat)))),"*",attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")),sum)))
                    }
                    rm(fff)
                    rm(m)
                    rm(mat)
            }
            else {
                    if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
                    {
                    if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                        preChisqind_kfolds[[nnkk]][[ii]] <- (pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]])^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]]))
                    } else {
                        preChisqind_kfolds[[nnkk]][[ii]] <- attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]])^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]]))
                    }
                    }                
                    else {
                        if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                        preChisqind_kfolds[[nnkk]][[ii]] <- colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]]))) 
                        } else {
                        preChisqind_kfolds[[nnkk]][[ii]] <- colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]]))) 
                        }
                    }
            }
        }}
    }
rm(ii)
rm(nnkk)
return(preChisqind_kfolds)
}
