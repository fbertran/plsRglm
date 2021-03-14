#' Extracts and computes information criteria and fits statistics for k-fold
#' cross validated partial least squares models
#' 
#' This function extracts and computes information criteria and fits statistics
#' for k-fold cross validated partial least squares models for both formula or
#' classic specifications of the model.
#' 
#' The Mclassed option should only set to \code{TRUE} if the response is
#' binary.
#' 
#' @param pls_kfolds an object computed using \code{\link{PLS_lm_kfoldcv}}
#' @param MClassed should number of miss classed be computed
#' @param verbose should infos be displayed ?
#' @return \item{list}{table of fit statistics for first group partition}
#' \item{list()}{\dots{}} \item{list}{table of fit statistics for last group
#' partition}
#' @note Use \code{\link{summary}} and \code{\link{cv.plsR}} instead.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}},
#' \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}} and
#' \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold
#' cross-validation.
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
#' summary(cv.plsR(Y~.,data=Cornell,nt=10,K=6,verbose=FALSE))
#' 
#' 
#' data(pine)
#' summary(cv.plsR(x11~.,data=pine,nt=10,NK=3,verbose=FALSE),verbose=FALSE)
#' data(pineNAX21)
#' summary(cv.plsR(x11~.,data=pineNAX21,nt=10,NK=3,
#' verbose=FALSE),verbose=FALSE)
#' 
#' 
#' data(aze_compl)
#' summary(cv.plsR(y~.,data=aze_compl,nt=10,K=8,NK=3,
#' verbose=FALSE),MClassed=TRUE,verbose=FALSE)
#' }
#' 
#' @export kfolds2CVinfos_lm
kfolds2CVinfos_lm <- function(pls_kfolds,MClassed=FALSE,verbose=TRUE) {
if(!(match("dataY",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("dataY", "dataX", "nt", "limQ2set", "modele", "family", "scaleX", "scaleY", "weights", "method", "sparse", "naive", "verbose"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
mf$verbose<-verbose
if(is.null(mf$modele)){mf$modele<-"pls"}
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
(mf[[1]] <- as.name("PLS_lm"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
press_kfolds <- kfolds2Press(pls_kfolds)
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}
if (mf$modele == "pls") {

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "MissClassed", "CV_MissClassed", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}
attr(CVinfos[[nnkk]],"computed_nt") <- computed_nt
    }
}
return(CVinfos)
}



if(!(match("formula",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("formula", "data", "nt", "limQ2set", "modele", "scaleX", "scaleY","weights","subset","contrasts", "method", "sparse", "naive", "verbose"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
mf$verbose<-verbose
if(is.null(mf$modele)){mf$modele<-"pls"}
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
(mf[[1]] <- as.name("PLS_lm_formula"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
press_kfolds <- kfolds2Press(pls_kfolds)
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}
if (mf$modele == "pls") {

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "MissClassed", "CV_MissClassed", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}
attr(CVinfos[[nnkk]],"computed_nt") <- computed_nt
    }
}
if(verbose){cat("\n")};
return(CVinfos)
}
}
