#' Extracts and computes information criteria and fits statistics for k-fold
#' cross validated partial least squares glm models
#' 
#' This function extracts and computes information criteria and fits statistics
#' for k-fold cross validated partial least squares glm models for both formula
#' or classic specifications of the model.
#' 
#' The Mclassed option should only set to \code{TRUE} if the response is
#' binary.
#' 
#' @param pls_kfolds an object computed using \code{\link{cv.plsRglm}}
#' @param MClassed should number of miss classed be computed ?
#' @param verbose should infos be displayed ?
#' @return \item{list}{table of fit statistics for first group partition}
#' \item{list()}{\dots{}} \item{list}{table of fit statistics for last group
#' partition}
#' @note Use \code{\link{summary}} and \code{\link{cv.plsRglm}} instead.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
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
#' 
#' \donttest{
#' data(Cornell)
#' summary(cv.plsRglm(Y~.,data=Cornell,
#' nt=6,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls",verbose=FALSE),MClassed=TRUE)
#' 
#' 
#' data(aze_compl)
#' summary(cv.plsR(y~.,data=aze_compl,nt=10,K=8,modele="pls",verbose=FALSE),
#' MClassed=TRUE,verbose=FALSE)
#' summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,modele="pls",verbose=FALSE),
#' MClassed=TRUE,verbose=FALSE)
#' summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,
#' modele="pls-glm-family",
#' family=gaussian(),verbose=FALSE),
#' MClassed=TRUE,verbose=FALSE)
#' summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,
#' modele="pls-glm-logistic",
#' verbose=FALSE),MClassed=TRUE,verbose=FALSE)
#' summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,
#' modele="pls-glm-family",
#' family=binomial(),verbose=FALSE),
#' MClassed=TRUE,verbose=FALSE)
#' 
#' 
#' if(require(chemometrics)){
#' data(hyptis)
#' hyptis
#' yhyptis <- factor(hyptis$Group,ordered=TRUE)
#' Xhyptis <- as.data.frame(hyptis[,c(1:6)])
#' options(contrasts = c("contr.treatment", "contr.poly"))
#' modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
#' modpls2$Coeffsmodel_vals
#' modpls2$InfCrit
#' modpls2$Coeffs
#' modpls2$std.coeffs
#' 
#' table(yhyptis,predict(modpls2$FinalModel,type="class"))
#' 
#' modpls3 <- PLS_glm(yhyptis[-c(1,2,3)],Xhyptis[-c(1,2,3),],3,modele="pls-glm-polr",
#' dataPredictY=Xhyptis[c(1,2,3),],verbose=FALSE)
#' 
#' summary(cv.plsRglm(factor(Group,ordered=TRUE)~.,data=hyptis[,-c(7,8)],nt=4,K=10,
#' random=TRUE,modele="pls-glm-polr",keepcoeffs=TRUE,verbose=FALSE),
#' MClassed=TRUE,verbose=FALSE)
#' }
#' }
#' 
#' @export kfolds2CVinfos_glm
kfolds2CVinfos_glm <- function(pls_kfolds,MClassed=FALSE,verbose=TRUE) {
if(!(match("dataY",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("dataY", "dataX", "nt", "limQ2set", "modele", "family", "scaleX", "scaleY", "weights", "method", "sparse", "naive", "verbose"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
mf$verbose<-verbose
if(is.null(mf$modele)){mf$modele<-"pls"}
if (is.null(mf$family)) {
  if (mf$modele=="pls") {mf$family<-NULL}
  if (mf$modele=="pls-glm-Gamma") {mf$family<-Gamma(link = "inverse")}
  if (mf$modele=="pls-glm-gaussian") {mf$family<-gaussian(link = "identity")}
  if (mf$modele=="pls-glm-inverse.gaussian") {mf$family<-inverse.gaussian(link = "1/mu^2")}
  if (mf$modele=="pls-glm-logistic") {mf$family<-binomial(link = "logit")}
  if (mf$modele=="pls-glm-poisson") {mf$family<-poisson(link = "log")}
  if (mf$modele=="pls-glm-polr") {mf$family<-NULL}
}
if (!is.null(mf$family)) {
  if (is.character(mf$family)) {mf$family <- get(mf$family, mode = "function", envir = parent.frame())}
  if (is.function(mf$family)) {mf$family <- mf$family()}
  if (is.language(mf$family)) {mf$family <- eval(mf$family)}
}
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
if (!is.null(mf$family)) {mf$modele <- "pls-glm-family"}
(mf[[1]] <- as.name("PLS_glm"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}

if (mf$modele == "pls") {
press_kfolds <- kfolds2Press(pls_kfolds)
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}
    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
            Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2_2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
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


    }
class(CVinfos) <- "summary.cv.plsRmodel"
}







if (as.character(pls_kfolds$call["modele"]) %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
press_kfolds <- kfolds2Press(pls_kfolds)
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
class(CVinfos) <- "summary.cv.plsRglmmodel"
}

if (as.character(pls_kfolds$call["modele"]) == "pls-glm-polr") {
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}



if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
class(CVinfos) <- "summary.cv.plsRglmmodel"
}
return(CVinfos)
}

if(!(match("formula",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("formula", "data", "nt", "limQ2set", "modele", "family", "scaleX", "scaleY", "weights","subset","start","etastart","mustart","offset","control","method","contrasts", "sparse", "naive", "verbose"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
mf$verbose<-verbose
if(is.null(mf$modele)){mf$modele<-"pls"}
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
if (mf$modele %in% c("pls","pls-glm-logistic","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson","pls-glm-polr")){mf$family <- NULL}
(mf[[1]] <- as.name("PLS_glm_formula"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}

if (mf$modele == "pls") {
press_kfolds <- kfolds2Press(pls_kfolds)
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2_2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
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


    }
class(CVinfos) <- "summary.cv.plsRmodel"
}

if (as.character(pls_kfolds$call["modele"]) %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
press_kfolds <- kfolds2Press(pls_kfolds)
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
class(CVinfos) <- "summary.cv.plsRglmmodel"
}

if (as.character(pls_kfolds$call["modele"]) == "pls-glm-polr") {
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}



if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
class(CVinfos) <- "summary.cv.plsRglmmodel"
}
if(verbose){cat("\n")};
return(CVinfos)
}
}
