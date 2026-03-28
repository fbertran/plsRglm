#' @rdname plsRglm
#' @aliases plsRglm
#' @export PLS_glm

PLS_glm <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,method,sparse=FALSE,sparseStop=FALSE,naive=FALSE,fit_backend="stats",verbose=TRUE) {  


##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

if(verbose){cat("____************************************************____\n")}
weights_missing <- missing(weights)
method_missing <- missing(method)
naive_missing <- missing(naive)
weights_value <- NULL
method_value <- NULL
if (!weights_missing) {
weights_value <- weights
}
if (!method_missing) {
method_value <- method
}

validation <- pls_glm_validate_inputs(
dataY = dataY,
dataX = dataX,
dataPredictY = dataPredictY,
weights = weights_value,
weights_missing = weights_missing,
method = method_value,
method_missing = method_missing,
modele = modele,
family = family,
sparse = sparse,
sparseStop = sparseStop,
pvals.expli = pvals.expli,
naive = naive,
naive_missing = naive_missing,
fit_backend = fit_backend,
verbose = verbose,
eval_env = parent.frame()
)

PredYisdataX <- validation$PredYisdataX
NoWeights <- validation$NoWeights
method <- validation$method
na.miss.X <- validation$na.miss.X
na.miss.Y <- validation$na.miss.Y
na.miss.PredictY <- validation$na.miss.PredictY
naive <- validation$naive
sparse <- validation$sparse
sparseStop <- validation$sparseStop
pvals.expli <- validation$pvals.expli
dataX <- validation$dataX
modele <- validation$modele
family <- validation$family
fit_backend <- validation$fit_backend

preprocessed <- pls_glm_preprocess_inputs(
dataY = dataY,
dataX = dataX,
dataPredictY = dataPredictY,
scaleX = scaleX,
scaleY = scaleY,
weights = weights_value,
NoWeights = NoWeights,
PredYisdataX = PredYisdataX,
modele = modele
)

scaleY <- preprocessed$scaleY
RepY <- preprocessed$RepY
ExpliX <- preprocessed$ExpliX
PredictY <- preprocessed$PredictY
XXNA <- preprocessed$XXNA
YNA <- preprocessed$YNA
PredictYNA <- preprocessed$PredictYNA
XXwotNA <- preprocessed$XXwotNA
YwotNA <- preprocessed$YwotNA
PredictYwotNA <- preprocessed$PredictYwotNA

if (modele == "pls-glm-polr") {
dataY <- pls_glm_as_polr_response(dataY)
YwotNA <- pls_glm_as_polr_response(YwotNA)}

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,RSS=rep(NA,nt),RSSresidY=rep(NA,nt),R2=rep(NA,nt),R2residY=rep(NA,nt),press.ind=NULL,press.tot=NULL,Q2cum=rep(NA, nt),family=family,fit_backend=fit_backend,ttPredictY = NULL,typeVC=typeVC,dataX=dataX,dataY=dataY) 
if(NoWeights){res$weights<-rep(1L,res$nr)} else {res$weights<-weights}
res$temppred <- NULL

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (scaleY) {res$YChapeau=rep(attr(RepY,"scaled:center"),nrow(ExpliX))
res$residYChapeau=rep(0,nrow(ExpliX))}
else
{res$YChapeau=rep(mean(RepY),nrow(ExpliX))
res$residYChapeau=rep(mean(RepY),nrow(ExpliX))}
}





################################################
################################################
##                                            ##
##  Beginning of the loop for the components  ##
##                                            ##
################################################
################################################

res$computed_nt <- 0
break_nt <- FALSE
break_nt_sparse <- FALSE
break_nt_sparse1 <- FALSE
break_nt_vc <- FALSE

for (kk in 1:nt) {
XXwotNA <- as.matrix(res$residXX)
XXwotNA[!XXNA] <- 0
if (modele == "pls-glm-polr") {
YwotNA <- pls_glm_as_polr_response(res$residY)
} else {
YwotNA <- as.matrix(res$residY)
YwotNA[!YNA] <- 0
}
tempww <- rep(0,res$nc)


temptest <- sqrt(colSums(res$residXX^2, na.rm=TRUE))
if(any(temptest<tol_Xi)) {
break_nt <- TRUE
if (is.null(names(which(temptest<tol_Xi)))) {
  if(verbose){cat(paste("Warning : ",paste(names(which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))}
} else {
  if(verbose){cat(paste("Warning : ",paste((which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))}
}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
rm(temptest)
break
}

res$computed_nt <- kk

##############################################
#                                            #
#     Weight computation for each model      #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if(NoWeights){
tempww <- t(XXwotNA)%*%YwotNA/(t(XXNA)%*%YwotNA^2)
}
if(!NoWeights){
tempww <- t(XXwotNA*weights)%*%YwotNA/(t(XXNA*weights)%*%YwotNA^2)
}
if (pvals.expli) {
sparse_step <- pls_glm_apply_sparse_filter(
tempww = tempww,
tempvalpvalstep = 2 * pnorm(-abs(tempww)),
alpha.pvals.expli = alpha.pvals.expli,
sparse = sparse,
sparseStop = sparseStop
)
tempww <- sparse_step$tempww
break_nt_sparse <- break_nt_sparse | sparse_step$break_nt_sparse
res$valpvalstep <- cbind(res$valpvalstep,sparse_step$tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,sparse_step$temppvalstep)
}
}

##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
glm_step <- pls_glm_compute_glm_weights(
modele = modele,
XXwotNA = XXwotNA,
XXNA = XXNA,
YwotNA = YwotNA,
tt = res$tt,
family = family,
method = method,
kk = kk,
pvals.expli = pvals.expli,
alpha.pvals.expli = alpha.pvals.expli,
sparse = sparse,
sparseStop = sparseStop,
fit_backend = fit_backend
)
tempww <- glm_step$tempww
break_nt_sparse <- break_nt_sparse | glm_step$break_nt_sparse
if (!is.null(glm_step$tempvalpvalstep)) {
res$valpvalstep <- cbind(res$valpvalstep,glm_step$tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,glm_step$temppvalstep)
}
}

##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
polr_step <- pls_glm_compute_glm_weights(
modele = modele,
XXwotNA = XXwotNA,
XXNA = XXNA,
YwotNA = YwotNA,
tt = res$tt,
family = family,
method = method,
kk = kk,
pvals.expli = pvals.expli,
alpha.pvals.expli = alpha.pvals.expli,
sparse = sparse,
sparseStop = sparseStop,
fit_backend = fit_backend
)
tempww <- polr_step$tempww
break_nt_sparse <- break_nt_sparse | polr_step$break_nt_sparse
if (!is.null(polr_step$tempvalpvalstep)) {
res$valpvalstep <- cbind(res$valpvalstep,polr_step$tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,polr_step$temppvalstep)
}
}


##############################################
#                                            #
# Computation of the components (model free) #
#                                            #
##############################################
component_step <- pls_glm_extract_component(
res = res,
tempww = tempww,
XXwotNA = XXwotNA,
XXNA = XXNA,
PredictYNA = PredictYNA,
PredictYwotNA = PredictYwotNA,
PredYisdataX = PredYisdataX,
na.miss.X = na.miss.X,
na.miss.Y = na.miss.Y,
na.miss.PredictY = na.miss.PredictY,
tol_Xi = tol_Xi,
kk = kk,
sparse = sparse,
break_nt_sparse = break_nt_sparse,
break_nt_sparse1 = break_nt_sparse1,
alpha.pvals.expli = alpha.pvals.expli,
verbose = verbose
)
res <- component_step$res
break_nt <- component_step$break_nt
break_nt_sparse1 <- component_step$break_nt_sparse1
if (component_step$should_break) {
break
}
tempwwnorm <- component_step$tempwwnorm
temptt <- component_step$temptt
temppp <- component_step$temppp




##############################################
#                                            #
#      Computation of the coefficients       #
#      of the model with kk components       #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (kk==1) {
tempCoeffC <- solve(t(res$tt[YNA])%*%res$tt[YNA])%*%t(res$tt[YNA])%*%YwotNA[YNA]
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- 0
} else {
if (!(na.miss.X | na.miss.Y)) {
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
else
{
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- diag(res$CoeffCFull)
res$CoeffConstante <- tempCoeffConstante
res$Std.Coeffs <- rbind(tempCoeffConstante,res$wwetoile%*%res$CoeffC)
rownames(res$Std.Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if (kk==1) {
tempconstglm <- pls_glm_fit_intercept_model(
y = YwotNA,
family = family,
fit_backend = fit_backend
)
res$AIC <- AIC(tempconstglm)
res$BIC <- AIC(tempconstglm, k = log(res$nr))
res$Coeffsmodel_vals <- rbind(summary(tempconstglm)$coefficients, matrix(rep(NA, 4 * nt), ncol = 4))
res$ChisqPearson <- crossprod(residuals(tempconstglm, type = "pearson"))
res$MissClassed <- sum(unclass(res$RepY) != ifelse(tempconstglm$fitted.values < 0.5, 0, 1))
}

tempregglm <- pls_glm_fit_score_model(
y = YwotNA,
tt = res$tt,
family = family,
fit_backend = fit_backend
)
res$AIC <- cbind(res$AIC, AIC(tempregglm))
res$BIC <- cbind(res$BIC, AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(
res$Coeffsmodel_vals,
rbind(summary(tempregglm)$coefficients, matrix(rep(NA, 4 * (nt - kk)), ncol = 4))
)
res$ChisqPearson <- c(res$ChisqPearson, crossprod(residuals(tempregglm, type = "pearson")))
res$MissClassed <- cbind(
res$MissClassed,
sum(unclass(res$RepY) != ifelse(tempregglm$fitted.values < 0.5, 0, 1))
)
tempCoeffC <- as.vector(coef(tempregglm))
if (kk == 1) {
res$CoeffCFull <- matrix(c(tempCoeffC, rep(NA, nt - kk)), ncol = 1)
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- tempCoeffConstante
} else {
res$CoeffCFull <- cbind(res$CoeffCFull, c(tempCoeffC, rep(NA, nt - kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante, tempCoeffConstante)
}
tempCoeffC <- tempCoeffC[-1]

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- rbind(tempCoeffConstante,res$wwetoile%*%res$CoeffC)
rownames(res$Std.Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######           PLS-GLM-POLR           ######
##############################################


if (modele %in% c("pls-glm-polr")) {
            Varyy <- function(piVaryy) {
            diag(piVaryy[-length(piVaryy)])-piVaryy[-length(piVaryy)]%*%t(piVaryy[-length(piVaryy)])
            }
            Chisqcomp <- function(yichisq,pichisq) {
            t(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])%*%MASS::ginv(Varyy(pichisq))%*%(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])
            }
            Chiscompmatrix <- function(rowspi,rowsyi) {
            sum(mapply(Chisqcomp,rowsyi,rowspi))
            }
if (kk==1) {
tempconstpolr <- MASS::polr(YwotNA~1,na.action=na.exclude,Hess=TRUE,method=method)
res$AIC <- AIC(tempconstpolr)
res$BIC <- AIC(tempconstpolr, k = log(res$nr))
res$MissClassed <- sum(!(unclass(predict(tempconstpolr,type="class"))==unclass(res$RepY)))
res$Coeffsmodel_vals <- rbind(summary(tempconstpolr)$coefficients,matrix(rep(NA,3*nt),ncol=3))
tempmodord <- predict(tempconstpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempconstpolr,type="probs")))),as.list(as.data.frame(t(tempmat)))))
suppressWarnings(rm(tempconstpolr))
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
tempregpolr <- MASS::polr(YwotNA~.,data=tttrain,na.action=na.exclude,Hess=TRUE,method=method)
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- matrix(c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)),ncol=1)
res$CoeffConstante <- tempCoeffConstante
} else {
if (!(na.miss.X | na.miss.Y)) {
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
tempregpolr <- MASS::polr(YwotNA~.,data=tttrain,na.action=na.exclude,Hess=TRUE,method=method)
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
else
{
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
tempregpolr <- MASS::polr(YwotNA~.,data=tttrain,na.action=na.exclude,Hess=TRUE,method=method)
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- as.matrix(rbind(as.matrix(tempCoeffConstante),res$wwetoile%*%res$CoeffC))
rownames(res$Std.Coeffs) <- c(names(tempregpolr$zeta),colnames(ExpliX))
}



##############################################
#                                            #
#       Prediction of the components         #
#     as if missing values (model free)      #
#       For cross-validating the GLM         #
#                                            #
##############################################





if (!(na.miss.X | na.miss.Y)) {

##############################################
#                                            #
#             Cross validation               #
#           without missing value            #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}
##############################################


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$residYChapeau <- tempregglm$linear.predictors
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values          
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}


##############################################
######              PLS-GLM-POLR         ######
##############################################
if (modele %in% c("pls-glm-polr")) {
tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*tempCoeffConstante
res$Coeffs <- rbind(as.matrix(tempConstante),tempCoeffs)
rownames(res$Coeffs) <- rownames(res$Std.Coeffs)
}                                                                                        


                                     
##############################################
}

else {
if (na.miss.X & !na.miss.Y) {


##############################################
#                                            #
#             Cross validation               #
#           with missing value(s)            #
#                                            #
##############################################


if (kk==1) {
  if(verbose){cat("____There are some NAs in X but not in Y____\n")}
}

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}
##############################################



##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$residYChapeau <- tempregglm$linear.predictors
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values                      
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}


##############################################
######              PLS-GLM-POLR         ######
##############################################
if (modele %in% c("pls-glm-polr")) {
tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")* tempCoeffConstante
res$Coeffs <- rbind(as.matrix(tempConstante),tempCoeffs)
rownames(res$Coeffs) <- rownames(res$Std.Coeffs)
}



##############################################
}

else {
if (kk==1) {
  if(verbose){cat("____There are some NAs both in X and Y____\n")}
}
}
}


##############################################
#                                            #
#      Update and end of loop cleaning       #
#        (Especially useful for PLS)         #
#                                            #
##############################################


##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$uscores <- cbind(res$uscores,res$residY/res$CoeffC[kk])
res$residY <- res$residY - res$tt%*%tempCoeffC 
res$residusY <- cbind(res$residusY,res$residY)


if (kk==1) {
res$AIC.std <- AIC(lm(res$RepY~1,weights=res$weights))
res$AIC.std <- cbind(res$AIC.std,AICpls(kk,res$residY,weights=res$weights))
res$AIC <- AIC(lm(dataY~1))
res$AIC <- cbind(res$AIC,AICpls(kk,res$Yresidus,weights=res$weights))
if (MClassed) {
res$MissClassed <- sum(unclass(dataY)!=ifelse(predict(lm(dataY~1,weights=res$weights)) < 0.5, 0,1))
res$MissClassed <- cbind(res$MissClassed,sum(unclass(dataY)!=ifelse(res$YChapeau < 0.5, 0,1)))
tempprob <- res$Probs <- predict(lm(dataY~1,weights=res$weights))
tempprob <- ifelse(tempprob<0,0,tempprob)
res$Probs.trc <- ifelse(tempprob>1,1,tempprob)
res$Probs <- cbind(res$Probs,res$YChapeau)
tempprob <- ifelse(res$YChapeau<0,0,res$YChapeau)
tempprob <- ifelse(tempprob>1,1,tempprob)
res$Probs.trc <- cbind(res$Probs.trc,tempprob)
}
} else {
res$AIC.std <- cbind(res$AIC.std,AICpls(kk,res$residY,weights=res$weights))
res$AIC <- cbind(res$AIC,AICpls(kk,res$Yresidus,weights=res$weights))
if (MClassed) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(dataY)!=ifelse(res$YChapeau < 0.5, 0,1)))
res$Probs <- cbind(res$Probs,res$YChapeau)
tempprob <- ifelse(res$YChapeau<0,0,res$YChapeau)
tempprob <- ifelse(tempprob>1,1,tempprob)
res$Probs.trc <- cbind(res$Probs.trc,tempprob)
}
}


rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}

##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}


if(verbose){cat("____Component____",kk,"____\n")}
}




##############################################
##############################################
##                                          ##
##    End of the loop on the components     ##
##                                          ##
##############################################
##############################################

if(res$computed_nt==0){
cat("No component could be extracted please check the data for NA only lines or columns\n"); stop()
}


if (pvals.expli&!(modele=="pls")) {
res$Coeffsmodel_vals<-res$Coeffsmodel_vals[1:(dim(res$Coeffsmodel_vals)[1]-(nt-res$computed_nt)),]
}


##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.PredictY | na.miss.Y)) {
  if(verbose){cat("____Predicting X without NA neither in X nor in Y____\n")}
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- NULL}
else {
if (na.miss.PredictY & !na.miss.Y) {
  if(verbose){cat("____Predicting X with NA in X and not in Y____\n")}
for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%res$pp[PredictYNA[ii,],,drop=FALSE])%*%t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}
colnames(res$ttPredictY) <- NULL
}
else {
  if(verbose){cat("____There are some NAs both in X and Y____\n")}
}
}




##############################################
#                                            #
#          Computing RSS, PRESS,             #
#           Chi2, Q2 and Q2cum               #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$R2residY <- 1-res$RSSresidY[2:(res$computed_nt+1)]/res$RSSresidY[1]
res$R2 <- 1-res$RSS[2:(res$computed_nt+1)]/res$RSS[1]
if (MClassed==FALSE) {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "AIC.std"))
res$ic.dof<-infcrit.dof(res,naive=naive)
res$InfCrit <- cbind(res$InfCrit,res$ic.dof)
} else {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(NA,res$R2), res$MissClassed, c(NA,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "AIC.std"))
res$ic.dof<-infcrit.dof(res,naive=naive)
res$InfCrit <- cbind(res$InfCrit,res$ic.dof)
}
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$R2residY <- 1-res$RSSresidY[2:(res$computed_nt+1)]/res$RSSresidY[1]
res$R2 <- 1-res$RSS[2:(res$computed_nt+1)]/res$RSS[1]
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC, res$ChisqPearson, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}
if ((modele %in% c("pls-glm-logistic"))|(family$family=="binomial")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC, res$MissClassed, res$ChisqPearson, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC", "Missclassed", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele == "pls-glm-polr") {

res$InfCrit <- t(rbind(res$AIC, res$BIC, res$MissClassed, res$ChisqPearson))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC", "Missclassed", "Chi2_Pearson_Y"))
}




##########################################
#                                        #
#          Predicting responses          #
#                                        #
##########################################


##############################################
######               PLS                ######
##############################################
if (modele == "pls") {
res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
rownames(res$YChapeau) <- rownames(ExpliX)

res$Std.ValsPredictY <- res$ttPredictY%*%res$CoeffC
res$ValsPredictY <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$ttPredictY%*%res$CoeffC

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
if (EstimXNA) {
res$XChapeau <- sweep(sweep(res$Std.XChapeau,2,attr(res$ExpliX,"scaled:scale"),FUN="*"),2,attr(res$ExpliX,"scaled:center"),FUN="+")
rownames(res$XChapeau) <- rownames(ExpliX)
colnames(res$XChapeau) <- colnames(ExpliX)

res$XChapeauNA <- sweep(sweep(res$Std.XChapeau,2,attr(res$ExpliX,"scaled:scale"),FUN="*"),2,attr(res$ExpliX,"scaled:center"),FUN="+")*!XXNA
rownames(res$XChapeau) <- rownames(ExpliX)
colnames(res$XChapeau) <- colnames(ExpliX)
}
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$YChapeau <- as.matrix(tempregglm$fitted.values)            
rownames(res$YChapeau) <- rownames(ExpliX)

res$Std.ValsPredictY <- pls_glm_predict_score_model(tempregglm, res$ttPredictY, type = "link")
res$ValsPredictY <- pls_glm_predict_score_model(tempregglm, res$ttPredictY, type = "response")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
res$FinalModel <- pls_glm_refit_compatible_model(YwotNA, res$tt, family)
attr(res$FinalModel, "fit_backend") <- fit_backend
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
res$YChapeau <- tempregpolr$fitted.values
res$YChapeauCat <- predict(tempregpolr,type="class")
rownames(res$YChapeau) <- rownames(ExpliX)

ttpred <- data.frame(tt=res$ttPredictY)
res$ValsPredictY <- predict(tempregpolr, newdata=ttpred,type="probs")
res$ValsPredictYCat <- predict(tempregpolr, newdata=ttpred,type="class")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt,sep="")
res$FinalModel <- tempregpolr
}


colnames(res$ttPredictY) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$pp) <- colnames(ExpliX)
colnames(res$pp) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$ww) <- colnames(ExpliX)
colnames(res$ww) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$wwnorm) <- colnames(ExpliX)
colnames(res$wwnorm) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$wwetoile) <- colnames(ExpliX)
colnames(res$wwetoile) <- paste("Coord_Comp_",1:res$computed_nt,sep="")
rownames(res$tt) <- rownames(ExpliX)
colnames(res$tt) <- paste("Comp_",1:res$computed_nt,sep="")
res$XXwotNA <- XXwotNA
if(verbose){cat("****________________________________________________****\n")}
if(verbose){cat("\n")}
return(res)
}
