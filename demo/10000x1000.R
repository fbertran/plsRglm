setwd("~/")
Genex <- read.csv("sample_factorial_n10000_x1000_gene.csv")
Genex <-Genex[,-(1:2)]
Geney <- read.csv("sample_factorial_n10000_x1000_link.csv")
Geney <-Geney[,-2]
Geneybin <- ifelse(Geney>=0,1,0)
Geneyrate <- exp(Geney)/(1+exp(Geney))
#Geneyrate <- pnorm(Geney)
str(Genex)
str(Geney)
install.packages("plsRglm")
library(plsRglm)


###Linear Regression####

Modlm1 <-PLS_lm_wvc(Geney,Genex,6)
ModFitted<-plsR(Geney,Genex,6)
ModFitted$AIC 
ModFitted$uscores
ModFitted$pp
ModFitted$Coeffs
ModFitted$InfCrit
bbb <- PLS_lm_kfoldcv(dataY=Geney,dataX=Genex,nt=5,keepcoeffs=TRUE)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb,MClassed=TRUE)
cv.modpls<-cv.plsR(Geney,Genex,nt=6,K=6)
system.time(PLS_lm_wvc(Geney,Genex,6),gcFirst=TRUE)


###GLM###

Modlm2 <-PLS_glm_wvc(Geneybin,Genex,6,modele="pls-glm-logistic")
GlmModFitted<-plsRglm(Geneybin,Genex,6)$AIC 
GlmModFitted$uscores
GlmModFitted$pp
GlmModFitted$Coeffs
GlmModFitted$InfCrit
bbb <- PLS_glm_kfoldcv(dataY=Geneybin,dataX=data.frame(scale(as.matrix(Genex))[,]),nt=6,K=6,NK=3,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
kfolds2CVinfos_glm(bbb)
cv.modpls<-cv.plsR(Geneybin,Genex,nt=6,K=6)
system.time(PLS_glm_wvc(Geneybin,Genex,6),gcFirst=TRUE)


### Beta Regression####
#install.packages("plasRbeta")
library(plsRbeta)
Modlm3 <-PLS_beta_wvc(Geneyrate,Genex,6)
plsbetamodel <- plsRbeta(Geneyrate,Genex,6)
plsbetamodel$AIC 
plsbetamodel$uscores
plsbetamodel$pp
plsbetamodel$Coeffs
plsbetamodel$InfCrit
bbb <- PLS_beta_kfoldcv(dataY=Geneyrate,dataX=data.frame(scale(as.matrix(Genex))[,]),nt=6,K=6,NK=3,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
kfolds2CVinfos_beta(bbb)
cv.modpls<-cv.plsR(Geneybin,Genex,nt=6,K=6)
system.time(PLS_beta_wvc(Geneybin,Genex,6),gcFirst=TRUE)






																###GPU####
															
setwd("~/")															
Genex <- read.csv("sample_factorial_n10000_x1000_gene.csv")
Genex <-Genex[,-(1:2)]
Genex <-t(Genex)
Geney <- read.csv("sample_factorial_n10000_x1000_link.csv")
Geney <-Geney[,-2]
Geneybin <- ifelse(Geney>=0,1,0)
Geneyrate <- exp(Geney)/(1+exp(Geney))
str(Genex)
str(Geney)
##install.packages("plsRglm")
library(plsRglm)
source('PLS_lm_wvc.r')
source('PLS_glm_wvc.r')
assignInNamespace('PLS_lm_wvc',PLS_lm_wvc,'plsRglm')
assignInNamespace('PLS_glm_wvc',PLS_glm_wvc,'plsRglm')


###Linear Regression####

Modlm1 <-PLS_lm_wvc(Geney,Genex,6,usegpu=TRUE)
#ModFitted<-plsR(Geney,Genex,6,usegpu=TRUE)
#ModFitted$AIC 
#ModFitted$uscores
#ModFitted$pp
#ModFitted$Coeffs
#ModFitted$InfCrit
#bbb <- PLS_lm_kfoldcv(dataY=Geney,dataX=Genex,nt=5,keepcoeffs=TRUE,usegpu=TRUE)
#kfolds2CVinfos_lm(bbb)
#kfolds2CVinfos_lm(bbb,MClassed=TRUE)
#cv.modpls<-cv.plsR(Geney,Genex,nt=6,K=6,usegpu=TRUE)
system.time(PLS_lm_wvc(Geney,Genex,6,usegpu=TRUE))


###GLM###

Modlm2 <-PLS_glm_wvc(Geneybin,Genex,6,usegpu=TRUE)
#GlmModFitted<-plsRglm(Geneybin,Genex,6,usegpu=TRUE)
#GlmModFitted$AIC 
#GlmModFitted$uscores
#GlmModFitted$pp
#GlmModFitted$Coeffs
#GlmModFitted$InfCrit
#bbb <- PLS_glm_kfoldcv(dataY=Geneybin,dataX=data.frame(scale(as.matrix(Genex))[,]),nt=6,K=6,NK=3,keepfolds=FALSE,keepdataY=TRUE,modele="pls",usegpu=TRUE)
#kfolds2CVinfos_glm(bbb)
#cv.modpls<-cv.plsR(Geneybin,Genex,nt=6,K=6,usegpu=TRUE)
system.time(PLS_glm_wvc(Geneybin,Genex,6,usegpu=TRUE),gcFirst=TRUE)


### Beta Regression####
##install.packages("plsRbeta")
library(plsRbeta)
Modlm3 <-PLS_beta_wvc(Geneyrate,Genex,6,usegpu=TRUE)
#plsbetamodel<-plsRbeta(Geneyrate,Genex,6,usegpu=TRUE)
#plsbetamodel$AIC 
#plsbetamodel$uscores
#plsbetamodel$pp
#plsbetamodel$Coeffs
#plsbetamodel$InfCrit
#bbb <- PLS_beta_kfoldcv(dataY=Geneyrate,dataX=data.frame(scale(as.matrix(Genex))[,]),nt=6,K=6,NK=3,keepfolds=FALSE,keepdataY=TRUE,modele="pls",usegpu=TRUE)
#kfolds2CVinfos_beta(bbb)
#cv.modpls<-cv.plsR(Geneyrate,Genex,nt=6,K=6,usegpu=TRUE)
system.time(PLS_beta_wvc(Geneyrate,Genex,6,usegpu=TRUE))

