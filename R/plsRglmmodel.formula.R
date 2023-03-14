#' @rdname plsRglm
#' @aliases plsRglm
#' @export

plsRglmmodel.formula <- function(object,data=NULL,nt=2,limQ2set=.0975,dataPredictY,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,subset,start=NULL,etastart,mustart,offset,method="glm.fit",control= list(),contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE,verbose=TRUE,...)
{
  if(missing(modele)){modele="pls"}
  if (typeVC!="none") {stop("Use plsRglm_kfoldcv for applying kfold cross validation to glms")}
if (missing(data)) {data <- environment(object)}
mf0 <- match.call(expand.dots = FALSE)
m0 <- match(c("object","data","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","subset","start","etastart","mustart","offset","method","control","contrasts","sparse","sparseStop","naive","verbose"), names(mf0), 0L)
mf0$formula <- mf0$object
m <- match(c("formula","data","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","subset","start","etastart","mustart","offset","method","control","contrasts","sparse","sparseStop","naive","verbose"), names(mf0), 0L)
mf <- mf0[c(1L, m)]
if(is.null(mf$modele)){mf$modele<-"pls"}
mf[[1L]] <- as.name("PLS_glm_formula")
estmodel <- eval(mf, parent.frame())

  callf0 <- match.call()
  callf0$dataY <- mf0$object
  call0 <- c(toString(callf0[[1]]),names(callf0))
  call1 <- call0[!(call0=="") & !(call0=="object")]
  estmodel$call <- callf0[call1]
  estmodel$call[[1L]] <- as.name(toString(callf0[[1]]))

  class(estmodel) <- "plsRglmmodel"
  estmodel
}
