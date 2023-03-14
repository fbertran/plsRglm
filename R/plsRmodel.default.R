#' @rdname plsR
#' @export
plsRmodel.default <- function(object,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,sparse=FALSE,sparseStop=TRUE,naive=FALSE,verbose=TRUE,...)
{
  if(missing(modele)){modele="pls"}
  if (!(modele %in% c("pls"))) {stop("Use plsRglm for applying PLSR to glms")}
  mf0 <- match.call(expand.dots = FALSE)
  m0 <- match(c("object","dataX","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","sparse","sparseStop","naive","verbose"), names(mf0), 0L)
  mf0$dataY <- mf0$object
  m <- match(c("dataY","dataX","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","sparse","sparseStop","naive","verbose"), names(mf0), 0L)
  mf <- mf0[c(1L, m)]
  if(is.null(mf$modele)){mf$modele<-"pls"}
  mf[[1L]] <- as.name("PLS_lm")
  
  estmodel <- eval(mf, parent.frame())
  
  callf0 <- match.call()
  callf0$dataY <- mf0$object
  call0 <- c(toString(callf0[[1]]),names(callf0))
  call1 <- call0[!(call0=="") & !(call0=="object")]
  estmodel$call <- callf0[call1]
  estmodel$call[[1L]] <- as.name(toString(callf0[[1]]))
  
  class(estmodel) <- "plsRmodel"
  return(estmodel)
}
