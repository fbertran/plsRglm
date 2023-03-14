#' @rdname cv.plsR
#' @export
cv.plsRmodel.formula <- function(object,data=NULL,nt=2,limQ2set=.0975,modele="pls", K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL,verbose=TRUE,...) 
{
if (!(modele %in% c("pls"))) {stop("Use cv.plsRglm to cross-validate PLSRGLRs")}
if (missing(data)) {data <- environment(object)}
mf0 <- match.call(expand.dots = FALSE)
m0 <- match(c("object","data","nt","limQ2set","modele","K","NK","grouplist","random","scaleX","scaleY","keepcoeffs","keepfolds","keepdataY","keepMclassed","tol_Xi","weights","subset","contrasts","verbose"), names(mf0), 0L)
mf0$formula <- mf0$object
m <- match(c("formula","data","nt","limQ2set","modele","K","NK","grouplist","random","scaleX","scaleY","keepcoeffs","keepfolds","keepdataY","keepMclassed","tol_Xi","weights","subset","contrasts","verbose"), names(mf0), 0L)
mf <- mf0[c(1L, m)]
mf[[1L]] <- as.name("PLS_lm_kfoldcv_formula")
cvmodel <- eval(mf, parent.frame())

callf0 <- match.call()
callf0$formula <- mf0$object
call0 <- c(toString(callf0[[1]]),names(callf0))
call1 <- call0[!(call0=="") & !(call0=="object")]
cvmodel$call <- callf0[call1]
cvmodel$call[[1L]] <- as.name(toString(callf0[[1]]))

  class(cvmodel) <- "cv.plsRmodel"
  return(cvmodel)
}
