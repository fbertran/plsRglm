cvtable <- function(x,...)
{
  if(class(x) %in% c("summary.cv.plsRmodel","summary.cv.plsRglmmodel")){
  if(class(x)=="summary.cv.plsRmodel"){return(cvtable.plsR(x))} 
  if(class(x)=="summary.cv.plsRglmmodel"){return(cvtable.plsRglm(x))}} else {
  stop("cvtable must be applied to a summary.cv.plsRmodel or summary.cv.plsRglmmodel object")
  }
}
