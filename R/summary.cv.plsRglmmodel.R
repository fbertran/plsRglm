summary.cv.plsRglmmodel <- function(object, ...)
{
  res <- kfolds2CVinfos_glm(object, ...)
  return(res)
}
