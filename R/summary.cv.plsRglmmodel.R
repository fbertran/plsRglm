summary.cv.plsRglmmodel <- function(object, ...)
{
  res <- kfolds2CVinfos_glm(object, ...)
  class(res) <- "summary.cv.plsRglmmodel"
  res
}
