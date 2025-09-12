test_that("registered S3 methods are available", {
  library(utils)
  pairs <- list(
    c("plsRmodel", "default"),
    c("plsRmodel", "formula"),
    c("plsRglmmodel", "default"),
    c("plsRglmmodel", "formula"),
    c("cv.plsRmodel", "default"),
    c("cv.plsRmodel", "formula"),
    c("cv.plsRglmmodel", "default"),
    c("cv.plsRglmmodel", "formula"),
    c("coef", "plsRmodel"),
    c("coef", "plsRglmmodel"),
    c("plot", "table.summary.cv.plsRglmmodel"),
    c("plot", "table.summary.cv.plsRmodel"),
    c("print", "coef.plsRmodel"),
    c("print", "coef.plsRglmmodel"),
    c("print", "cv.plsRmodel"),
    c("print", "cv.plsRglmmodel"),
    c("summary", "cv.plsRmodel"),
    c("summary", "cv.plsRglmmodel"),
    c("predict", "plsRmodel"),
    c("predict", "plsRglmmodel"),
    c("print", "plsRmodel"),
    c("print", "plsRglmmodel"),
    c("summary", "plsRmodel"),
    c("summary", "plsRglmmodel"),
    c("print", "summary.plsRmodel"),
    c("print", "summary.plsRglmmodel")
  )
  for (pc in pairs) {
    gen <- gsub('"', "", pc[1])
    cls <- gsub('"', "", pc[2])
    fn <- getS3method(gen, cls, optional = TRUE)
    expect_true(is.function(fn), info = paste(gen, cls, sep = "."))
  }
})
