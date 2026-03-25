test_that("registered S3 methods are available", {
  library(utils)
  ns <- asNamespace("plsRglm")
  s3_table <- get(".__S3MethodsTable__.", envir = ns)
  lookup_s3_method <- function(gen, cls) {
    fn <- getS3method(gen, cls, optional = TRUE)
    if (is.function(fn)) {
      return(fn)
    }

    key <- paste(gen, cls, sep = ".")
    if (exists(key, envir = s3_table, inherits = FALSE)) {
      get(key, envir = s3_table, inherits = FALSE)
    } else {
      NULL
    }
  }

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
    fn <- lookup_s3_method(gen, cls)
    expect_true(is.function(fn), info = paste(gen, cls, sep = "."))
  }
})
