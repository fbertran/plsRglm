test_that("glm on a small sample of fowlkes runs", {
  skip_on_cran()
  data(fowlkes, package = "plsRglm")
  df <- fowlkes[sample(1:nrow(fowlkes),500,replace=FALSE), ]
  Xfowlkes <- df[,2:13]
  yfowlkes <- df[,1]
  set.seed(7)
  fit <- plsRglm(yfowlkes, Xfowlkes, nt = 2, modele = "pls-glm-logistic", 
                 verbose = FALSE)
  expect_s3_class(fit, "plsRglmmodel")
})
