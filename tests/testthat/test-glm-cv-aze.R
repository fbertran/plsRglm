test_that("cv.plsRglm + helpers on aze run and produce coherent objects", {
  skip_on_cran()
  data(aze, package = "plsRglm")
  set.seed(123)
  cvfit <- cv.plsRglm(y ~ ., data = aze, nt = 2,
                      modele = "pls-glm-family", 
                      family = stats::binomial(),
                      K = 3, NK = 1, random = TRUE, verbose = FALSE)
  expect_s3_class(cvfit, "cv.plsRglmmodel")
  s <- summary(cvfit)
  expect_s3_class(s, "summary.cv.plsRglmmodel")

  # cvtable should return the specialized table class
  tab <- cvtable(s, verbose = FALSE)
  expect_s3_class(tab, "table.summary.cv.plsRglmmodel")
  expect_true(length(tab) >= 1)

  # PRESS and MissClassed summaries must compute
  pr <- kfolds2Press(cvfit)
  expect_true(is.list(pr) || is.numeric(pr))
  mc <- kfolds2Mclassed(cvfit)
  expect_true(is.list(mc) || is.numeric(mc))
})
