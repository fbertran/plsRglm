test_that("plsRglm fits a simple logistic model on 'aze' dataset", {
  skip_on_cran()
  skip_on_ci()
  data(aze, package = "plsRglm")
  expect_true(is.data.frame(aze))
  expect_true(is.logical(aze$y) || all(aze$y %in% c(0,1)))
  # Fit with 1-2 components to keep it fast
  set.seed(123)
  fit <- plsRglm(y ~ ., data = aze, nt = 2, modele = "pls-glm-logistic",
                verbose = FALSE)
  expect_s3_class(fit, "plsRglmmodel")
  pr <- predict(fit, type = "response")
  expect_equal(length(pr), nrow(aze))
  expect_true(all(is.finite(pr)))
})
