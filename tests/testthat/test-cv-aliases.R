test_that("cv.plsR accepts dataY as an alias for object", {
  set.seed(1)
  X <- data.frame(matrix(rnorm(60), nrow = 15))
  y <- rnorm(15)

  fit_alias <- cv.plsR(
    dataY = y,
    dataX = X,
    nt = 2,
    K = 3,
    NK = 1,
    random = FALSE,
    keepfolds = TRUE,
    verbose = FALSE
  )
  fit_object <- cv.plsR(
    object = y,
    dataX = X,
    nt = 2,
    K = 3,
    NK = 1,
    random = FALSE,
    keepfolds = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit_alias, "cv.plsRmodel")
  expect_equal(fit_alias$results_kfolds, fit_object$results_kfolds)
  expect_equal(fit_alias$folds, fit_object$folds)
  expect_equal(fit_alias$dataY_kfolds, fit_object$dataY_kfolds)
})

test_that("cv.plsRglm accepts dataY as an alias for object", {
  set.seed(1)
  X <- data.frame(matrix(rnorm(60), nrow = 15))
  y <- rnorm(15)

  fit_alias <- cv.plsRglm(
    dataY = y,
    dataX = X,
    nt = 2,
    K = 3,
    NK = 1,
    random = FALSE,
    modele = "pls-glm-family",
    family = stats::gaussian(),
    keepfolds = TRUE,
    verbose = FALSE
  )
  fit_object <- cv.plsRglm(
    object = y,
    dataX = X,
    nt = 2,
    K = 3,
    NK = 1,
    random = FALSE,
    modele = "pls-glm-family",
    family = stats::gaussian(),
    keepfolds = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit_alias, "cv.plsRglmmodel")
  expect_equal(fit_alias$results_kfolds, fit_object$results_kfolds)
  expect_equal(fit_alias$folds, fit_object$folds)
  expect_equal(fit_alias$dataY_kfolds, fit_object$dataY_kfolds)
})
