test_that("fastglm backend matches stats on supported gaussian and logistic fits", {
  skip_if_not_installed("fastglm")

  data(Cornell, package = "plsRglm")
  fit_stats_gaussian <- plsRglm(
    Y ~ ., data = Cornell, nt = 2,
    modele = "pls-glm-gaussian",
    fit_backend = "stats",
    verbose = FALSE
  )
  fit_fast_gaussian <- plsRglm(
    Y ~ ., data = Cornell, nt = 2,
    modele = "pls-glm-gaussian",
    fit_backend = "fastglm",
    verbose = FALSE
  )

  expect_identical(fit_fast_gaussian$fit_backend, "fastglm")
  expect_equal(fit_fast_gaussian$CoeffC, fit_stats_gaussian$CoeffC, tolerance = 1e-6)
  expect_equal(fit_fast_gaussian$Coeffs, fit_stats_gaussian$Coeffs, tolerance = 1e-6)
  expect_equal(
    as.numeric(fit_fast_gaussian$ValsPredictY),
    as.numeric(fit_stats_gaussian$ValsPredictY),
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(predict(fit_fast_gaussian, type = "response")),
    as.numeric(predict(fit_stats_gaussian, type = "response")),
    tolerance = 1e-6
  )

  data(aze_compl, package = "plsRglm")
  fit_stats_logistic <- plsRglm(
    y ~ ., data = aze_compl, nt = 2,
    modele = "pls-glm-logistic",
    fit_backend = "stats",
    verbose = FALSE
  )
  fit_fast_logistic <- plsRglm(
    y ~ ., data = aze_compl, nt = 2,
    modele = "pls-glm-logistic",
    fit_backend = "fastglm",
    verbose = FALSE
  )

  expect_identical(fit_fast_logistic$fit_backend, "fastglm")
  expect_equal(fit_fast_logistic$CoeffC, fit_stats_logistic$CoeffC, tolerance = 1e-5)
  expect_equal(fit_fast_logistic$Coeffs, fit_stats_logistic$Coeffs, tolerance = 1e-5)
  expect_equal(
    as.numeric(fit_fast_logistic$ValsPredictY),
    as.numeric(fit_stats_logistic$ValsPredictY),
    tolerance = 1e-5
  )
  expect_equal(
    as.numeric(predict(fit_fast_logistic, type = "response")),
    as.numeric(predict(fit_stats_logistic, type = "response")),
    tolerance = 1e-5
  )
})

test_that("fastglm backend falls back to stats on unsupported fits", {
  skip_if_not_installed("fastglm")

  data(aze, package = "plsRglm")
  expect_warning(
    fit_missing <- plsRglm(
      y ~ ., data = aze, nt = 1,
      modele = "pls-glm-logistic",
      fit_backend = "fastglm",
      verbose = FALSE
    ),
    "Falling back"
  )
  expect_identical(fit_missing$fit_backend, "stats")

  data(Cornell, package = "plsRglm")
  obs_weights <- rep(1, nrow(Cornell))
  obs_weights[1] <- 2
  expect_warning(
    fit_weighted <- plsRglm(
      Cornell$Y, Cornell[, 1:7], nt = 1,
      modele = "pls-glm-gaussian",
      weights = obs_weights,
      fit_backend = "fastglm",
      verbose = FALSE
    ),
    "Falling back"
  )
  expect_identical(fit_weighted$fit_backend, "stats")

  expect_warning(
    fit_pvals <- plsRglm(
      Y ~ ., data = Cornell, nt = 1,
      modele = "pls-glm-gaussian",
      pvals.expli = TRUE,
      fit_backend = "fastglm",
      verbose = FALSE
    ),
    "Falling back"
  )
  expect_identical(fit_pvals$fit_backend, "stats")

  data(bordeaux, package = "plsRglm")
  expect_warning(
    fit_polr <- plsRglm(
      factor(bordeaux$Quality, ordered = TRUE),
      bordeaux[, 1:4],
      nt = 1,
      modele = "pls-glm-polr",
      fit_backend = "fastglm",
      verbose = FALSE
    ),
    "Falling back"
  )
  expect_identical(fit_polr$fit_backend, "stats")
})

test_that("cv.plsRglm stores the backend that actually runs", {
  skip_if_not_installed("fastglm")

  data(Cornell, package = "plsRglm")
  set.seed(123)
  cv_fast <- cv.plsRglm(
    Y ~ ., data = Cornell, nt = 2,
    modele = "pls-glm-gaussian",
    K = 3, NK = 1, random = TRUE,
    fit_backend = "fastglm",
    verbose = FALSE
  )
  expect_identical(cv_fast$fit_backend, "fastglm")
  expect_s3_class(summary(cv_fast), "summary.cv.plsRglmmodel")

  data(aze, package = "plsRglm")
  set.seed(123)
  cv_fallback <- cv.plsRglm(
    y ~ ., data = aze, nt = 1,
    modele = "pls-glm-logistic",
    K = 3, NK = 1, random = TRUE,
    fit_backend = "fastglm",
    verbose = FALSE
  )
  expect_identical(cv_fallback$fit_backend, "stats")
})

test_that("bootplsglm reuses the stored backend for built-in statistics", {
  skip_if_not_installed("fastglm")
  skip_on_cran()

  data(Cornell, package = "plsRglm")
  fit <- plsRglm(
    Y ~ ., data = Cornell, nt = 1,
    modele = "pls-glm-gaussian",
    fit_backend = "fastglm",
    verbose = FALSE
  )

  set.seed(123)
  bt <- bootplsglm(fit, typeboot = "plsmodel", R = 2, verbose = FALSE)
  expect_s3_class(bt, "boot")
  expect_true(length(bt$t0) > 0)

  set.seed(123)
  bt_raw <- bootplsglm(
    fit,
    typeboot = "plsmodel",
    statistic = coefs.plsRglm.raw,
    R = 2,
    verbose = FALSE
  )
  expect_s3_class(bt_raw, "boot")
  expect_true(length(bt_raw$t0) > 0)
})
