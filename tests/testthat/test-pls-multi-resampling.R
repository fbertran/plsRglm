make_pls_multi_resampling_fixture <- function() {
  set.seed(202)
  X <- matrix(rnorm(80 * 4), nrow = 80, ncol = 4)
  colnames(X) <- paste0("x", 1:4)

  B <- matrix(
    c(
      1.0, -0.2,
      -0.4, 0.8,
      0.6, 0.3,
      -0.1, 0.9
    ),
    nrow = 4,
    ncol = 2,
    byrow = TRUE
  )
  Y <- X %*% B + matrix(rnorm(80 * 2, sd = 0.05), nrow = 80, ncol = 2)
  colnames(Y) <- c("y1", "y2")

  df <- data.frame(Y, X, check.names = FALSE)
  list(X = X, Y = Y, df = df)
}

test_that("cv.plsRmulti matrix interface returns coherent summary objects", {
  fixture <- make_pls_multi_resampling_fixture()

  cvfit <- cv.plsRmulti(
    fixture$Y,
    fixture$X,
    nt = 3,
    K = 4,
    NK = 2,
    keepcoeffs = TRUE,
    verbose = FALSE
  )

  expect_s3_class(cvfit, "cv.plsRmultiModel")
  expect_length(cvfit$results_kfolds, 2L)
  expect_length(cvfit$results_kfolds[[1]], 4L)
  expect_equal(dim(cvfit$results_kfolds[[1]][[1]])[2], ncol(fixture$Y))
  expect_equal(dim(cvfit$results_kfolds[[1]][[1]])[3], 3L)
  expect_equal(length(cvfit$coeffs_kfolds[[1]][[1]]), (ncol(fixture$X) + 1) * ncol(fixture$Y))

  s <- summary(cvfit, verbose = FALSE)
  expect_s3_class(s, "summary.cv.plsRmultiModel")
  expect_true(inherits(s, "summary.cv.plsRmodel"))
  expect_true(all(c("Q2_Y", "PRESS_Y", "RSS_Y") %in% colnames(s[[1]])))

  tab <- cvtable(s, verbose = FALSE)
  expect_s3_class(tab, "table.summary.cv.plsRmultiModel")
  expect_true(inherits(tab, "table.summary.cv.plsRmodel"))
  expect_true(all(c("CVQ2", "CVPress") %in% names(tab)))
})

test_that("cv.plsRmulti formula interface matches the matrix interface summary surface", {
  fixture <- make_pls_multi_resampling_fixture()

  cv_matrix <- cv.plsRmulti(
    fixture$Y,
    fixture$X,
    nt = 2,
    K = 5,
    NK = 1,
    verbose = FALSE
  )
  cv_formula <- cv.plsRmulti(
    cbind(y1, y2) ~ .,
    data = fixture$df,
    nt = 2,
    K = 5,
    NK = 1,
    verbose = FALSE
  )

  s_matrix <- summary(cv_matrix, verbose = FALSE)
  s_formula <- summary(cv_formula, verbose = FALSE)

  expect_equal(dim(s_formula[[1]]), dim(s_matrix[[1]]))
  expect_equal(rownames(s_formula[[1]]), rownames(s_matrix[[1]]))
  expect_equal(colnames(s_formula[[1]])[1:8], colnames(s_matrix[[1]])[1:8])
})

test_that("bootpls supports plsRmultiModel with YX and YT resampling", {
  fixture <- make_pls_multi_resampling_fixture()
  fit <- plsRmulti(fixture$Y, fixture$X, nt = 2, verbose = FALSE)

  bt_yx <- bootpls(fit, R = 4, verbose = FALSE)
  expect_s3_class(bt_yx, "boot")
  expect_equal(ncol(bt_yx$t), (ncol(fixture$X) + 1) * ncol(fixture$Y))
  expect_true(all(grepl("^y[12]:", rownames(bt_yx$t0))))

  bt_yt <- bootpls(fit, typeboot = "fmodel_np", R = 4, verbose = FALSE)
  expect_s3_class(bt_yt, "boot")
  expect_equal(dim(bt_yt$t), dim(bt_yx$t))

  ci <- suppressWarnings(confints.bootpls(bt_yx, indices = 1:2, typeBCa = FALSE))
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2L)
})

test_that("cv.plsRmulti rejects unsupported options", {
  fixture <- make_pls_multi_resampling_fixture()

  expect_error(
    cv.plsRmulti(fixture$Y, fixture$X, nt = 2, weights = rep(1, nrow(fixture$X)), verbose = FALSE),
    "does not support weights"
  )
  expect_error(
    cv.plsRmulti(fixture$Y, fixture$X, nt = 2, family = gaussian(), verbose = FALSE),
    "does not support GLM families"
  )
  expect_error(
    cv.plsRmulti(fixture$Y, fixture$X, nt = 2, keepMclassed = TRUE, verbose = FALSE),
    "does not support classification diagnostics"
  )
  expect_error(
    cv.plsRmulti(fixture$Y, replace(fixture$X, 1, NA), nt = 2, verbose = FALSE),
    "missing values in X"
  )
})
