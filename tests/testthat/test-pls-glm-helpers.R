validate_pls_glm <- function(...) {
getFromNamespace("pls_glm_validate_inputs", "plsRglm")(...)
}

preprocess_pls_glm <- function(...) {
getFromNamespace("pls_glm_preprocess_inputs", "plsRglm")(...)
}

glm_weights_pls_glm <- function(...) {
getFromNamespace("pls_glm_compute_glm_weights", "plsRglm")(...)
}

extract_component_pls_glm <- function(...) {
getFromNamespace("pls_glm_extract_component", "plsRglm")(...)
}

test_that("pls_glm_validate_inputs flags missing-data constraints", {
  dataX <- data.frame(x1 = c(1, NA, 3), x2 = c(2, 4, 6))
  out <- validate_pls_glm(
    dataY = c(1, 2, NA),
    dataX = dataX,
    dataPredictY = dataX,
    weights = c(1, 2, 1),
    weights_missing = FALSE,
    method = NULL,
    method_missing = TRUE,
    modele = "pls",
    family = NULL,
    sparse = TRUE,
    sparseStop = FALSE,
    pvals.expli = FALSE,
    naive = FALSE,
    naive_missing = TRUE,
    verbose = FALSE
  )

  expect_true(out$PredYisdataX)
  expect_false(out$NoWeights)
  expect_true(out$na.miss.X)
  expect_true(out$na.miss.Y)
  expect_true(out$naive)
  expect_false(out$sparse)
  expect_false(out$sparseStop)
  expect_false(out$pvals.expli)
  expect_identical(out$method, "logistic")
  expect_s3_class(out$dataX, "data.frame")
})

test_that("pls_glm_validate_inputs rejects all-missing rows", {
  expect_error(
    validate_pls_glm(
      dataY = c(1, 2),
      dataX = data.frame(x1 = c(1, NA), x2 = c(2, NA)),
      dataPredictY = data.frame(x1 = c(1, 0), x2 = c(2, 1)),
      weights = NULL,
      weights_missing = TRUE,
      method = NULL,
      method_missing = TRUE,
      modele = "pls",
      family = NULL,
      sparse = FALSE,
      sparseStop = FALSE,
      pvals.expli = FALSE,
      naive = FALSE,
      naive_missing = TRUE,
      verbose = FALSE
    ),
    "One of the rows of dataX is completely filled with missing data"
  )
})

test_that("PLS_glm surfaces informative errors for fully missing dataX columns and rows", {
  expect_error(
    PLS_glm(
      dataY = c(1, 2),
      dataX = data.frame(x1 = c(NA, NA), x2 = c(1, 2)),
      nt = 1,
      verbose = FALSE
    ),
    "One of the columns of dataX is completely filled with missing data"
  )

  expect_error(
    PLS_glm(
      dataY = c(1, 2),
      dataX = data.frame(x1 = c(1, NA), x2 = c(2, NA)),
      nt = 1,
      verbose = FALSE
    ),
    "One of the rows of dataX is completely filled with missing data"
  )
})

test_that("pls_glm_preprocess_inputs applies weighted scaling and NA masks", {
  weights <- c(1, 2, 1)
  dataX <- data.frame(a = c(1, 2, 4), b = c(2, 6, 8))
  dataPredictY <- data.frame(a = c(2, NA), b = c(4, 10))
  out <- preprocess_pls_glm(
    dataY = c(1, 3, 5),
    dataX = dataX,
    dataPredictY = dataPredictY,
    scaleX = TRUE,
    scaleY = NULL,
    weights = weights,
    NoWeights = FALSE,
    PredYisdataX = FALSE,
    modele = "pls"
  )

  expected_center_x <- apply(dataX, 2, weighted.mean, weights)
  expected_scale_x <- sqrt(
    ((length(weights) - 1) / length(weights)) *
      apply((sweep(dataX, 2, expected_center_x))^2, 2, weighted.mean, weights)
  )
  expected_center_y <- weighted.mean(c(1, 3, 5), weights)
  expected_scale_y <- sqrt(
    ((length(weights) - 1) / length(weights)) *
      weighted.mean((c(1, 3, 5) - expected_center_y)^2, weights)
  )

  expect_true(out$scaleY)
  expect_equal(attr(out$ExpliX, "scaled:center"), expected_center_x)
  expect_equal(attr(out$ExpliX, "scaled:scale"), expected_scale_x)
  expect_equal(attr(out$RepY, "scaled:center"), expected_center_y)
  expect_equal(attr(out$RepY, "scaled:scale"), expected_scale_y)
  expect_false(out$PredictYNA[2, 1])
  expect_equal(unname(out$PredictYwotNA[2, 1]), 0)
})

test_that("pls_glm_compute_glm_weights handles sparse filtering and sparse stop", {
  XXwotNA <- cbind(
    x1 = 1:8,
    x2 = c(2, 5, 6, 9, 10, 13, 14, 17),
    x3 = c(1, -1, 1, -1, 1, -1, 1, -1)
  )
  XXNA <- matrix(TRUE, nrow = nrow(XXwotNA), ncol = ncol(XXwotNA))
  YwotNA <- c(1, 2, 3, 4, 5, 6, 7, 8)

  sparse_fit <- glm_weights_pls_glm(
    modele = "pls-glm-gaussian",
    XXwotNA = XXwotNA,
    XXNA = XXNA,
    YwotNA = YwotNA,
    tt = NULL,
    family = stats::gaussian(link = "identity"),
    method = "glm.fit",
    kk = 1,
    pvals.expli = TRUE,
    alpha.pvals.expli = 0.05,
    sparse = TRUE,
    sparseStop = FALSE
  )

  expect_true(sum(sparse_fit$temppvalstep) >= 2)
  expect_equal(sparse_fit$tempww[3], 0)
  expect_false(sparse_fit$break_nt_sparse)

  sparse_stop_fit <- glm_weights_pls_glm(
    modele = "pls-glm-gaussian",
    XXwotNA = XXwotNA,
    XXNA = XXNA,
    YwotNA = YwotNA,
    tt = NULL,
    family = stats::gaussian(link = "identity"),
    method = "glm.fit",
    kk = 1,
    pvals.expli = TRUE,
    alpha.pvals.expli = 0,
    sparse = TRUE,
    sparseStop = TRUE
  )

  expect_true(sparse_stop_fit$break_nt_sparse)
})

test_that("pls_glm_extract_component appends model-free component pieces", {
  XXwotNA <- matrix(c(
    1, 2,
    2, 1,
    3, 1,
    4, 2
  ), ncol = 2, byrow = TRUE)
  XXNA <- matrix(TRUE, nrow = 4, ncol = 2)
  res <- list(
    nr = 4,
    nc = 2,
    computed_nt = 1,
    pp = NULL,
    ww = NULL,
    wwnorm = NULL,
    tt = NULL,
    residXX = XXwotNA
  )

  out <- extract_component_pls_glm(
    res = res,
    tempww = c(2, 1),
    XXwotNA = XXwotNA,
    XXNA = XXNA,
    PredictYNA = XXNA,
    PredictYwotNA = XXwotNA,
    PredYisdataX = TRUE,
    na.miss.X = FALSE,
    na.miss.Y = FALSE,
    na.miss.PredictY = FALSE,
    tol_Xi = 1e-12,
    kk = 1,
    sparse = FALSE,
    break_nt_sparse = FALSE,
    break_nt_sparse1 = FALSE,
    alpha.pvals.expli = 0.05,
    verbose = FALSE
  )

  expect_false(out$should_break)
  expect_equal(drop(crossprod(out$tempwwnorm)), 1)
  expect_equal(dim(out$res$ww), c(2, 1))
  expect_equal(dim(out$res$wwnorm), c(2, 1))
  expect_equal(dim(out$res$tt), c(4, 1))
  expect_equal(dim(out$res$pp), c(2, 1))
  expect_equal(out$res$residXX, XXwotNA - out$temptt %*% out$temppp)
})

test_that("pls_glm_extract_component stops after sparse-stop signal", {
  XXwotNA <- matrix(1:8, ncol = 2)
  XXNA <- matrix(TRUE, nrow = 4, ncol = 2)
  res <- list(
    nr = 4,
    nc = 2,
    computed_nt = 2,
    pp = matrix(c(1, 0), ncol = 1),
    ww = matrix(c(1, 0), ncol = 1),
    wwnorm = matrix(c(1, 0), ncol = 1),
    tt = matrix(1:4, ncol = 1),
    residXX = XXwotNA
  )

  out <- extract_component_pls_glm(
    res = res,
    tempww = c(1, 1),
    XXwotNA = XXwotNA,
    XXNA = XXNA,
    PredictYNA = XXNA,
    PredictYwotNA = XXwotNA,
    PredYisdataX = TRUE,
    na.miss.X = FALSE,
    na.miss.Y = FALSE,
    na.miss.PredictY = FALSE,
    tol_Xi = 1e-12,
    kk = 2,
    sparse = TRUE,
    break_nt_sparse = TRUE,
    break_nt_sparse1 = TRUE,
    alpha.pvals.expli = 0.05,
    verbose = FALSE
  )

  expect_true(out$should_break)
  expect_equal(out$res$computed_nt, 1)
})
