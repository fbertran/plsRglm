make_pls_multi_fixture <- function() {
  set.seed(42)
  X <- matrix(rnorm(90 * 5), nrow = 90, ncol = 5)
  colnames(X) <- paste0("x", 1:5)

  B <- matrix(
    c(
      1.2, -0.4, 0.3,
      -0.8, 0.7, 0.5,
      0.5, 1.1, -0.6,
      0.0, 0.4, 0.9,
      -0.3, 0.2, 1.0
    ),
    nrow = 5,
    ncol = 3,
    byrow = TRUE
  )
  Y <- X %*% B + matrix(rnorm(90 * 3, sd = 0.05), nrow = 90, ncol = 3)
  colnames(Y) <- paste0("y", 1:3)

  df <- data.frame(Y, X, check.names = FALSE)
  list(X = X, Y = Y, df = df)
}

reference_pls2 <- function(X, Y, nt, scaleX = TRUE, scaleY = TRUE, tol = 1e-12) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (scaleX) {
    x_center <- colMeans(X)
    Xs <- scale(X)
  } else {
    x_center <- rep(0, ncol(X))
    Xs <- X
    attr(Xs, "scaled:center") <- x_center
    attr(Xs, "scaled:scale") <- rep(1, ncol(X))
  }

  if (scaleY) {
    y_center <- colMeans(Y)
    Ys <- scale(Y)
  } else {
    y_center <- rep(0, ncol(Y))
    Ys <- Y
    attr(Ys, "scaled:center") <- y_center
    attr(Ys, "scaled:scale") <- rep(1, ncol(Y))
  }

  Xres <- Xs
  Yres <- Ys
  W <- matrix(numeric(0), nrow = ncol(Xs), ncol = 0L)
  Tt <- matrix(numeric(0), nrow = nrow(Xs), ncol = 0L)
  P <- matrix(numeric(0), nrow = ncol(Xs), ncol = 0L)
  C <- matrix(numeric(0), nrow = ncol(Ys), ncol = 0L)

  for (k in seq_len(nt)) {
    S <- crossprod(Xres, Yres)
    if (sum(S^2) <= tol) {
      break
    }
    sv <- svd(S, nu = 0L, nv = 1L)
    w_raw <- drop(S %*% sv$v[, 1, drop = FALSE])
    w <- w_raw / sqrt(sum(w_raw^2))
    idx <- which.max(abs(w))
    if (w[idx] < 0) {
      w <- -w
    }
    t_score <- drop(Xres %*% w)
    tt_norm <- drop(crossprod(t_score))
    p <- drop(crossprod(Xres, t_score) / tt_norm)
    c <- drop(crossprod(Yres, t_score) / tt_norm)
    W <- cbind(W, w)
    Tt <- cbind(Tt, t_score)
    P <- cbind(P, p)
    C <- cbind(C, c)
    Xres <- Xres - tcrossprod(t_score, p)
    Yres <- Yres - tcrossprod(t_score, c)
  }

  B_scaled <- W %*% solve(crossprod(P, W)) %*% t(C)
  Yhat_scaled <- Xs %*% B_scaled
  Yhat <- sweep(Yhat_scaled, 2, attr(Ys, "scaled:scale"), "*")
  Yhat <- sweep(Yhat, 2, attr(Ys, "scaled:center"), "+")

  list(tt = Tt, pp = P, CoeffC = C, fitted = Yhat)
}

test_that("plsRmulti matrix interface returns coherent multivariate structures", {
  fixture <- make_pls_multi_fixture()
  fit <- plsRmulti(fixture$Y, fixture$X, nt = 3, verbose = FALSE)

  expect_s3_class(fit, "plsRmultiModel")
  expect_equal(dim(fit$tt), c(nrow(fixture$X), 3))
  expect_equal(dim(fit$pp), c(ncol(fixture$X), 3))
  expect_equal(dim(fit$CoeffC), c(ncol(fixture$Y), 3))
  expect_equal(dim(fit$Coeffs), c(ncol(fixture$X), ncol(fixture$Y)))
  expect_length(fit$CoeffConstante, ncol(fixture$Y))
  expect_equal(dim(fit$YChapeau), dim(fixture$Y))
  expect_equal(dim(predict(fit, type = "response")), dim(fixture$Y))
  expect_equal(dim(predict(fit, type = "scores")), c(nrow(fixture$X), 3))
})

test_that("plsRmulti formula and matrix interfaces agree", {
  fixture <- make_pls_multi_fixture()
  fit_matrix <- plsRmulti(fixture$Y[, 1:2, drop = FALSE], fixture$X, nt = 2, verbose = FALSE)
  fit_formula <- plsRmulti(cbind(y1, y2) ~ ., data = fixture$df[, c("y1", "y2", colnames(fixture$X))], nt = 2, verbose = FALSE)

  expect_equal(unname(fit_formula$tt), unname(fit_matrix$tt), tolerance = 1e-8)
  expect_equal(unname(fit_formula$pp), unname(fit_matrix$pp), tolerance = 1e-8)
  expect_equal(unname(fit_formula$CoeffC), unname(fit_matrix$CoeffC), tolerance = 1e-8)
  expect_equal(unname(fit_formula$Coeffs), unname(fit_matrix$Coeffs), tolerance = 1e-8)
})

test_that("plsRmulti predict supports newdata for response and scores", {
  fixture <- make_pls_multi_fixture()
  fit <- plsRmulti(cbind(y1, y2) ~ ., data = fixture$df[, c("y1", "y2", colnames(fixture$X))], nt = 2, verbose = FALSE)
  newdata <- fixture$df[1:7, colnames(fixture$X)]

  pred_response <- predict(fit, newdata = newdata, type = "response")
  pred_scores <- predict(fit, newdata = newdata, type = "scores")

  expect_equal(dim(pred_response), c(7, 2))
  expect_equal(dim(pred_scores), c(7, 2))
})

test_that("plsRmulti satisfies core PLS2 numerical identities", {
  fixture <- make_pls_multi_fixture()
  fit <- plsRmulti(fixture$Y, fixture$X, nt = 3, verbose = FALSE)

  score_crossprod <- crossprod(fit$tt)
  expect_equal(score_crossprod[upper.tri(score_crossprod)], rep(0, 3), tolerance = 1e-8)

  expect_equal(
    fit$wwetoile %*% t(fit$CoeffC),
    fit$Std.Coeffs,
    tolerance = 1e-8
  )

  fitted_from_original <- sweep(fixture$X %*% fit$Coeffs, 2, fit$CoeffConstante, "+")
  expect_equal(unname(fitted_from_original), unname(fit$YChapeau), tolerance = 1e-8)
})

test_that("plsRmulti matches an independent complete-case reference implementation", {
  fixture <- make_pls_multi_fixture()
  fit <- plsRmulti(fixture$Y, fixture$X, nt = 3, verbose = FALSE)
  ref <- reference_pls2(fixture$X, fixture$Y, nt = 3)

  expect_equal(unname(fit$tt), unname(ref$tt), tolerance = 1e-8)
  expect_equal(unname(fit$pp), unname(ref$pp), tolerance = 1e-8)
  expect_equal(unname(fit$CoeffC), unname(ref$CoeffC), tolerance = 1e-8)
  expect_equal(unname(fit$YChapeau), unname(ref$fitted), tolerance = 1e-8)
})

test_that("plsRmulti rejects unsupported and invalid inputs", {
  fixture <- make_pls_multi_fixture()

  expect_error(
    plsRmulti(fixture$Y[, 1], fixture$X, nt = 2, verbose = FALSE),
    "at least two response columns"
  )
  expect_error(
    plsRmulti(data.frame(a = letters[1:nrow(fixture$X)], b = letters[1:nrow(fixture$X)]), fixture$X, nt = 2, verbose = FALSE),
    "numeric multivariate response"
  )
  expect_error(
    plsRmulti(fixture$Y, replace(fixture$X, 1, NA), nt = 2, verbose = FALSE),
    "missing values in X"
  )
  expect_error(
    plsRmulti(replace(fixture$Y, 1, NA), fixture$X, nt = 2, verbose = FALSE),
    "missing values in Y"
  )
  expect_error(
    plsRmulti(fixture$Y, fixture$X, nt = 2, weights = rep(1, nrow(fixture$X)), verbose = FALSE),
    "does not support weights"
  )
  expect_error(
    plsRmulti(fixture$Y, fixture$X, nt = 2, family = gaussian(), verbose = FALSE),
    "does not support GLM families"
  )
  expect_error(
    plsRmulti(fixture$Y, fixture$X, nt = 2, typeVC = "standard", verbose = FALSE),
    "does not support cross-validation"
  )
  expect_error(
    plsRmulti(fixture$Y, fixture$X, nt = 2, dataPredictY = fixture$X, verbose = FALSE),
    "does not support dataPredictY"
  )
})
