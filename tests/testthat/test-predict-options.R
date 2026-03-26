test_that("predict.plsRglmmodel supports link/response/scores", {
  skip_on_cran()
  data(aze, package = "plsRglm")
  set.seed(42)
  fit <- plsRglm(y ~ ., data = aze, nt = 2, modele = "pls-glm-family",
                 family = stats::binomial(), verbose = FALSE)

  # link & response
  eta <- predict(fit, type = "link")
  mu  <- predict(fit, type = "response")
  expect_equal(length(eta), nrow(aze))
  expect_equal(length(mu), nrow(aze))
  expect_true(all(is.finite(eta)))
  expect_true(all(mu >= 0 & mu <= 1))

  # scores
  sc <- predict(fit, type = "scores")
  expect_true(is.matrix(sc) || is.data.frame(sc))
  expect_equal(nrow(sc), nrow(aze))
  expect_true(ncol(sc) >= 1)
})

test_that("predict.plsRglmmodel supports pls mode and dataY alias", {
  skip_on_cran()
  data(Cornell, package = "plsRglm")
  X <- Cornell[, 1:7]
  y <- Cornell$Y

  fit <- plsRglm(dataY = y, dataX = X, nt = 2, modele = "pls",
                 verbose = FALSE)

  eta <- predict(fit, type = "link")
  mu <- predict(fit, type = "response")
  mu_new <- predict(fit, newdata = X, type = "response")
  sc <- predict(fit, type = "scores")

  expect_s3_class(fit, "plsRglmmodel")
  expect_true(is.matrix(eta) || is.data.frame(eta))
  expect_true(is.matrix(mu) || is.data.frame(mu))
  expect_equal(nrow(mu), nrow(Cornell))
  expect_equal(unname(drop(mu)), unname(drop(fit$YChapeau)))
  expect_equal(unname(drop(mu_new)), unname(drop(fit$ValsPredictY)))
  expect_true(is.matrix(sc) || is.data.frame(sc))
  expect_equal(nrow(sc), nrow(Cornell))
})
