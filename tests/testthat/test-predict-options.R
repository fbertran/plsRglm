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
