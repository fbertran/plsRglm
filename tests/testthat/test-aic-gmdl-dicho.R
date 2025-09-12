test_that("AICpls and info-criteria helpers return numeric values", {
  set.seed(1)
  # Simple sanity check with lm residuals
  y <- rnorm(20)
  x <- rnorm(20)
  r <- residuals(lm(y ~ x))
  expect_type(AICpls(1, r), "double")
  # DoF-based helpers on a small plsR model
  data(Cornell, package = "plsRglm")
  mod <- plsR(Cornell$Y, Cornell[,1:7], nt = 1, verbose = FALSE)
  # infcrit.dof already validated in another test
  ic <- infcrit.dof(mod)
  expect_true(is.matrix(ic))
})

test_that("dicho behaves as expected on numeric input", {
  x <- c(-1, 0, 0.2, 0.8, 1.5)
  y <- dicho(x)
  expect_true(all(y %in% c(0,1)))
  expect_equal(which(y == 1), which(x > 0))
})
