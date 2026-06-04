test_that("plsR falls back to naive information criteria when plsdof is unavailable", {
  skip_on_cran()
  local_mocked_bindings(
    plsR_has_plsdof = function() FALSE,
    .package = "plsRglm"
  )

  data(Cornell, package = "plsRglm")
  mod <- plsR(Cornell$Y, Cornell[, 1:7], nt = 2, verbose = FALSE)

  expect_s3_class(mod, "plsRmodel")
  expect_true(all(is.na(mod$ic.dof[, c("DoF.dof", "sigmahat.dof", "AIC.dof", "BIC.dof", "GMDL.dof")])))
  expect_false(any(is.na(mod$ic.dof[, c("DoF.naive", "sigmahat.naive", "AIC.naive", "BIC.naive", "GMDL.naive")])))
})

test_that("cv.plsR summary works when plsdof is unavailable", {
  skip_on_cran()
  local_mocked_bindings(
    plsR_has_plsdof = function() FALSE,
    .package = "plsRglm"
  )

  data(Cornell, package = "plsRglm")
  cvfit <- cv.plsR(Y ~ ., data = Cornell, nt = 2, K = 5, NK = 1, verbose = FALSE)
  s <- summary(cvfit, verbose = FALSE)

  expect_s3_class(cvfit, "cv.plsRmodel")
  expect_s3_class(s, "summary.cv.plsRmodel")
  expect_true(all(is.na(s[[1]][, c("DoF.dof", "sigmahat.dof", "AIC.dof", "BIC.dof", "GMDL.dof")])))
  expect_false(any(is.na(s[[1]][, c("DoF.naive", "sigmahat.naive", "AIC.naive", "BIC.naive", "GMDL.naive")])))
})
