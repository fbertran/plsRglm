test_that("PLS (lm) fits on Cornell and infcrit.dof returns expected structure", {
  skip_on_cran()
  data(Cornell, package = "plsRglm")
  X <- Cornell[,1:7]
  y <- Cornell$Y
  # Fit 2 comps for speed
  set.seed(1)
  mod <- plsR(y, X, nt = 2, verbose = FALSE)
  expect_s3_class(mod, "plsRmodel")
  pv <- predict(mod)
  expect_equal(length(pv), length(y))

  ic <- infcrit.dof(mod)
  expect_true(is.matrix(ic))
  expect_true(all(c("AIC.dof","BIC.dof","GMDL.dof","AIC.naive","BIC.naive","GMDL.naive") %in% colnames(ic)))
  # rows are Nb_Comp_0..nt
  expect_equal(rownames(ic)[1], "Nb_Comp_0")
})
