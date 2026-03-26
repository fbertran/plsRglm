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

test_that("summary.cv.plsRmodel handles symbolic nt values", {
  skip_on_cran()
  data(Cornell, package = "plsRglm")

  nt_val <- 5
  cvfit <- cv.plsR(Y ~ ., data = Cornell, nt = nt_val, K = 5, NK = 1,
                   verbose = FALSE)
  rm(nt_val)

  s <- summary.cv.plsRmodel(cvfit, verbose = FALSE)

  expect_s3_class(cvfit, "cv.plsRmodel")
  expect_s3_class(s, "summary.cv.plsRmodel")
  expect_length(s, 1L)
  expect_equal(nrow(s[[1]]), 6L)
})
