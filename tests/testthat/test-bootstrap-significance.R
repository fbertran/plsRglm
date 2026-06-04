test_that("confints.bootpls labels CI columns and preserves predictor names", {
  skip_on_cran()
  data(Cornell, package = "plsRglm")

  set.seed(1)
  mod <- plsR(Y ~ ., data = Cornell, nt = 1, verbose = FALSE)
  bt <- bootpls(mod, R = 20, verbose = FALSE)
  ci <- suppressWarnings(confints.bootpls(bt, indices = 2:4, typeBCa = FALSE))

  expect_equal(
    colnames(ci),
    c("Normal.Lower", "Normal.Upper",
      "Basic.Lower", "Basic.Upper",
      "Percentile.Lower", "Percentile.Upper")
  )
  expect_equal(rownames(ci), c("X1", "X2", "X3"))
  expect_false(isTRUE(attr(ci, "typeBCa")))
})

test_that("confints2signifind defaults to BCa when available", {
  ci <- matrix(
    c(-1, 1, -1, 1, -1, 1, 0.2, 0.8,
      -1, 1, -1, 1, -1, 1, -0.4, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("X1", "X2"), NULL)
  )
  colnames(ci) <- c("Normal.Lower", "Normal.Upper",
                    "Basic.Lower", "Basic.Upper",
                    "Percentile.Lower", "Percentile.Upper",
                    "BCa.Lower", "BCa.Upper")
  attr(ci, "typeBCa") <- TRUE

  ind <- confints2signifind(ci)

  expect_identical(ind, c(X1 = TRUE, X2 = FALSE))
})

test_that("confints2signifind falls back to percentile intervals when BCa is absent", {
  ci <- matrix(
    c(-0.7, 0.2, -0.4, 0.1, 0.2, 0.9,
      -0.4, 0.5, -0.3, 0.4, -0.2, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("X1", "X2"), NULL)
  )
  colnames(ci) <- c("Normal.Lower", "Normal.Upper",
                    "Basic.Lower", "Basic.Upper",
                    "Percentile.Lower", "Percentile.Upper")
  attr(ci, "typeBCa") <- FALSE

  ind <- confints2signifind(ci)

  expect_identical(ind, c(X1 = TRUE, X2 = FALSE))
})

test_that("confints2signifind rejects BCa when those limits are unavailable", {
  ci <- matrix(
    c(-0.7, 0.2, -0.4, 0.1, 0.2, 0.9),
    nrow = 1,
    dimnames = list("X1", c("Normal.Lower", "Normal.Upper",
                             "Basic.Lower", "Basic.Upper",
                             "Percentile.Lower", "Percentile.Upper"))
  )
  attr(ci, "typeBCa") <- FALSE

  expect_error(
    confints2signifind(ci, typeIC = "BCa"),
    "BCa intervals were not computed"
  )
})
