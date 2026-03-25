test_that("PLS_glm and PLS_glm_formula share weighted preprocessing behavior", {
  weights <- c(1, 2, 1, 3)
  dat <- data.frame(
    y = c(1, 3, 5, 7),
    x1 = c(1, 2, 4, 8),
    x2 = c(2, 6, 8, 10)
  )

  direct <- PLS_glm(
    dataY = dat$y,
    dataX = dat[c("x1", "x2")],
    nt = 1,
    modele = "pls",
    weights = weights,
    verbose = FALSE
  )
  formula <- PLS_glm_formula(
    y ~ .,
    data = data.frame(
      y = c(1, 3, 5, 7),
      x1 = c(1, 2, 4, 8),
      x2 = c(2, 6, 8, 10)
    ),
    nt = 1,
    modele = "pls",
    weights = c(1, 2, 1, 3),
    verbose = FALSE
  )

  expect_equal(direct$weights, formula$weights)
  expect_equal(unname(direct$RepY), unname(formula$RepY))
  expect_equal(unname(as.matrix(direct$ExpliX)), unname(as.matrix(formula$ExpliX)))
  expect_equal(unname(as.matrix(direct$PredictY)), unname(as.matrix(formula$PredictY)))
  expect_equal(unname(as.matrix(direct$XXNA)), unname(as.matrix(formula$XXNA)))
  expect_equal(unname(direct$YNA), unname(formula$YNA))
  expect_equal(attr(direct$RepY, "scaled:center"), attr(formula$RepY, "scaled:center"))
  expect_equal(attr(direct$RepY, "scaled:scale"), attr(formula$RepY, "scaled:scale"))
  expect_equal(attr(direct$ExpliX, "scaled:center"), attr(formula$ExpliX, "scaled:center"))
  expect_equal(attr(direct$ExpliX, "scaled:scale"), attr(formula$ExpliX, "scaled:scale"))
})

test_that("PLS_glm and PLS_glm_formula share missing-data and prediction preparation", {
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x1 = c(1, NA, 3, 4),
    x2 = c(2, 3, 4, 5)
  )
  pred <- data.frame(
    x1 = c(2, NA),
    x2 = c(3, 6)
  )

  direct <- PLS_glm(
    dataY = dat$y,
    dataX = dat[c("x1", "x2")],
    nt = 1,
    dataPredictY = pred,
    modele = "pls",
    verbose = FALSE
  )
  formula <- PLS_glm_formula(
    y ~ .,
    data = data.frame(
      y = c(1, 2, 3, 4),
      x1 = c(1, NA, 3, 4),
      x2 = c(2, 3, 4, 5)
    ),
    nt = 1,
    dataPredictY = pred,
    modele = "pls",
    verbose = FALSE
  )

  expect_true(direct$na.miss.X)
  expect_true(formula$na.miss.X)
  expect_equal(direct$na.miss.X, formula$na.miss.X)
  expect_equal(direct$na.miss.Y, formula$na.miss.Y)
  expect_equal(direct$XXNA, formula$XXNA)
  expect_equal(direct$ExpliX, formula$ExpliX)
  expect_equal(direct$PredictY, formula$PredictY)
  expect_equal(direct$ttPredictY, formula$ttPredictY)
})

test_that("PLS_glm and PLS_glm_formula share mode selection and family setup", {
  dat <- data.frame(
    y = c(0, 0, 1, 1, 0, 1, 0, 1),
    x1 = c(1, 2, 1, 2, 3, 4, 3, 4),
    x2 = c(1, 1, 2, 2, 3, 3, 4, 4)
  )

  direct <- PLS_glm(
    dataY = dat$y,
    dataX = dat[c("x1", "x2")],
    nt = 1,
    modele = NULL,
    family = stats::binomial(),
    verbose = FALSE
  )
  formula <- PLS_glm_formula(
    y ~ .,
    data = data.frame(
      y = c(0, 0, 1, 1, 0, 1, 0, 1),
      x1 = c(1, 2, 1, 2, 3, 4, 3, 4),
      x2 = c(1, 1, 2, 2, 3, 3, 4, 4)
    ),
    nt = 1,
    modele = NULL,
    family = stats::binomial(),
    verbose = FALSE
  )

  expect_identical(direct$family$family, "binomial")
  expect_identical(formula$family$family, "binomial")
  expect_identical(direct$family$link, formula$family$link)
  expect_equal(direct$ExpliX, formula$ExpliX)
  expect_equal(direct$PredictY, formula$PredictY)
  expect_equal(direct$computed_nt, formula$computed_nt)
})
