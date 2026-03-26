test_that("classbiplot colors grouped individuals on a biplot", {
  data(aze)

  fit <- plsRglm(y ~ ., data = aze, nt = 2, modele = "pls-glm-logistic")

  outfile <- tempfile(fileext = ".pdf")
  grDevices::pdf(outfile)
  dev_id <- grDevices::dev.cur()
  on.exit({
    devs <- grDevices::dev.list()
    if (!is.null(devs) && dev_id %in% devs) {
      grDevices::dev.off(which = dev_id)
    }
    unlink(outfile)
  }, add = TRUE)

  res <- classbiplot(fit, group = aze$y,
                     col = c("firebrick3", "steelblue3"))

  expect_type(res, "list")
  expect_equal(names(res),
               c("scores", "loadings", "colours", "group", "ratio"))
  expect_equal(dim(res$scores), c(nrow(aze), 2L))
  expect_equal(ncol(res$loadings), 2L)
  expect_equal(length(res$colours), nrow(aze))
  expect_true(is.factor(res$group))
  expect_true(is.numeric(res$ratio))
  expect_gt(res$ratio, 0)
})

test_that("classbiplot validates the color specification", {
  data(Cornell)

  fit <- plsR(Y ~ ., data = Cornell, nt = 2)
  grp <- factor(Cornell$Y > median(Cornell$Y))

  expect_error(
    classbiplot(fit, group = grp, col = c("red", "blue", "green")),
    "'col' should have length 1, the number of groups, or the number of individuals."
  )
})
