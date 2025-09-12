test_that("bundled datasets load with expected shape and types", {
  # bordeaux & XbordeauxNA
  data(bordeaux, package = "plsRglm")
  expect_s3_class(bordeaux, "data.frame")
  expect_equal(dim(bordeaux), c(34L, 5L))
  expect_true(is.ordered(bordeaux$Quality))
  expect_equal(levels(bordeaux$Quality), c("1","2","3"))
  expect_true(is.numeric(bordeaux$Temperature))
  expect_true(is.numeric(bordeaux$Sunshine))
  expect_true(is.numeric(bordeaux$Heat))
  expect_true(is.numeric(bordeaux$Rain))

  data(XbordeauxNA, package = "plsRglm")
  expect_equal(dim(XbordeauxNA), c(34L, 4L))
  # Documentation mentions the first Temperature was removed on purpose
  expect_true(is.na(XbordeauxNA$Temperature[1]))
  expect_true(sum(is.na(XbordeauxNA$Temperature)) >= 1)

  # Pine family
  data(pine, package = "plsRglm")
  expect_equal(dim(pine), c(33L, 11L))
  expect_true(all(startsWith(names(pine), "x")))
  expect_true(is.numeric(pine$x11))

  data(pine_sup, package = "plsRglm")
  expect_equal(dim(pine_sup), c(25L, 11L))

  data(pineNAX21, package = "plsRglm")
  expect_equal(dim(pineNAX21), c(33L, 11L))
  expect_true(sum(is.na(pineNAX21)) >= 1)

  data(XpineNAX21, package = "plsRglm")
  expect_equal(dim(XpineNAX21), c(33L, 10L))
  expect_true(sum(is.na(XpineNAX21)) >= 1)

  # Cornell
  data(Cornell, package = "plsRglm")
  expect_equal(dim(Cornell), c(12L, 8L))
  expect_true(is.numeric(Cornell$Y))

  # CorMat
  data(CorMat, package = "plsRglm")
  expect_equal(dim(CorMat), c(17L, 17L))
  expect_true(all(vapply(CorMat, is.numeric, logical(1))))
})
