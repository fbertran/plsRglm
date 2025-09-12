test_that("package loads and has a namespace", {
  expect_silent(library(plsRglm))
  expect_true("plsRglm" %in% loadedNamespaces())
  ns <- asNamespace("plsRglm")
  expect_true(is.environment(ns))
})
