test_that("package provides a citation", {
  cit <- utils::citation("plsRglm")
  # citation() returns a 'bibentry' or list of bibentry
  expect_true(inherits(cit, "bibentry") || (is.list(cit) && all(vapply(cit, inherits, logical(1), "bibentry"))))
})
