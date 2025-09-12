test_that("exported objects exist in the namespace", {
  # Get exports from the current namespace (works for source or installed)
  ns <- asNamespace("plsRglm")
  ex <- getNamespaceExports("plsRglm")
  expect_gt(length(ex), 0L)
  missing <- character(0)
  for (sym in ex) {
    if (!exists(sym, envir = ns, inherits = FALSE)) missing <- c(missing, sym)
  }
  if (length(missing)) {
    fail(sprintf("Missing exported symbols: %s", paste(missing, collapse = ", ")))
  } else {
    succeed()
  }
})
