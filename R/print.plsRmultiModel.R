#' @export
print.plsRmultiModel <- function(x, ...) {
  cat("Experimental multivariate PLS2 model\n")
  cat("Number of responses:\n")
  print(x$ny)
  cat("Number of required components:\n")
  print(x$nt)
  cat("Number of successfully computed components:\n")
  print(x$computed_nt)
  cat("Coefficient matrix:\n")
  print(x$Coeffs)
  invisible(x)
}
