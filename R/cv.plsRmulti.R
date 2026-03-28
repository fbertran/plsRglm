#' Experimental k-fold cross-validation for multivariate-response PLS2
#'
#' This function performs repeated k-fold cross-validation for the experimental
#' complete-case linear \code{\link{plsRmulti}} workflow.
#'
#' @name cv.plsRmulti
#' @aliases cv.plsRmulti cv.plsRmultiModel.default cv.plsRmultiModel.formula
#' @param object For the default method, a numeric multivariate response matrix
#' or data frame with at least two columns. For the formula method, a formula
#' of the form \code{cbind(y1, y2, ...) ~ .}.
#' @param ... Additional arguments passed to the selected method.
#' @return An object of class \code{"cv.plsRmultiModel"}.
#' @export
cv.plsRmulti <- function(object, ...) UseMethod("cv.plsRmultiModel")

#' @rdname cv.plsRmulti
#' @aliases cv.plsRmulti
#' @export cv.plsRmultiModel
cv.plsRmultiModel <- cv.plsRmulti
