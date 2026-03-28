#' Experimental multivariate-response PLS2 models
#'
#' This function implements an experimental complete-case linear PLS2 workflow
#' for multivariate numeric responses. It is intentionally separate from
#' \code{\link{plsR}} so the current PLS1 API remains unchanged.
#'
#' @name plsRmulti
#' @aliases plsRmulti plsRmultiModel.default plsRmultiModel.formula
#' @param object a multivariate response matrix/data frame for the default
#' method, or a model formula such as \code{cbind(y1, y2) ~ .} for the formula
#' method.
#' @param ... additional arguments passed to the selected method.
#' @return An object of class \code{"plsRmultiModel"}.
#' @export
plsRmulti <- function(object, ...) UseMethod("plsRmultiModel")

#' @rdname plsRmulti
#' @aliases plsRmulti
#' @export plsRmultiModel
plsRmultiModel <- plsRmulti
