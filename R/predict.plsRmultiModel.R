#' Predict method for experimental plsRmulti models
#'
#' @param object An object of class \code{"plsRmultiModel"}.
#' @param newdata Optional predictor values for out-of-sample prediction.
#' @param comps Number of components to use for prediction.
#' @param type Either \code{"response"} or \code{"scores"}.
#' @param verbose Should informational messages be displayed?
#' @param ... not used.
#' @return A matrix of predicted responses or score coordinates.
#' @export
predict.plsRmultiModel <- function(object, newdata,
                                   comps = object$computed_nt,
                                   type = c("response", "scores"),
                                   verbose = TRUE, ...) {
  dots <- list(...)
  pls_multi_check_dots(dots)

  if (!inherits(object, "plsRmultiModel")) {
    stop("Primary argument much be a plsRmultiModel object")
  }

  type <- match.arg(type)
  if (length(comps) != 1L || is.na(comps)) {
    pls_multi_stop("comps must be a single positive integer")
  }
  comps <- as.integer(comps)
  if (comps < 1L || comps > object$computed_nt) {
    pls_multi_stop("Cannot predict using more components than extracted")
  }

  if (missing(newdata) || is.null(newdata)) {
    if (type == "scores") {
      return(object$tt[, seq_len(comps), drop = FALSE])
    }

    fitted_scaled <- pls_multi_predict_response_scaled(
      x_scaled = object$ExpliX,
      object = object,
      comps = comps
    )
    return(pls_multi_backtransform_response(fitted_scaled, object))
  }

  if (verbose && anyNA(newdata)) {
    pls_multi_stop("plsRmulti predict() does not support missing values in this experimental release")
  }

  newdata_matrix <- pls_multi_prepare_newdata(object, newdata)
  scaled_newdata <- pls_multi_scale_newdata(newdata_matrix, object)
  tt_new <- pls_multi_predict_scores_matrix(
    x_scaled = scaled_newdata,
    object = object,
    comps = comps
  )
  colnames(tt_new) <- paste("Comp_", seq_len(comps), sep = "")

  if (type == "scores") {
    return(tt_new)
  }

  pred_scaled <- pls_multi_predict_response_scaled(
    x_scaled = scaled_newdata,
    object = object,
    comps = comps
  )
  pred <- pls_multi_backtransform_response(pred_scaled, object)
  rownames(pred) <- rownames(newdata_matrix)
  pred
}
