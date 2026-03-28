#' Summary method for experimental multivariate PLS2 CV models
#'
#' @param object an object of class \code{"cv.plsRmultiModel"}
#' @param verbose should progress information be displayed?
#' @param ... further arguments passed to methods.
#' @return A list of per-partition CV summary tables.
#' @export
summary.cv.plsRmultiModel <- function(object, verbose = TRUE, ...) {
  computed_nt <- min(
    pls_multi_cv_min_components(object),
    object$reference_fit$computed_nt
  )
  press <- pls_multi_cv_press(object, max_comps = computed_nt)

  res <- vector("list", length(object$results_kfolds))
  for (nnkk in seq_along(object$results_kfolds)) {
    if (verbose) {
      if (nnkk %% 10 == 1) {
        cat("\n")
        cat(paste("NK:", nnkk))
      } else {
        cat(paste(", ", nnkk))
      }
    }
    res[[nnkk]] <- pls_multi_cv_summary_matrix(
      press_total = press$total[[nnkk]],
      press_by_response = press$by_response[[nnkk]],
      reference_fit = object$reference_fit,
      limQ2set = object$limQ2set
    )
  }

  if (verbose) {
    cat("\n")
  }

  class(res) <- c("summary.cv.plsRmultiModel", "summary.cv.plsRmodel")
  res
}
