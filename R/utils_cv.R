kfolds_requested_nt <- function(pls_kfolds) {
  if (!length(pls_kfolds) || !length(pls_kfolds[[1L]])) {
    return(0L)
  }

  max(vapply(pls_kfolds[[1L]], function(one_nk) {
    if (!length(one_nk)) {
      return(0L)
    }
    max(vapply(one_nk, NCOL, integer(1L)))
  }, integer(1L)))
}

kfolds_resolve_numeric_arg <- function(arg, default = NA_real_, envir = parent.frame()) {
  if (is.null(arg)) {
    return(default)
  }

  value <- try(eval(arg, envir = envir), silent = TRUE)
  if (inherits(value, "try-error") || length(value) != 1L || !is.numeric(value)) {
    return(default)
  }

  as.numeric(value)
}
