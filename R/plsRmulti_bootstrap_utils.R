coefs.plsRmulti <- function(dataset, ind, nt, ny, scaleX, scaleY, tol_Xi,
                            maxcoefvalues, ifbootfail, response_names,
                            predictor_names, verbose, ...) {
  tempfit <- try(
    plsRmulti(
      object = dataset[ind, seq_len(ny), drop = FALSE],
      dataX = dataset[ind, -(seq_len(ny)), drop = FALSE],
      nt = nt,
      scaleX = scaleX,
      scaleY = scaleY,
      tol_Xi = tol_Xi,
      verbose = verbose
    ),
    silent = TRUE
  )

  if (inherits(tempfit, "try-error")) {
    return(ifbootfail)
  }

  tempcoefs <- pls_multi_standardized_boot_vector(tempfit)
  cond <- FALSE
  try(cond <- is.numeric(tempcoefs) && all(abs(tempcoefs) < maxcoefvalues), silent = TRUE)
  if (cond) {
    return(tempcoefs)
  }

  ifbootfail
}

permcoefs.plsRmulti <- function(...) {
  coefs.plsRmulti(...)
}

coefs.plsRmulti_np <- function(dataRepYtt, ind, nt, ny, maxcoefvalues,
                               wwetoile, ifbootfail, response_names,
                               predictor_names, verbose = TRUE, ...) {
  y_boot <- as.matrix(dataRepYtt[ind, seq_len(ny), drop = FALSE])
  t_boot <- as.matrix(dataRepYtt[ind, ny + seq_len(nt), drop = FALSE])
  tempcb <- try(
    solve(crossprod(t_boot), crossprod(t_boot, y_boot)),
    silent = TRUE
  )

  if (inherits(tempcb, "try-error")) {
    return(ifbootfail)
  }

  tempcoefs <- pls_multi_vectorize_coefficients(
    coef_matrix = wwetoile[, seq_len(nt), drop = FALSE] %*% tempcb,
    intercept = rep(0, ny),
    response_names = response_names,
    predictor_names = predictor_names
  )
  cond <- FALSE
  try(cond <- is.numeric(tempcoefs) && all(abs(tempcoefs) < maxcoefvalues), silent = TRUE)
  if (cond) {
    return(tempcoefs)
  }

  ifbootfail
}

permcoefs.plsRmulti_np <- function(...) {
  coefs.plsRmulti_np(...)
}

bootpls_multi <- function(object, typeboot = "plsmodel", R = 250,
                          statistic = NULL, sim = "ordinary", stype = "i",
                          stabvalue = 1e6, verbose = TRUE, ...) {
  nt <- object$computed_nt
  tol_Xi <- kfolds_resolve_numeric_arg(object$call$tol_Xi, default = 10^(-12))
  base_vector <- pls_multi_standardized_boot_vector(object)
  maxcoefvalues <- pls_multi_boot_thresholds(base_vector, stabvalue = stabvalue)
  ifbootfail <- matrix(NA_real_, nrow = nrow(base_vector), ncol = 1L,
                       dimnames = dimnames(base_vector))

  if (typeboot == "plsmodel") {
    dataset <- cbind(object$dataY, object$dataX)
    colnames(dataset) <- c(object$response_names, object$predictor_names)
    if (is.null(statistic)) {
      statistic <- if (identical(sim, "permutation")) permcoefs.plsRmulti else coefs.plsRmulti
    }
    temp.bootplsR <- boot(
      data = dataset,
      statistic = statistic,
      sim = sim,
      stype = stype,
      R = R,
      nt = nt,
      ny = object$ny,
      scaleX = object$scaleX,
      scaleY = object$scaleY,
      tol_Xi = tol_Xi,
      maxcoefvalues = maxcoefvalues,
      ifbootfail = ifbootfail,
      response_names = object$response_names,
      predictor_names = object$predictor_names,
      verbose = verbose,
      ...
    )
  } else if (typeboot == "fmodel_np") {
    dataRepYtt <- cbind(object$RepY, object$tt)
    colnames(dataRepYtt) <- c(object$response_names, colnames(object$tt))
    if (is.null(statistic)) {
      statistic <- if (identical(sim, "permutation")) permcoefs.plsRmulti_np else coefs.plsRmulti_np
    }
    temp.bootplsR <- boot(
      data = dataRepYtt,
      statistic = statistic,
      sim = sim,
      stype = stype,
      R = R,
      nt = nt,
      ny = object$ny,
      maxcoefvalues = maxcoefvalues,
      wwetoile = object$wwetoile,
      ifbootfail = ifbootfail,
      response_names = object$response_names,
      predictor_names = object$predictor_names,
      verbose = verbose,
      ...
    )
  } else if (typeboot == "fmodel_par") {
    dataset <- cbind(object$dataY, object$dataX)
    colnames(dataset) <- c(object$response_names, object$predictor_names)
    if (is.null(statistic)) {
      statistic <- if (identical(sim, "permutation")) permcoefs.plsRmulti else coefs.plsRmulti
    }
    temp.bootplsR <- boot(
      data = dataset,
      statistic = statistic,
      sim = sim,
      stype = stype,
      R = R,
      nt = nt,
      ny = object$ny,
      scaleX = object$scaleX,
      scaleY = object$scaleY,
      tol_Xi = tol_Xi,
      maxcoefvalues = maxcoefvalues,
      ifbootfail = ifbootfail,
      response_names = object$response_names,
      predictor_names = object$predictor_names,
      verbose = verbose,
      ...
    )
  } else {
    pls_multi_stop("Unsupported typeboot for plsRmulti bootstrap")
  }

  valid_rows <- stats::complete.cases(temp.bootplsR$t)
  temp.bootplsR$t <- temp.bootplsR$t[valid_rows, , drop = FALSE]
  temp.bootplsR$R <- sum(valid_rows)
  temp.bootplsR$call$R <- sum(valid_rows)
  temp.bootplsR
}
