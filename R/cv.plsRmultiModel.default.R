#' @rdname cv.plsRmulti
#' @method cv.plsRmultiModel default
#' @export
cv.plsRmultiModel.default <- function(object, dataX, nt = 2, limQ2set = .0975,
                                      modele = "pls", family = NULL,
                                      K = 5, NK = 1, grouplist = NULL,
                                      random = TRUE, scaleX = TRUE,
                                      scaleY = NULL, keepcoeffs = FALSE,
                                      keepfolds = FALSE, keepdataY = TRUE,
                                      keepMclassed = FALSE,
                                      EstimXNA = FALSE,
                                      pvals.expli = FALSE,
                                      alpha.pvals.expli = .05,
                                      MClassed = FALSE,
                                      tol_Xi = 10^(-12), weights,
                                      sparse = FALSE, sparseStop = FALSE,
                                      naive = FALSE, verbose = TRUE, ...) {
  dots <- list(...)
  pls_multi_check_dots(dots)
  pls_multi_reject_cv_unsupported(
    modele = modele,
    family = family,
    weights_supplied = !missing(weights) && !is.null(weights),
    keepMclassed = keepMclassed,
    EstimXNA = EstimXNA,
    pvals.expli = pvals.expli,
    alpha.pvals.expli = alpha.pvals.expli,
    MClassed = MClassed,
    sparse = sparse,
    sparseStop = sparseStop,
    naive = naive
  )

  if (is.null(scaleY)) {
    scaleY <- TRUE
  }

  dataY <- pls_multi_prepare_response(object, tol_Xi = tol_Xi)
  dataX <- pls_multi_prepare_predictors(dataX, n_expected = nrow(dataY), tol_Xi = tol_Xi)
  nt <- pls_multi_validate_nt(nt)
  K <- pls_multi_validate_positive_integer(K, "K")
  NK <- pls_multi_validate_positive_integer(NK, "NK")
  limQ2set <- pls_multi_validate_limited_numeric(limQ2set, "limQ2set")

  if (K > nrow(dataX)) {
    K <- nrow(dataX)
    random <- FALSE
  }

  groups_all <- pls_multi_make_cv_groups(
    n = nrow(dataX),
    K = K,
    NK = NK,
    random = isTRUE(random),
    grouplist = grouplist
  )

  results_kfolds <- vector("list", NK)
  dataY_kfolds <- vector("list", NK)
  coeffs_kfolds <- if (keepcoeffs) vector("list", NK) else NULL
  folds_kfolds <- if (keepfolds) vector("list", NK) else NULL

  for (nnkk in seq_len(NK)) {
    results_kfolds[[nnkk]] <- vector("list", K)
    dataY_kfolds[[nnkk]] <- vector("list", K)
    if (keepcoeffs) {
      coeffs_kfolds[[nnkk]] <- vector("list", K)
    }
    if (keepfolds) {
      folds_kfolds[[nnkk]] <- vector("list", K)
    }

    for (ii in seq_len(K)) {
      fold_res <- pls_multi_cv_fit_fold(
        dataY = dataY,
        dataX = dataX,
        holdout_idx = groups_all[[nnkk]][[ii]],
        nt = nt,
        scaleX = scaleX,
        scaleY = scaleY,
        tol_Xi = tol_Xi,
        keepcoeffs = keepcoeffs,
        verbose = verbose
      )

      results_kfolds[[nnkk]][[ii]] <- fold_res$predictions
      dataY_kfolds[[nnkk]][[ii]] <- fold_res$observed
      if (keepcoeffs) {
        coeffs_kfolds[[nnkk]][[ii]] <- fold_res$coeffs
      }
      if (keepfolds) {
        folds_kfolds[[nnkk]][[ii]] <- fold_res$train_idx
      }
    }
  }

  reference_fit <- plsRmulti(
    object = dataY,
    dataX = dataX,
    nt = nt,
    scaleX = scaleX,
    scaleY = scaleY,
    tol_Xi = tol_Xi,
    verbose = verbose
  )

  cvmodel <- list(
    results_kfolds = results_kfolds,
    dataY_kfolds = dataY_kfolds,
    coeffs_kfolds = coeffs_kfolds,
    folds = folds_kfolds,
    call = NULL,
    reference_fit = reference_fit,
    limQ2set = limQ2set,
    scaleX = scaleX,
    scaleY = scaleY,
    nt = nt,
    ny = reference_fit$ny,
    response_names = reference_fit$response_names,
    predictor_names = reference_fit$predictor_names
  )

  callf0 <- match.call(expand.dots = FALSE)
  callf0$... <- NULL
  callf0$dataY <- callf0$object
  call0 <- c(toString(callf0[[1]]), names(callf0))
  call1 <- call0[!(call0 == "") & !(call0 == "object")]
  cvmodel$call <- callf0[call1]
  cvmodel$call[[1L]] <- as.name(toString(callf0[[1]]))

  class(cvmodel) <- "cv.plsRmultiModel"
  cvmodel
}
