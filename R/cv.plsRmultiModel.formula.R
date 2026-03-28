#' @rdname cv.plsRmulti
#' @method cv.plsRmultiModel formula
#' @export
cv.plsRmultiModel.formula <- function(object, data = NULL, nt = 2,
                                      limQ2set = .0975, modele = "pls",
                                      family = NULL, K = 5, NK = 1,
                                      grouplist = NULL, random = TRUE,
                                      scaleX = TRUE, scaleY = NULL,
                                      keepcoeffs = FALSE, keepfolds = FALSE,
                                      keepdataY = TRUE, keepMclassed = FALSE,
                                      EstimXNA = FALSE,
                                      pvals.expli = FALSE,
                                      alpha.pvals.expli = .05,
                                      MClassed = FALSE,
                                      tol_Xi = 10^(-12), weights = NULL,
                                      subset = NULL, contrasts = NULL,
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
  if (missing(data) || is.null(data)) {
    data <- environment(object)
  }

  built <- pls_multi_build_matrix_from_formula(
    formula = object,
    data = data,
    subset = subset,
    weights = if (missing(weights)) NULL else weights,
    contrasts = contrasts
  )

  cvmodel <- cv.plsRmulti(
    object = built$dataY,
    dataX = built$dataX,
    nt = nt,
    limQ2set = limQ2set,
    modele = modele,
    family = family,
    K = K,
    NK = NK,
    grouplist = grouplist,
    random = random,
    scaleX = scaleX,
    scaleY = scaleY,
    keepcoeffs = keepcoeffs,
    keepfolds = keepfolds,
    keepdataY = keepdataY,
    keepMclassed = keepMclassed,
    EstimXNA = EstimXNA,
    pvals.expli = pvals.expli,
    alpha.pvals.expli = alpha.pvals.expli,
    MClassed = MClassed,
    tol_Xi = tol_Xi,
    weights = weights,
    sparse = sparse,
    sparseStop = sparseStop,
    naive = naive,
    verbose = verbose
  )

  callf0 <- match.call(expand.dots = FALSE)
  callf0$... <- NULL
  callf0$formula <- object
  call0 <- c(toString(callf0[[1]]), names(callf0))
  call1 <- call0[!(call0 == "") & !(call0 == "object")]
  cvmodel$call <- callf0[call1]
  cvmodel$call[[1L]] <- as.name(toString(callf0[[1]]))
  cvmodel$call$formula <- object
  cvmodel$formula <- object
  cvmodel$terms_x <- built$terms_x
  cvmodel$contrasts <- built$contrasts
  cvmodel$xlevels <- built$xlevels

  class(cvmodel) <- "cv.plsRmultiModel"
  cvmodel
}
