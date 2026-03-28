#' @rdname plsRmulti
#' @method plsRmultiModel formula
#' @export
plsRmultiModel.formula <- function(object, data, nt = 2, limQ2set = .0975,
                                   modele = "pls", family = NULL,
                                   typeVC = "none", EstimXNA = FALSE,
                                   scaleX = TRUE, scaleY = NULL,
                                   pvals.expli = FALSE,
                                   alpha.pvals.expli = .05,
                                   MClassed = FALSE, tol_Xi = 10^(-12),
                                   weights = NULL, subset = NULL, contrasts = NULL,
                                   sparse = FALSE, sparseStop = FALSE,
                                   naive = FALSE, verbose = TRUE, ...) {
  dots <- list(...)
  pls_multi_check_dots(dots)
  pls_multi_reject_unsupported(
    modele = modele,
    family = family,
    typeVC = typeVC,
    weights_supplied = !missing(weights) && !is.null(weights),
    limQ2set = limQ2set,
    dataPredictY_supplied = FALSE,
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
  if (missing(data)) {
    data <- environment(object)
  }

  built <- pls_multi_build_matrix_from_formula(
    formula = object,
    data = data,
    subset = subset,
    weights = if (missing(weights)) NULL else weights,
    contrasts = contrasts
  )

  estmodel <- pls_multi_fit_core(
    dataY = built$dataY,
    dataX = built$dataX,
    nt = nt,
    scaleX = scaleX,
    scaleY = scaleY,
    tol_Xi = tol_Xi,
    verbose = verbose
  )

  callf0 <- match.call(expand.dots = FALSE)
  callf0$... <- NULL
  callf0$dataY <- callf0$object
  call0 <- c(toString(callf0[[1]]), names(callf0))
  call1 <- call0[!(call0 == "") & !(call0 == "object")]
  estmodel$call <- callf0[call1]
  estmodel$call[[1L]] <- as.name(toString(callf0[[1]]))
  estmodel$call$formula <- object
  estmodel$scaleX <- scaleX
  estmodel$scaleY <- scaleY
  estmodel$terms_x <- built$terms_x
  estmodel$contrasts <- built$contrasts
  estmodel$xlevels <- built$xlevels

  class(estmodel) <- "plsRmultiModel"
  estmodel
}
