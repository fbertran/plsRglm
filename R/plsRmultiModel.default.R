#' @rdname plsRmulti
#' @method plsRmultiModel default
#' @export
plsRmultiModel.default <- function(object, dataX, nt = 2, limQ2set = .0975,
                                   dataPredictY, modele = "pls",
                                   family = NULL, typeVC = "none",
                                   EstimXNA = FALSE, scaleX = TRUE,
                                   scaleY = NULL, pvals.expli = FALSE,
                                   alpha.pvals.expli = .05,
                                   MClassed = FALSE, tol_Xi = 10^(-12),
                                   weights, sparse = FALSE,
                                   sparseStop = FALSE, naive = FALSE,
                                   verbose = TRUE, ...) {
  dots <- list(...)
  pls_multi_check_dots(dots)
  pls_multi_reject_unsupported(
    modele = modele,
    family = family,
    typeVC = typeVC,
    weights_supplied = !missing(weights) && !is.null(weights),
    limQ2set = limQ2set,
    dataPredictY_supplied = !missing(dataPredictY) && !is.null(dataPredictY),
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

  estmodel <- pls_multi_fit_core(
    dataY = object,
    dataX = dataX,
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
  estmodel$scaleX <- scaleX
  estmodel$scaleY <- scaleY

  class(estmodel) <- "plsRmultiModel"
  estmodel
}
