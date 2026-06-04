#' Detect significant predictors from bootstrap confidence intervals
#' 
#' This function converts the matrix returned by \code{\link{confints.bootpls}}
#' into a logical indicator of significance. A predictor is marked as
#' significant when the selected confidence interval lies entirely on one side
#' of the reference value.
#' 
#' @param ic_bootobject a matrix returned by \code{\link{confints.bootpls}}.
#' @param typeIC type of confidence interval to use. Defaults to \code{"BCa"}
#' when BCa limits are available and to \code{"Percentile"} otherwise.
#' @param threshold reference value to test against. Defaults to \code{0}.
#' @return Named logical vector. \code{TRUE} means that the selected
#' confidence interval excludes \code{threshold}; \code{FALSE} means that it
#' contains \code{threshold}; \code{NA} means that the interval could not be
#' evaluated.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{confints.bootpls}} and \code{\link{signpred}}.
#' @keywords regression models
#' @examples
#' 
#' \donttest{
#' data(Cornell)
#' set.seed(250)
#' modpls <- plsR(Y ~ ., data = Cornell, 3)
#' Cornell.bootYT <- bootpls(modpls, typeboot = "fmodel_np", R = 250, verbose = FALSE)
#' temp.ci <- confints.bootpls(Cornell.bootYT, indices = 2:8, typeBCa = FALSE)
#' confints2signifind(temp.ci)
#' }
#' 
#' @export confints2signifind
confints2signifind <- function(ic_bootobject, typeIC, threshold = 0) {
  if (!(is.matrix(ic_bootobject) || is.data.frame(ic_bootobject))) {
    stop("'ic_bootobject' must be a matrix or data frame returned by confints.bootpls().")
  }
  
  ic_bootobject <- as.matrix(ic_bootobject)
  if (missing(typeIC)) {
    typeIC <- if (isTRUE(attr(ic_bootobject, "typeBCa")) || ncol(ic_bootobject) >= 8L) {
      "BCa"
    } else {
      "Percentile"
    }
  }
  typeIC <- match.arg(typeIC, c("Normal", "Basic", "Percentile", "BCa"))
  ci_columns <- confints_bootpls_columns(ic_bootobject, typeIC)
  lower <- ic_bootobject[, ci_columns[1L], drop = TRUE]
  upper <- ic_bootobject[, ci_columns[2L], drop = TRUE]
  
  signifind <- (lower < threshold & upper < threshold) |
    (lower > threshold & upper > threshold)
  signifind[is.na(lower) | is.na(upper)] <- NA
  names(signifind) <- rownames(ic_bootobject)
  signifind
}

confints_bootpls_columns <- function(ic_bootobject, typeIC) {
  named_columns <- list(
    Normal = c("Normal.Lower", "Normal.Upper"),
    Basic = c("Basic.Lower", "Basic.Upper"),
    Percentile = c("Percentile.Lower", "Percentile.Upper"),
    BCa = c("BCa.Lower", "BCa.Upper")
  )
  positional_columns <- list(
    Normal = c(1L, 2L),
    Basic = c(3L, 4L),
    Percentile = c(5L, 6L),
    BCa = c(7L, 8L)
  )
  
  if (identical(typeIC, "BCa") && !isTRUE(attr(ic_bootobject, "typeBCa")) && ncol(ic_bootobject) < 8L) {
    stop("BCa intervals were not computed, hence cannot be used to detect significance.")
  }
  
  cn <- colnames(ic_bootobject)
  expected_names <- named_columns[[typeIC]]
  if (!is.null(cn) && all(expected_names %in% cn)) {
    return(match(expected_names, cn))
  }
  
  expected_positions <- positional_columns[[typeIC]]
  if (max(expected_positions) > ncol(ic_bootobject)) {
    stop("The supplied confidence-interval matrix does not contain the requested interval type.")
  }
  expected_positions
}
