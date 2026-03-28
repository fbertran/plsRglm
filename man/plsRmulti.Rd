\name{plsRmulti}
\alias{plsRmulti}
\alias{plsRmultiModel}
\alias{plsRmultiModel.default}
\alias{plsRmultiModel.formula}
\title{Experimental multivariate-response PLS2 models}
\usage{
plsRmulti(object, ...)

\method{plsRmultiModel}{default}(
  object,
  dataX,
  nt = 2,
  limQ2set = 0.0975,
  dataPredictY,
  modele = "pls",
  family = NULL,
  typeVC = "none",
  EstimXNA = FALSE,
  scaleX = TRUE,
  scaleY = NULL,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights,
  sparse = FALSE,
  sparseStop = FALSE,
  naive = FALSE,
  verbose = TRUE,
  ...
)

\method{plsRmultiModel}{formula}(
  object,
  data,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  typeVC = "none",
  EstimXNA = FALSE,
  scaleX = TRUE,
  scaleY = NULL,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights,
  subset,
  contrasts = NULL,
  sparse = FALSE,
  sparseStop = FALSE,
  naive = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{object}{For the default method, a numeric multivariate response matrix
  or data frame with at least two columns. For the formula method, a formula of
  the form \code{cbind(y1, y2, ...) ~ .}.}

  \item{...}{Not used. Extra arguments are rejected in this experimental
  release.}

  \item{dataX}{Numeric predictor matrix or data frame.}

  \item{nt}{Number of components to extract.}

  \item{limQ2set}{Kept for interface compatibility. Not supported by
  \code{plsRmulti}.}

  \item{dataPredictY}{Kept for interface compatibility. Not supported by
  \code{plsRmulti}; fit first and then use \code{\link{predict}}.}

  \item{modele}{Only \code{"pls"} is supported.}

  \item{family}{Not supported in this experimental release.}

  \item{typeVC}{Only \code{"none"} is supported.}

  \item{EstimXNA}{Not supported in this experimental release.}

  \item{scaleX}{Should predictors be scaled?}

  \item{scaleY}{Should responses be scaled? Defaults to \code{TRUE}.}

  \item{pvals.expli}{Not supported in this experimental release.}

  \item{alpha.pvals.expli}{Not supported in this experimental release.}

  \item{MClassed}{Not supported in this experimental release.}

  \item{tol_Xi}{Tolerance used for degeneracy checks during component
  extraction.}

  \item{weights}{Not supported in this experimental release.}

  \item{sparse}{Not supported in this experimental release.}

  \item{sparseStop}{Not supported in this experimental release.}

  \item{naive}{Not supported in this experimental release.}

  \item{verbose}{Should informational messages be displayed?}

  \item{data}{An optional data frame for the formula method.}

  \item{subset}{An optional subset for the formula method.}

  \item{contrasts}{Optional contrasts for the formula method.}
}
\value{
An object of class \code{"plsRmultiModel"} with multivariate analogues of the
linear \code{plsR} outputs, including the extracted scores \code{tt}, X
loadings \code{pp}, response score coefficients \code{CoeffC}, coefficient
matrix \code{Coeffs}, intercept vector \code{CoeffConstante}, scaled response
matrix \code{RepY}, and fitted response matrices \code{YChapeau},
\code{Std.ValsPredictY}, and \code{ValsPredictY}.
}
\description{
\code{plsRmulti()} implements an experimental complete-case linear PLS2 fit for
multivariate numeric responses. It is intentionally separate from
\code{\link{plsR}} so the current PLS1 API remains unchanged.
}
\details{
This experimental release supports complete-case linear PLS2 fitting,
prediction, repeated k-fold cross-validation via \code{\link{cv.plsRmulti}},
and bootstrap resampling via \code{\link{bootpls}}. It still does not support
missing values, weights, sparse extraction, classification diagnostics, or GLM
families.
}
\examples{
set.seed(123)
X <- matrix(rnorm(60 * 4), ncol = 4)
Y <- cbind(
  y1 = X[, 1] - 0.5 * X[, 2] + rnorm(60, sd = 0.1),
  y2 = 0.3 * X[, 2] + X[, 3] + rnorm(60, sd = 0.1)
)

fit <- plsRmulti(Y, X, nt = 2, verbose = FALSE)
fit
head(predict(fit))
}
\seealso{
\code{\link{predict.plsRmultiModel}}, \code{\link{cv.plsRmulti}},
\code{\link{bootpls}}, \code{\link{plsR}}
}
\keyword{models}
\keyword{regression}
