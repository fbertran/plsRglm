\name{cv.plsRmulti}
\alias{cv.plsRmulti}
\alias{cv.plsRmultiModel}
\alias{cv.plsRmultiModel.default}
\alias{cv.plsRmultiModel.formula}
\title{Experimental k-fold cross-validation for multivariate-response PLS2}
\usage{
cv.plsRmulti(object, ...)

\method{cv.plsRmultiModel}{default}(
  object,
  dataX,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  K = 5,
  NK = 1,
  grouplist = NULL,
  random = TRUE,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepfolds = FALSE,
  keepdataY = TRUE,
  keepMclassed = FALSE,
  EstimXNA = FALSE,
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

\method{cv.plsRmultiModel}{formula}(
  object,
  data = NULL,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  K = 5,
  NK = 1,
  grouplist = NULL,
  random = TRUE,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepfolds = FALSE,
  keepdataY = TRUE,
  keepMclassed = FALSE,
  EstimXNA = FALSE,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights = NULL,
  subset = NULL,
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

\item{\dots}{Not used. Extra arguments are rejected in this experimental
release.}

\item{dataX}{Numeric predictor matrix or data frame.}

\item{nt}{Number of components to extract in each fold fit.}

\item{limQ2set}{Threshold used by \code{\link{cvtable}} for the aggregated
\code{Q2} criterion.}

\item{modele}{Only \code{"pls"} is supported.}

\item{family}{Not supported in this experimental release.}

\item{K}{Number of groups for each partition.}

\item{NK}{Number of repeated partitions.}

\item{grouplist}{Optional user-supplied partitions.}

\item{random}{Should the folds be generated randomly?}

\item{scaleX}{Should predictors be scaled?}

\item{scaleY}{Should responses be scaled? Defaults to \code{TRUE}.}

\item{keepcoeffs}{Should standardized coefficient vectors be stored for each
fold fit?}

\item{keepfolds}{Should training indices be stored for each fold fit?}

\item{keepdataY}{Kept for interface compatibility. Observed fold responses are
stored so that summaries can be computed.}

\item{keepMclassed}{Not supported in this experimental release.}

\item{EstimXNA}{Not supported in this experimental release.}

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
An object of class \code{"cv.plsRmultiModel"} with repeated fold predictions,
observed fold responses, optional coefficient vectors and fold indices, and the
reference full-data \code{"plsRmultiModel"} fit used for aggregated summary
metrics.
}
\description{
\code{cv.plsRmulti()} performs repeated k-fold cross-validation for the
experimental complete-case linear \code{\link{plsRmulti}} workflow.
}
\details{
Only the linear multivariate-response PLS2 mode is supported here. Missing
values, weights, sparse extraction options, classification diagnostics, and
GLM families remain out of scope for this experimental API.
}
\examples{
set.seed(123)
X <- matrix(rnorm(60 * 4), ncol = 4)
Y <- cbind(
  y1 = X[, 1] - 0.5 * X[, 2] + rnorm(60, sd = 0.1),
  y2 = 0.3 * X[, 2] + X[, 3] + rnorm(60, sd = 0.1)
)

cv_fit <- cv.plsRmulti(Y, X, nt = 2, K = 3, NK = 1, verbose = FALSE)
summary(cv_fit, verbose = FALSE)
}
\seealso{
\code{\link{plsRmulti}}, \code{\link{summary.cv.plsRmultiModel}},
\code{\link{cvtable}}, \code{\link{bootpls}}
}
\keyword{models}
\keyword{regression}
