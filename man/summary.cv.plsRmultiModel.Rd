\name{summary.cv.plsRmultiModel}
\alias{summary.cv.plsRmultiModel}
\title{Summary method for experimental multivariate PLS2 CV models}
\usage{
\method{summary}{cv.plsRmultiModel}(object, verbose = TRUE, ...)
}
\arguments{
\item{object}{An object of class \code{"cv.plsRmultiModel"}.}

\item{verbose}{Should progress information be displayed?}

\item{\dots}{Further arguments passed to methods.}
}
\value{
A list of per-partition summary matrices with the same aggregate columns used
by \code{\link{summary.cv.plsRmodel}} for \code{Q2}, \code{PRESS}, and
\code{RSS}, plus response-specific \code{PRESS}, \code{RSS}, \code{Q2}, and
\code{R2} columns.
}
\description{
Summarizes repeated k-fold cross-validation results from
\code{\link{cv.plsRmulti}}.
}
\details{
The returned object inherits from \code{"summary.cv.plsRmodel"} so that
\code{\link{cvtable}} and the existing plot method can be reused for the
aggregated multivariate criteria.
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
\code{\link{cv.plsRmulti}}, \code{\link{cvtable}}
}
\keyword{methods}
\keyword{print}
