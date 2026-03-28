\name{predict.plsRmultiModel}
\alias{predict.plsRmultiModel}
\title{Predict method for experimental plsRmulti models}
\usage{
\method{predict}{plsRmultiModel}(
  object,
  newdata,
  comps = object$computed_nt,
  type = c("response", "scores"),
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{object}{An object of class \code{"plsRmultiModel"}.}

  \item{newdata}{Optional predictor matrix or data frame. For formula-fitted
  models, a data frame with the predictor variables used at fit time.}

  \item{comps}{Number of extracted components to use.}

  \item{type}{Either \code{"response"} or \code{"scores"}.}

  \item{verbose}{Should informational messages be displayed?}

  \item{...}{Not used.}
}
\value{
If \code{type = "response"}, a matrix of predicted responses. If
\code{type = "scores"}, a matrix of latent score coordinates.
}
\description{
Prediction method for \code{"plsRmultiModel"} objects.
}
\examples{
set.seed(123)
X <- matrix(rnorm(40 * 3), ncol = 3)
Y <- cbind(
  y1 = X[, 1] + rnorm(40, sd = 0.1),
  y2 = X[, 2] - X[, 3] + rnorm(40, sd = 0.1)
)
fit <- plsRmulti(Y, X, nt = 2, verbose = FALSE)

predict(fit, type = "response")
predict(fit, type = "scores")
}
\seealso{
\code{\link{plsRmulti}}
}
\keyword{methods}
\keyword{predict}
