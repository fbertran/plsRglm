\name{confints.bootpls}
\alias{confints.bootpls}
\title{Bootstrap confidence intervals}
\description{
This function is a wrapper for \code{\link[boot:boot.ci]{boot.ci}} to derive bootstrap-based confidence intervals from a \code{"boot"} object.
}
\usage{
confints.bootpls(bootobject, indices = NULL, typeBCa=TRUE)
}
\arguments{
  \item{bootobject}{an object of class \code{"boot"}}
  \item{indices}{the indices of the predictor for which CIs should be calculated. Defaults to \code{NULL}: all the predictors will be used.}
  \item{typeBCa}{shall BCa bootstrap based CI derived ? Defaults to \code{TRUE}. This is a safety option since sometimes computing BCa bootstrap based CI fails whereas the other types of CI can still be derived.}
}
\value{Matrix with the limits of bootstrap based CI for all (defaults) or only the selected predictors (\code{indices} option). The limits are given in that order: Normal Lower then Upper Limit, Basic Lower then Upper Limit, Percentile Lower then Upper Limit, BCa Lower then Upper Limit.}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{See also \code{\link{bootpls}} and \code{\link{bootplsglm}}.}
\examples{
\donttest{
data(Cornell)

#Lazraq-Cleroux PLS (Y,X) bootstrap
set.seed(250)
modpls <- plsR(Y~.,data=Cornell,3)
Cornell.bootYX <- bootpls(modpls, R=250, verbose=FALSE)
confints.bootpls(Cornell.bootYX,2:8)
confints.bootpls(Cornell.bootYX,2:8,typeBCa=FALSE)
}
}\keyword{regression}
\keyword{models}
