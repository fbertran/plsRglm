\name{permcoefs.plsRglm}
\alias{permcoefs.plsRglm}
\title{Coefficients for permutation bootstrap computations of PLSGLR models}
\description{
A function passed to \code{boot} to perform bootstrap.
}
\usage{
permcoefs.plsRglm(dataset, ind, nt, modele, family=NULL, maxcoefvalues,ifbootfail,verbose)
}
\arguments{
  \item{dataset}{dataset to resample}
  \item{ind}{indices for resampling}
  \item{nt}{number of components to use}
  \item{modele}{type of modele to use, see \link{plsRglm}}
  \item{family}{glm family to use, see \link{plsRglm}}
  \item{maxcoefvalues}{maximum values allowed for the estimates of the coefficients to discard those coming from singular bootstrap samples}
  \item{ifbootfail}{value to return if the estimation fails on a bootstrap sample}
  \item{verbose}{should info messages be displayed ?}
}
\value{estimates on a bootstrap sample or \code{ifbootfail} value if the bootstrap computation fails.}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{See also \code{\link{bootplsglm}}.}
\examples{
data(Cornell)

# (Y,X) bootstrap of a PLSGLR model
# statistic=coefs.plsRglm is the default for (Y,X) bootstrap of a PLSGLR models.
set.seed(250)
modplsglm <- plsRglm(Y~.,data=Cornell,1,modele="pls-glm-family",family=gaussian)
Cornell.bootYX <- bootplsglm(modplsglm, R=250, typeboot="plsmodel", 
sim="permutation", statistic=permcoefs.plsRglm)
}
\keyword{models}
