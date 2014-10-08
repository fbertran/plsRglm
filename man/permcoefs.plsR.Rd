\name{permcoefs.plsR}
\alias{permcoefs.plsR}
\title{Coefficients for permutation bootstrap computations of PLSR models}
\description{
A function passed to \code{boot} to perform bootstrap.
}
\usage{
permcoefs.plsR(dataset,ind,nt,modele, maxcoefvalues,ifbootfail,verbose)
}
\arguments{
  \item{dataset}{dataset to resample}
  \item{ind}{indices for resampling}
  \item{nt}{number of components to use}
  \item{modele}{type of modele to use, see \link{plsR}}
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
\seealso{See also \code{\link{bootpls}}.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS (Y,X) bootstrap
# statistic=permcoefs.plsR is the default for (Y,X) permutation resampling of PLSR models.
set.seed(250)
modpls <- plsR(yCornell,XCornell,1)
Cornell.bootYX <- bootpls(modpls, sim="permutation", R=250, statistic=permcoefs.plsR)
}
\keyword{models}
