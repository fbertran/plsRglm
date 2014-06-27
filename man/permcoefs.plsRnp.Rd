\name{permcoefs.plsRnp}
\alias{permcoefs.plsRnp}
\title{Coefficients computation for permutation bootstrap}
\description{
A function passed to \code{boot} to perform bootstrap.
}
\usage{
permcoefs.plsRnp(dataRepYtt,ind,nt,modele, maxcoefvalues,wwetoile,ifbootfail)
}
\arguments{
  \item{dataRepYtt}{components' coordinates to bootstrap}
  \item{ind}{indices for resampling}
  \item{nt}{number of components to use}
  \item{modele}{type of modele to use, see \link{plsRglm}}
  \item{maxcoefvalues}{maximum values allowed for the estimates of the coefficients to discard those coming from singular bootstrap samples}
  \item{wwetoile}{values of the Wstar matrix in the original fit}  
  \item{ifbootfail}{value to return if the estimation fails on a bootstrap sample}
}
\value{estimates on a bootstrap sample or \code{ifbootfail} value if the bootstrap computation fails.}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{See also \code{\link{bootpls}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS (Y,X) bootstrap
# statistic=coefs.plsR is the default for (Y,X) resampling of PLSR models.
set.seed(250)
modpls <- plsR(yCornell,XCornell,1)
Cornell.bootYT <- bootpls(modpls, R=250, typeboot="fmodel_np", sim="permutation",
statistic=permcoefs.plsRnp)
}
\keyword{models}
