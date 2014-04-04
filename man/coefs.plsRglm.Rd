\name{coefs.plsRglm}
\alias{coefs.plsRglm}
\title{Coefficients for bootstrap computations}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
coefs.plsRglm(dataset, ind, nt, modele, family=NULL, maxcoefvalues,ifbootfail)
}
\arguments{
  \item{dataset}{ ~~Explain \code{dataset} here~~ }
  \item{ind}{ ~~Explain \code{ind} here~~ }
  \item{nt}{ ~~Explain \code{nt} here~~ }
  \item{modele}{ ~~Explain \code{modele} here~~ }
  \item{family}{ ~~Explain \code{family} here~~ }
  \item{maxcoefvalues}{ ~~Explain \code{maxcoefvalues} here~~ }
  \item{ifbootfail}{ ~~Explain \code{ifbootfail} here~~ }
}
\details{
  ~~ More details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
}
\keyword{models}
