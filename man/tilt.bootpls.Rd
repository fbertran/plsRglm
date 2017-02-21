\name{tilt.bootpls}
\alias{tilt.bootpls}
\title{Tilted bootstrap for PLS models}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
tilt.bootpls(object, typeboot="plsmodel", statistic=coefs.plsR, 
R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", 
stype="i", index=1, stabvalue=1e6,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Explain \code{object} here~~ }
  \item{typeboot}{ ~~Explain \code{typeboot} here~~ }
  \item{statistic}{ ~~Explain \code{statistic} here~~ }
  \item{R}{ ~~Explain \code{R} here~~ }
  \item{alpha}{ ~~Explain \code{alpha} here~~ }
  \item{sim}{ ~~Explain \code{sim} here~~ }
  \item{stype}{ ~~Explain \code{stype} here~~ }
  \item{index}{ ~~Explain \code{index} here~~ }
  \item{stabvalue}{ ~~Explain \code{stabvalue} here~~ }
  \item{\dots}{ ~~Explain \code{\dots} here~~ }
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
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[boot:boot]{tilt.boot}}}
\examples{
\dontrun{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

set.seed(1385)
Cornell.tilt.boot <- tilt.bootpls(plsR(yCornell,XCornell,3), statistic=coefs.plsR, 
R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
Cornell.tilt.boot
str(Cornell.tilt.boot)

boxplots.bootpls(Cornell.tilt.boot,indices=2:7)

rm(Cornell.tilt.boot)
}
}
\keyword{models}
