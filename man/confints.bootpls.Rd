\name{confints.bootpls}
\alias{confints.bootpls}
\title{Bootstrap confidence intervals}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
confints.bootpls(bootobject, indices = NULL, typeBCa=TRUE)
}
\arguments{
  \item{bootobject}{ ~~Explain \code{bootobject} here~~ }
  \item{indices}{ ~~Explain \code{indices} here~~ }
  \item{typeBCa}{ ~~Explain \code{typeBCa} here~~ }
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
\seealso{\code{\link{bootpls}}}
\examples{
\donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)

(temp.ci <- confints.bootpls(Cornell.boot,2:8))
plots.confints.bootpls(temp.ci)
(temp.ci <- confints.bootpls(Cornell.boot,2:8,typeBCa=FALSE))
plots.confints.bootpls(temp.ci)
(temp.ci <- confints.bootpls(Cornell.boot,c(2,4,6)))
plots.confints.bootpls(temp.ci)

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

#should be run with R=1000 but takes much longer time
aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,
modele="pls-glm-logistic"), sim="ordinary", stype="i", R=250)
(temp.ci <- confints.bootpls(aze_compl.boot3))
(temp.ci <- confints.bootpls(aze_compl.boot3,1:34,typeBCa=FALSE))
(temp.ci <- confints.bootpls(aze_compl.boot3,c(2,4,6)))
}
}\keyword{regression}
\keyword{models}
