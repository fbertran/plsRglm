\name{boxplots.bootpls}
\alias{boxplots.bootpls}
\title{Boxplot bootstrap distributions}
\description{
Boxplot plots for bootstrap distributions.
}
\usage{
boxplots.bootpls(bootobject, indices = NULL, prednames = TRUE, 
articlestyle = TRUE, xaxisticks=TRUE, ranget0= FALSE, ...)
}
\arguments{
  \item{bootobject}{ ~~Explain \code{bootobject} here~~ }
  \item{indices}{ ~~Explain \code{indices} here~~ }
  \item{prednames}{ ~~Explain \code{prednames} here~~ }
  \item{articlestyle}{ ~~Explain \code{articlestyle} here~~ }
  \item{xaxisticks}{ ~~Explain \code{articlestyle} here~~ }
  \item{ranget0}{ ~~Explain \code{articlestyle} here~~ }
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
\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 

}
\seealso{\code{\link{bootpls}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplots.bootpls(Cornell.boot,indices=2:8)

\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

#should be run with R=1000 but takes much longer time
aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,
modele="pls-glm-logistic"), sim="ordinary", stype="i", R=250)
boxplots.bootpls(aze_compl.boot3)
boxplots.bootpls(aze_compl.boot3,prednames=FALSE)
boxplots.bootpls(aze_compl.boot3,prednames=FALSE,articlestyle=FALSE,
main="Bootstrap distribution for the bj")
boxplots.bootpls(aze_compl.boot3,indices=1:34,prednames=FALSE)
boxplots.bootpls(aze_compl.boot3,indices=c(2,4,6),prednames=FALSE)
}
}
\keyword{regression}
\keyword{models}
