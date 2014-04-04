\name{plots.confints.bootpls}
\alias{plots.confints.bootpls}
\title{Plot bootstrap confidence intervals}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
plots.confints.bootpls(ic_bootobject, indices = NULL, legendpos = "topleft", 
prednames = TRUE, articlestyle = TRUE, xaxisticks=TRUE, ltyIC=c(2, 4, 5, 1), 
colIC=c("darkgreen", "blue", "red", "black"), typeIC, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ic_bootobject}{ ~~Explain \code{ic_bootobject} here~~ }
  \item{indices}{ ~~Explain \code{indices} here~~ }
  \item{legendpos}{ ~~Explain \code{legendpos} here~~ }
  \item{prednames}{ ~~Explain \code{prednames} here~~ }
  \item{articlestyle}{ ~~Explain \code{articlestyle} here~~ }
  \item{xaxisticks}{ ~~Explain \code{articlestyle} here~~ }
  \item{ltyIC}{ ~~Explain \code{ltyIC} here~~ }
  \item{colIC}{ ~~Explain \code{colIC} here~~ }
  \item{typeIC}{ ~~Explain \code{typeIC} here~~ }
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
\seealso{\code{\link{confints.bootpls}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), typeboot="plsmodel", 
sim="ordinary", stype="i", R=250)
temp.ci <- confints.bootpls(Cornell.boot,2:8)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:3,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomright")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
main="Bootstrap confidence intervals for some of the bj")

temp.ci <- confints.bootpls(Cornell.boot,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,prednames=TRUE)
plots.confints.bootpls(temp.ci,prednames=FALSE)


# Bastien CSDA 2005 bootstrap

Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), typeboot="fmodel_np", 
sim="ordinary", stype="i", R=250)
temp.ci <- confints.bootpls(Cornell.boot,2:8)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:3,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomright")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
main="Bootstrap confidence intervals for some of the bj")

temp.ci <- confints.bootpls(Cornell.boot,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,prednames=TRUE)
plots.confints.bootpls(temp.ci,prednames=FALSE)


\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

# Lazraq-Cléroux PLS ordinary bootstrap
#should be run with R=1000 but takes much longer time
aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), 
typeboot="plsmodel", sim="ordinary", stype="i", R=250)
temp.ci <- confints.bootpls(aze_compl.boot3)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomleft")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
main="Bootstrap confidence intervals for some of the bj")
plots.confints.bootpls(temp.ci,indices=1:34,prednames=FALSE)
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE,ltyIC=1,colIC=c(1,2))
 
temp.ci <- confints.bootpls(aze_compl.boot3,1:34,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)




# Bastien CSDA 2005 Bootstrap
#should be run with R=1000 but takes much longer time
aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), 
typeboot="fmodel_np", sim="ordinary", stype="i", R=1000)
temp.ci <- confints.bootpls(aze_compl.boot3)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,typeIC="Normal")
plots.confints.bootpls(temp.ci,typeIC=c("Normal","Basic"))
plots.confints.bootpls(temp.ci,typeIC="BCa",legendpos="bottomleft")
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomleft")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
main="Bootstrap confidence intervals for some of the bj")
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE,ltyIC=c(2,1),colIC=c(1,2))
 
temp.ci <- confints.bootpls(aze_compl.boot3,1:34,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
}
}\keyword{regression}
\keyword{models}

