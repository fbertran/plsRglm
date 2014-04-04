\name{bootpls}
\alias{bootpls}
\title{Non-parametric Bootstrap for PLS models}
\description{
Provides a wrapper for the bootstrap function \code{boot} from the \code{boot} R package.\cr Implements non-parametric bootstrap for PLS generalized linear models by case resampling.
}
\usage{
bootpls(object, typeboot="plsmodel", R=250, statistic=coefs.plsR, sim="ordinary", 
stype="i", stabvalue=1e6,...)
}
\arguments{
  \item{object}{ ~~Explain \code{object} here~~ }
  \item{typeboot}{ ~~Explain \code{typeboot} here~~ }
  \item{R}{ ~~Explain \code{R} here~~ }
  \item{statistic}{ ~~Explain \code{statistic} here~~ }
  \item{sim}{ ~~Explain \code{sim} here~~ }
  \item{stype}{ ~~Explain \code{stype} here~~ }
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
\references{Lazraq et al.\cr Bastien et al. \cr Davison et al.}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{\code{\link[boot:boot]{boot}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)
boot::boot.array(Cornell.boot, indices=TRUE)

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(250,7))), 
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.boot,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.boot,indices=2:8))


\donttest{
library(boot)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=7)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=8)

plot(Cornell.boot,index=2)
boot::jack.after.boot(Cornell.boot, index=2, useJ=TRUE, nt=3)
plot(Cornell.boot,index=2,jack=TRUE)

car::data.ellipse(Cornell.boot$t[,2], Cornell.boot$t[,3], cex=.3, levels=c(.5, .95, .99), robust=T)

rm(Cornell.boot)


# PLS balanced bootstrap

set.seed(225)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="balanced", stype="i", R=250)
boot::boot.array(Cornell.boot, indices=TRUE)


# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(250,7))), 
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.boot,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.boot,indices=2:8))


library(boot)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=7)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=8)

plot(Cornell.boot,index=2)
boot::jack.after.boot(Cornell.boot, index=2, useJ=TRUE, nt=3)
plot(Cornell.boot,index=2,jack=TRUE)

rm(Cornell.boot)

# PLS permutation bootstrap

set.seed(500)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="permutation", stype="i", R=1000)
boot::boot.array(Cornell.boot, indices=TRUE)


# Graph of bootstrap distributions
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(1000,7))),
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)




library(boot)
plot(Cornell.boot,index=2)

qqnorm(Cornell.boot$t[,2],ylim=c(-1,1))
abline(h=Cornell.boot$t0[2],lty=2)
(sum(abs(Cornell.boot$t[,2])>=abs(Cornell.boot$t0[2]))+1)/(length(Cornell.boot$t[,2])+1)

rm(Cornell.boot)
}
}\keyword{models}
