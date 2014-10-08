\name{bootpls}
\alias{bootpls}
\title{Non-parametric Bootstrap for PLS models}
\description{
Provides a wrapper for the bootstrap function \code{boot} from the \code{boot} R package.\cr Implements non-parametric bootstraps for PLS Regression models by either (Y,X) or (Y,T) resampling.
}
\usage{
bootpls(object, typeboot="plsmodel", R=250, statistic=coefs.plsR, sim="ordinary", 
stype="i", stabvalue=1e6, verbose=TRUE,...)
}
\arguments{
  \item{object}{An object of class \code{plsRmodel} to bootstrap}
  \item{typeboot}{The type of bootstrap. Either (Y,X) boostrap (\code{typeboot="plsmodel"}) or (Y,T) bootstrap (\code{typeboot="fmodel_np"}). Defaults to (Y,X) resampling.}
  \item{R}{The number of bootstrap replicates. Usually this will be a single positive integer. For importance resampling, some resamples may use one set of weights and others use a different set of weights. In this case \code{R} would be a vector of integers where each component gives the number of resamples from each of the rows of weights.}
  \item{statistic}{A function which when applied to data returns a vector containing the statistic(s) of interest. \code{statistic} must take at least two arguments. The first argument passed will always be the original data. The second will be a vector of indices, frequencies or weights which define the bootstrap sample. Further, if predictions are required, then a third argument is required which would be a vector of the random indices used to generate the bootstrap predictions. Any further arguments can be passed to statistic through the \code{...} argument.}
  \item{sim}{A character string indicating the type of simulation required. Possible values are  \code{"ordinary"} (the default), \code{"balanced"}, \code{"permutation"}, or \code{"antithetic"}.}
  \item{stype}{A character string indicating what the second argument of \code{statistic} represents. Possible values of stype are \code{"i"} (indices - the default), \code{"f"} (frequencies), or \code{"w"} (weights).}
  \item{stabvalue}{A value to hard threshold bootstrap estimates computed from atypical resamplings. Especially useful for Generalized Linear Models.}
  \item{verbose}{should info messages be displayed ?}
  \item{\dots}{Other named arguments for \code{statistic} which are passed unchanged each time it is called. Any such arguments to \code{statistic} should follow the arguments which \code{statistic} is required to have for the simulation. Beware of partial matching to arguments of \code{boot} listed above.}
}
\details{
More details on bootstrap techniques are available in the help of the \code{\link[boot:boot]{boot}} function.
}
\value{
  An object of class \code{"boot"}. See the Value part of the help of the function \code{\link[boot:boot]{boot}}.
}
\references{
A. Lazraq, R. Cleroux, and J.-P. Gauchi. (2003). Selecting both latent and explanatory variables in the PLS1 regression model. \emph{Chemometrics and Intelligent Laboratory Systems}, 66(2):117-126.\cr 
P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS generalised linear regression. \emph{Computational Statistics & Data Analysis}, 48(1):17-46.\cr 
A. C. Davison and D. V. Hinkley. (1997). \emph{Bootstrap Methods and Their Applications}. Cambridge University Press, Cambridge.
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[boot:boot]{boot}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS ordinary bootstrap
set.seed(250)
modpls <- plsR(yCornell,XCornell,3)

#(Y,X) resampling
Cornell.bootYX <- bootpls(modpls, R=250)

#(Y,T) resampling
Cornell.bootYT <- bootpls(modpls, typeboot="fmodel_np", R=250)

# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.bootYX,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.bootYX,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.bootYX,indices=2:8))
# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.bootYX$t[,-1])~factor(rep(1:7,rep(250,7))), 
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.bootYX$t0[-1],col="red",pch=19)


\donttest{
library(boot)
boot.ci(Cornell.bootYX, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
plot(Cornell.bootYX,index=2)
jack.after.boot(Cornell.bootYX, index=2, useJ=TRUE, nt=3)
plot(Cornell.bootYX,index=2,jack=TRUE)

car::dataEllipse(Cornell.bootYX$t[,2], Cornell.bootYX$t[,3], cex=.3, 
levels=c(.5, .95, .99), robust=T)
rm(Cornell.bootYX)


# PLS balanced bootstrap

set.seed(225)
Cornell.bootYX <- bootpls(modpls, sim="balanced", R=250)
boot.array(Cornell.bootYX, indices=TRUE)

# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.bootYX,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.bootYX,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.bootYX,indices=2:8))
# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.bootYX$t[,-1])~factor(rep(1:7,rep(250,7))), 
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.bootYX$t0[-1],col="red",pch=19)


library(boot)
boot.ci(Cornell.bootYX, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
plot(Cornell.bootYX,index=2)
jack.after.boot(Cornell.bootYX, index=2, useJ=TRUE, nt=3)
plot(Cornell.bootYX,index=2,jack=TRUE)

rm(Cornell.bootYX)

# PLS permutation bootstrap

set.seed(500)
Cornell.bootYX <- bootpls(modpls, sim="permutation", R=1000)
boot.array(Cornell.bootYX, indices=TRUE)


# Graph of bootstrap distributions
boxplot(as.vector(Cornell.bootYX$t[,-1])~factor(rep(1:7,rep(1000,7))),
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.bootYX$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.bootYX,indices=2:8)


library(boot)
plot(Cornell.bootYX,index=2)

qqnorm(Cornell.bootYX$t[,2],ylim=c(-1,1))
abline(h=Cornell.bootYX$t0[2],lty=2)
(sum(abs(Cornell.bootYX$t[,2])>=abs(Cornell.bootYX$t0[2]))+1)/(length(Cornell.bootYX$t[,2])+1)

rm(Cornell.bootYX)
}
}\keyword{models}
