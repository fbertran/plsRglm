\name{tilt.bootplsglm}
\alias{tilt.bootplsglm}
\title{Tilted bootstrap for PLS models}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
tilt.bootplsglm(object, typeboot="fmodel_np", statistic=coefs.plsRglm, 
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
\examples{\donttest{

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

dataset <- cbind(y=yaze_compl,Xaze_compl)

library(boot)
# Lazraq-Cl?roux PLS bootstrap Classic

tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), alpha=c(0.025, 0.975), 
sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic", family=NULL)
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic", family=NULL), statistic=coefs.plsR, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-family", family=binomial), statistic=coefs.plsR, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)

boxplots.bootpls(aze_compl.tilt.boot,1:2)


# PLS bootstrap balanced

# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 250, 250), 
alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)


data(aze)
Xaze<-aze[,2:34]
yaze<-aze$y

dataset <- cbind(y=yaze,Xaze)

library(boot)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze.tilt.boot <- tilt.bootplsglm(plsRglm(yaze,Xaze,3, modele="pls-glm-logistic"), 
statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), 
sim="ordinary", stype="i", index=1)
}
}\keyword{models}
