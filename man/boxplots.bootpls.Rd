\name{boxplots.bootpls}
\alias{boxplots.bootpls}
\title{Boxplot bootstrap distributions}
\description{
Boxplots for bootstrap distributions.
}
\usage{
boxplots.bootpls(bootobject, indices = NULL, prednames = TRUE, 
articlestyle = TRUE, xaxisticks=TRUE, ranget0= FALSE, las = par("las"), 
mar, mgp, ...)
}
\arguments{
  \item{bootobject}{a object of class \code{"boot"}}
  \item{indices}{vector of indices of the variables to plot. Defaults to \code{NULL}: all the predictors will be used.}
  \item{prednames}{do the original names of the predictors shall be plotted ? Defaults to \code{TRUE}: the names are plotted.}
  \item{articlestyle}{do the extra blank zones of the margin shall be removed from the plot ? Defaults to \code{TRUE}: the margins are removed.}
  \item{xaxisticks}{do ticks for the x axis shall be plotted ? Defaults to \code{TRUE}: the ticks are plotted.}
  \item{ranget0}{does the vertival range of the plot shall be computed to include the initial estimates of the coefficients ? Defaults to \code{FALSE}: the vertical range is calculated only using the bootstrapped values of the statistics. Especially using for permutation bootstrap.}
  \item{las}{numeric in {0,1,2,3}; the style of axis labels. 0: always parallel to the axis [default], 1: always horizontal, 2: always perpendicular to the axis, 3: always vertical.}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)} which gives the number of lines of margin to be specified on the four sides of the plot. The default is \code{c(5, 4, 4, 2) + 0.1.}}
  \item{mgp}{The margin line (in mex units) for the axis title, axis labels and axis line. Note that \code{mgp[1]} affects title whereas \code{mgp[2:3]} affect axis. The default is \code{c(3, 1, 0)}.}
  \item{\dots}{further options to pass to the \code{\link[graphics:boxplot]{boxplot}} function.}
}
\value{\code{NULL}}
%\references{Chambers, J. M., Cleveland, W. S., Kleiner, B. and Tukey, P. A. (1983) Graphical Methods for Data Analysis. Wadsworth & Brooks/Cole.}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{bootpls}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS ordinary bootstrap
set.seed(250)
modpls <- plsR(yCornell,XCornell,3)
Cornell.bootYX <- bootpls(modpls, R=250)

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplots.bootpls(Cornell.bootYX,indices=2:8)

\donttest{
data(aze_compl)
modplsglm<-plsRglm(y~.,data=aze_compl,3,modele="pls-glm-logistic")
aze_compl.boot3 <- bootplsglm(modplsglm, R=250)
boxplots.bootpls(aze_compl.boot3)
boxplots.bootpls(aze_compl.boot3,las=3,mar=c(5,2,1,1))
boxplots.bootpls(aze_compl.boot3,indices=c(2,4,6),prednames=FALSE)
}
}
\keyword{regression}
\keyword{models}
