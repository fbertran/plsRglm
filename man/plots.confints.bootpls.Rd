\name{plots.confints.bootpls}
\alias{plots.confints.bootpls}
\title{Plot bootstrap confidence intervals}
\description{
This function plots the confidence intervals derived using the function \code{confints.bootpls} from from a \code{bootpls} based object.
}
\usage{
plots.confints.bootpls(ic_bootobject, indices = NULL, legendpos = "topleft", 
prednames = TRUE, articlestyle = TRUE, xaxisticks=TRUE, ltyIC=c(2, 4, 5, 1), 
colIC=c("darkgreen", "blue", "red", "black"), typeIC, las = par("las"), 
mar, mgp, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ic_bootobject}{an object created with the \code{confints.bootpls} function.}
  \item{indices}{vector of indices of the variables to plot. Defaults to \code{NULL}: all the predictors will be used.}
  \item{legendpos}{position of the legend as in \code{\link[graphics:legend]{legend}}, defaults to \code{"topleft"}}.
  \item{prednames}{do the original names of the predictors shall be plotted ? Defaults to \code{TRUE}: the names are plotted.}
  \item{articlestyle}{do the extra blank zones of the margin shall be removed from the plot ? Defaults to \code{TRUE}: the margins are removed.}
  \item{xaxisticks}{do ticks for the x axis shall be plotted ? Defaults to \code{TRUE}: the ticks are plotted.}
  \item{ltyIC}{lty as in \code{\link[graphics:plot]{plot}}}
  \item{colIC}{col as in \code{\link[graphics:plot]{plot}}}
  \item{typeIC}{type of CI to plot. Defaults to \code{typeIC=c("Normal", "Basic", "Percentile", "BCa")} if BCa intervals limits were computed and to \code{typeIC=c("Normal", "Basic", "Percentile")} otherwise.}
  \item{las}{numeric in {0,1,2,3}; the style of axis labels. 0: always parallel to the axis [default], 1: always horizontal, 2: always perpendicular to the axis, 3: always vertical.}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)} which gives the number of lines of margin to be specified on the four sides of the plot. The default is \code{c(5, 4, 4, 2) + 0.1.}}
  \item{mgp}{The margin line (in mex units) for the axis title, axis labels and axis line. Note that \code{mgp[1]} affects title whereas \code{mgp[2:3]} affect axis. The default is \code{c(3, 1, 0)}.}
  \item{\dots}{further options to pass to the \code{\link[graphics:plot]{plot}} function.}
}
\value{\code{NULL}}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{confints.bootpls}}}
\examples{
data(Cornell)
modpls <- plsR(Y~.,data=Cornell,3)

# Lazraq-Cleroux PLS (Y,X) bootstrap
set.seed(250)
Cornell.bootYX <- bootpls(modpls, R=250)
temp.ci <- confints.bootpls(Cornell.bootYX,2:8)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:3,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomright")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
main="Bootstrap confidence intervals for some of the bj")

temp.ci <- confints.bootpls(Cornell.bootYX,typeBCa=FALSE)
plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,2:8)
plots.confints.bootpls(temp.ci,prednames=FALSE)


# Bastien CSDA 2005 (Y,T) bootstrap
Cornell.boot <- bootpls(modpls, typeboot="fmodel_np", R=250)
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
plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,2:8)
plots.confints.bootpls(temp.ci,prednames=FALSE)


\donttest{
data(aze_compl)
modplsglm <- plsRglm(y~.,data=aze_compl,3,modele="pls-glm-logistic")

# Lazraq-Cleroux PLS (Y,X) bootstrap
# should be run with R=1000 but takes much longer time
aze_compl.bootYX3 <- bootplsglm(modplsglm, typeboot="plsmodel", R=250)
temp.ci <- confints.bootpls(aze_compl.bootYX3)

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
 
temp.ci <- confints.bootpls(aze_compl.bootYX3,1:34,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)


# Bastien CSDA 2005 (Y,T) Bootstrap
# much faster
aze_compl.bootYT3 <- bootplsglm(modplsglm, R=1000)
temp.ci <- confints.bootpls(aze_compl.bootYT3)

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
plots.confints.bootpls(temp.ci,prednames=FALSE,ltyIC=c(2,1),colIC=c(1,2))
 
temp.ci <- confints.bootpls(aze_compl.bootYT3,1:33,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE)
}
}\keyword{regression}
\keyword{models}

