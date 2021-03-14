#' Boxplot bootstrap distributions
#' 
#' Boxplots for bootstrap distributions.
#' 
#' 
#' @param bootobject a object of class \code{"boot"}
#' @param indices vector of indices of the variables to plot. Defaults to
#' \code{NULL}: all the predictors will be used.
#' @param prednames do the original names of the predictors shall be plotted ?
#' Defaults to \code{TRUE}: the names are plotted.
#' @param articlestyle do the extra blank zones of the margin shall be removed
#' from the plot ? Defaults to \code{TRUE}: the margins are removed.
#' @param xaxisticks do ticks for the x axis shall be plotted ? Defaults to
#' \code{TRUE}: the ticks are plotted.
#' @param ranget0 does the vertival range of the plot shall be computed to
#' include the initial estimates of the coefficients ? Defaults to
#' \code{FALSE}: the vertical range is calculated only using the bootstrapped
#' values of the statistics. Especially using for permutation bootstrap.
#' @param las numeric in 0,1,2,3; the style of axis labels. 0: always parallel
#' to the axis [default], 1: always horizontal, 2: always perpendicular to the
#' axis, 3: always vertical.
#' @param mar A numerical vector of the form \code{c(bottom, left, top, right)}
#' which gives the number of lines of margin to be specified on the four sides
#' of the plot. The default is \code{c(5, 4, 4, 2) + 0.1.}
#' @param mgp The margin line (in mex units) for the axis title, axis labels
#' and axis line. Note that \code{mgp[1]} affects title whereas \code{mgp[2:3]}
#' affect axis. The default is \code{c(3, 1, 0)}.
#' @param \dots further options to pass to the
#' \code{\link[graphics:boxplot]{boxplot}} function.
#' @return \code{NULL}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{bootpls}}
#' @keywords regression models
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' 
#' # Lazraq-Cleroux PLS ordinary bootstrap
#' set.seed(250)
#' modpls <- plsR(yCornell,XCornell,3)
#' Cornell.bootYX <- bootpls(modpls, R=250)
#' 
#' # Graph similar to the one of Bastien et al. in CSDA 2005
#' boxplots.bootpls(Cornell.bootYX,indices=2:8)
#' 
#' \donttest{
#' data(aze_compl)
#' modplsglm<-plsRglm(y~.,data=aze_compl,3,modele="pls-glm-logistic")
#' aze_compl.boot3 <- bootplsglm(modplsglm, R=250, verbose=FALSE)
#' boxplots.bootpls(aze_compl.boot3)
#' boxplots.bootpls(aze_compl.boot3,las=3,mar=c(5,2,1,1))
#' boxplots.bootpls(aze_compl.boot3,indices=c(2,4,6),prednames=FALSE)
#' }
#' 
#' @export boxplots.bootpls
boxplots.bootpls <- function(bootobject,indices=NULL,prednames=TRUE,articlestyle=TRUE,xaxisticks=TRUE,ranget0=FALSE, las=par("las"), mar, mgp,...){
nr <- length(bootobject$t0)
nboot <- dim(bootobject$t)[1]
if(is.null(indices)){indices <- 1:nr}
plotpos <- (1:nr)[1:length(indices)]
if(articlestyle){
  oldparmar <- par("mar")
  oldparmgp <- par("mgp")
  if(missing(mar)){mar=c(2, 2, 1, 1) + 0.1}
  if(missing(mgp)){mgp=c(2, 1, 0)}
  par(mar = mar); par(mgp = mgp)
}
if(!ranget0){boxplot(as.vector(bootobject$t[,indices])~factor(rep(1:length(indices),rep(nboot,length(indices)))),ylim=c(max(-5,min(as.vector(bootobject$t[,indices]))),min(5,max(as.vector(bootobject$t[,indices])))),xaxt="n",...)} else
{boxplot(as.vector(bootobject$t[,indices])~factor(rep(1:length(indices),rep(nboot,length(indices)))),ylim=c(max(min(min(bootobject$t0[indices]),min(as.vector(bootobject$t[,indices])))), min(max(max(bootobject$t0[indices]),max(as.vector(bootobject$t[, indices]))))),xaxt="n",...)}
#if(prednames){axis(1, at = plotpos+.225, labels = rownames(bootobject$t0)[indices])} else {axis(1, at = plotpos+.225, labels = paste("x",(1:nr)[indices],sep=""))} 
if(xaxisticks)
  {
  if(prednames)
    {
    axis(1, at = plotpos, labels = rownames(bootobject$t0)[indices], las=las)
    }
    else
    {
    axis(1, at = plotpos, labels = paste("x",(1:nr)[indices],sep=""), las=las)
    }
  }
  else
  {
  if(prednames)
    {
    axis(1, at = plotpos, labels = rownames(bootobject$t0)[indices],lwd.ticks=0, las=las)
    }
    else
    {
    axis(1, at = plotpos, labels = paste("x",(1:nr)[indices],sep=""),lwd.ticks=0, las=las)
    }
  } 
abline(h=0,lty=2,col="blue",lwd=2)
points(plotpos,bootobject$t0[indices],col="red",pch=19)
if(articlestyle){
  par(mar=oldparmar)
  par(mgp=oldparmgp)
}
}
