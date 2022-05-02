#' Plot bootstrap confidence intervals
#' 
#' This function plots the confidence intervals derived using the function
#' \code{confints.bootpls} from from a \code{bootpls} based object.
#' 
#' 
#' @param ic_bootobject an object created with the \code{confints.bootpls}
#' function.
#' @param indices vector of indices of the variables to plot. Defaults to
#' \code{NULL}: all the predictors will be used.
#' @param legendpos position of the legend as in
#' \code{\link[graphics:legend]{legend}}, defaults to \code{"topleft"}
#' @param prednames do the original names of the predictors shall be plotted ?
#' Defaults to \code{TRUE}: the names are plotted.
#' @param articlestyle do the extra blank zones of the margin shall be removed
#' from the plot ? Defaults to \code{TRUE}: the margins are removed.
#' @param xaxisticks do ticks for the x axis shall be plotted ? Defaults to
#' \code{TRUE}: the ticks are plotted.
#' @param ltyIC lty as in \code{\link[graphics:plot]{plot}}
#' @param colIC col as in \code{\link[graphics:plot]{plot}}
#' @param typeIC type of CI to plot. Defaults to \code{typeIC=c("Normal",
#' "Basic", "Percentile", "BCa")} if BCa intervals limits were computed and to
#' \code{typeIC=c("Normal", "Basic", "Percentile")} otherwise.
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
#' \code{\link[graphics:plot]{plot}} function.
#' @return \code{NULL}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{confints.bootpls}}
#' @keywords regression models
#' @examples
#' 
#' data(Cornell)
#' modpls <- plsR(Y~.,data=Cornell,3)
#' 
#' # Lazraq-Cleroux PLS (Y,X) bootstrap
#' set.seed(250)
#' Cornell.bootYX <- bootpls(modpls, R=250, verbose=FALSE)
#' temp.ci <- confints.bootpls(Cornell.bootYX,2:8)
#' 
#' plots.confints.bootpls(temp.ci)
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
#' main="Bootstrap confidence intervals for the bj")
#' plots.confints.bootpls(temp.ci,indices=1:3,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,c(2,4,6),"bottomright")
#' plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
#' main="Bootstrap confidence intervals for some of the bj")
#' 
#' temp.ci <- confints.bootpls(Cornell.bootYX,typeBCa=FALSE)
#' plots.confints.bootpls(temp.ci)
#' plots.confints.bootpls(temp.ci,2:8)
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' 
#' 
#' # Bastien CSDA 2005 (Y,T) bootstrap
#' Cornell.boot <- bootpls(modpls, typeboot="fmodel_np", R=250, verbose=FALSE)
#' temp.ci <- confints.bootpls(Cornell.boot,2:8)
#' 
#' plots.confints.bootpls(temp.ci)
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
#' main="Bootstrap confidence intervals for the bj")
#' plots.confints.bootpls(temp.ci,indices=1:3,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,c(2,4,6),"bottomright")
#' plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
#' main="Bootstrap confidence intervals for some of the bj")
#' 
#' temp.ci <- confints.bootpls(Cornell.boot,typeBCa=FALSE)
#' plots.confints.bootpls(temp.ci)
#' plots.confints.bootpls(temp.ci,2:8)
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' 
#' 
#' \donttest{
#' data(aze_compl)
#' modplsglm <- plsRglm(y~.,data=aze_compl,3,modele="pls-glm-logistic")
#' 
#' # Lazraq-Cleroux PLS (Y,X) bootstrap
#' # should be run with R=1000 but takes much longer time
#' aze_compl.bootYX3 <- bootplsglm(modplsglm, typeboot="plsmodel", R=250, verbose=FALSE)
#' temp.ci <- confints.bootpls(aze_compl.bootYX3)
#' 
#' plots.confints.bootpls(temp.ci)
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
#' main="Bootstrap confidence intervals for the bj")
#' plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,c(2,4,6),"bottomleft")
#' plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
#' main="Bootstrap confidence intervals for some of the bj")
#' plots.confints.bootpls(temp.ci,indices=1:34,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE,ltyIC=1,colIC=c(1,2))
#'  
#' temp.ci <- confints.bootpls(aze_compl.bootYX3,1:34,typeBCa=FALSE)
#' plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
#' 
#' 
#' # Bastien CSDA 2005 (Y,T) Bootstrap
#' # much faster
#' aze_compl.bootYT3 <- bootplsglm(modplsglm, R=1000, verbose=FALSE)
#' temp.ci <- confints.bootpls(aze_compl.bootYT3)
#' 
#' plots.confints.bootpls(temp.ci)
#' plots.confints.bootpls(temp.ci,typeIC="Normal")
#' plots.confints.bootpls(temp.ci,typeIC=c("Normal","Basic"))
#' plots.confints.bootpls(temp.ci,typeIC="BCa",legendpos="bottomleft")
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,
#' main="Bootstrap confidence intervals for the bj")
#' plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
#' plots.confints.bootpls(temp.ci,c(2,4,6),"bottomleft")
#' plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,
#' main="Bootstrap confidence intervals for some of the bj")
#' plots.confints.bootpls(temp.ci,prednames=FALSE,ltyIC=c(2,1),colIC=c(1,2))
#'  
#' temp.ci <- confints.bootpls(aze_compl.bootYT3,1:33,typeBCa=FALSE)
#' plots.confints.bootpls(temp.ci,prednames=FALSE)
#' }
#' 
#' @export plots.confints.bootpls
plots.confints.bootpls = function (ic_bootobject, indices = NULL, legendpos = "topleft", 
    prednames = TRUE, articlestyle = TRUE, xaxisticks=TRUE, ltyIC=c(2, 4, 5, 1), colIC=c("darkgreen", "blue", "red", "black"), typeIC, las=par("las"), mar, mgp, ...) 
{  
    if(missing(typeIC)){
      if(attr(ic_bootobject, "typeBCa")){
        typeIC <- c("Normal", "Basic", "Percentile", "BCa")
      }
      else {
        typeIC <- c("Normal", "Basic", "Percentile")
      }
    }
    if((!attr(ic_bootobject, "typeBCa"))&("BCa" %in% typeIC)){stop("BCa intervals were not computed, hence cannot be plotted.")}
    if(length(ltyIC)<length(typeIC)){ltyIC <- rep_len(ltyIC,length(typeIC))}
    if(length(colIC)<length(typeIC)){colIC <- rep_len(colIC,length(typeIC))}
    nr <- nrow(ic_bootobject)
    if (is.null(indices)) {
        indices <- 1:nr
    }
    plotpos <- (1:nr)[1:length(indices)]
    if (articlestyle) {
      oldparmar <- par("mar")
      oldparmgp <- par("mgp")
      if(missing(mar)){mar=c(2, 2, 1, 1) + 0.1}
        if(missing(mgp)){mgp=c(2, 1, 0)}
        par(mar = mar); par(mgp = mgp)
    }
    plot(c(1, 1), xlab = "", ylab = "", type = "n", xlim = c(1, 
        length(indices) + 0.5), ylim = c(min(ic_bootobject[indices,]), 
        max(ic_bootobject[indices,])), xaxt = "n", ...)
        legendtxt <- NULL
        indictypeIC <- rep(FALSE,4)
        nbIC <- 0
        if ("Normal" %in% typeIC){
        indictypeIC[1] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 1], plotpos, ic_bootobject[indices, 
            2], lend = "butt", lwd = 2, lty = ltyIC[1], col = colIC[1], 
            code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"Normal")
        nbIC <- nbIC+1
            }
        if ("Basic" %in% typeIC){        
        indictypeIC[2] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 3], plotpos + 
            nbIC*0.15, ic_bootobject[indices, 4], lend = "butt", lwd = 2, 
            lty = ltyIC[2], col = colIC[2], code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"Basic")
        nbIC <- nbIC+1
            }
        if ("Percentile" %in% typeIC){        
        indictypeIC[3] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 5], plotpos + 
            nbIC*0.15, ic_bootobject[indices, 6], lend = "butt", lwd = 2, 
            lty = ltyIC[3], col = colIC[3], code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"Percentile")
        nbIC <- nbIC+1
        }
        if (("BCa" %in% typeIC)&(attr(ic_bootobject, "typeBCa"))){
        indictypeIC[4] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 7], plotpos + 
            nbIC*0.15, ic_bootobject[indices, 8], lend = "butt", lwd = 2, 
            lty = ltyIC[4], col = colIC[4], code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"BCa")     
        nbIC <- nbIC+1
        }
        if (prednames) {
            if(xaxisticks){
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = rownames(ic_bootobject)[indices], las=las)
            }
            else
            {
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = rownames(ic_bootobject)[indices],lwd.ticks=0, las=las)
            }
        }
        else {
            if(xaxisticks){
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = paste("x", (1:nr)[indices], sep = ""), las=las)
            }
            else
            {
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = paste("x", (1:nr)[indices], sep = ""),lwd.ticks=0, las=las)
            }
        }
        abline(h = 0, lty = 3, lwd = 2)
        legend(legendpos, legend = legendtxt, lty = ltyIC[indictypeIC], col = colIC[indictypeIC], lwd = 2)
    if (articlestyle) {
      par(mar=oldparmar)
      par(mgp=oldparmgp)
    }
}
