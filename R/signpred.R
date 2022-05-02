#' Graphical assessment of the stability of selected variables
#' 
#' This fonctions plots, for each of the model, the
#' 
#' This function is based on the \code{\link[bipartite]{visweb}} function from
#' the bipartite package.
#' 
#' @param matbin Matrix with 0 or 1 entries. Each row per predictor and a
#' column for every model. 0 means the predictor is not significant in the
#' model and 1 that, on the contrary, it is significant.
#' @param pred.lablength Maximum length of the predictors labels. Defaults to
#' full label length.
#' @param labsize Size of the predictors labels.
#' @param plotsize Global size of the graph.
#' @return A plot window.
#' @author Bernd Gruber with minor modifications from
#' Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See Also \code{\link{visweb}}
#' @references Vazquez, P.D., Chacoff, N.,P. and Cagnolo, L. (2009) Evaluating
#' multiple determinants of the structure of plant-animal mutualistic networks.
#' \emph{Ecology}, 90:2039-2046.
#' @keywords hplot
#' @examples
#' 
#' signpred(matrix(rbinom(160,1,.2),ncol=8,dimnames=list(as.character(1:20),as.character(1:8))))
#' 
#' @export signpred
signpred <- function(matbin,pred.lablength=max(sapply(rownames(matbin),nchar)),labsize=1,plotsize = 12){
if(!any(matbin)){matbin<-!matbin;colbin="white"} else{colbin="grey25"}
if(is.null(rownames(matbin))){rownames(matbin) <- paste("x",1:nrow(matbin),sep="")}
lll=max(sapply(rownames(matbin),nchar))
text = "no"
ncol <- ncol(matbin)
nrow <- nrow(matbin)
plotsize = plotsize/2.54
mcol = max(matbin)
if (ncol > nrow) {
        wx <- plotsize
        wy <- (plotsize)/ncol * nrow
}
else {
        wy <- plotsize
        wx <- (plotsize)/nrow * ncol
}
m.colsize = max(strwidth(colnames(matbin), units = "inches"))
m.rowsize = max(strwidth(rownames(matbin), units = "inches"))
cellsize = wx/ncol
if (substr(text, 1, 1) == "i") 
        s <- as.character(max(matbin))
else s = "A"
lettersize = strwidth(s, units = "inches")
clratio = cellsize/lettersize
mm <- max(m.colsize, m.rowsize)
op <- par(las=3,mar = c(2, 2, 1, 1) + 0.1, mgp = c(2, 1, 0))
on.exit(par(op))
visweb(t(matbin),type="None",labsize=labsize,square="b",box.col=colbin,prednames=FALSE,clear=FALSE)
text(.5+0:(length(rownames(matbin))-1),-lll+1,rownames(matbin),cex=0.4 * labsize *clratio,srt=90)
}
