boxplots.bootpls <- function(bootobject,indices=NULL,prednames=TRUE,articlestyle=TRUE,xaxisticks=TRUE,ranget0=FALSE,...){
nr <- length(bootobject$t0)
nboot <- dim(bootobject$t)[1]
if(is.null(indices)){indices <- 1:nr}
plotpos <- (1:nr)[1:length(indices)]
if(articlestyle){oldpar <- par();par(mar = c(2, 2, 1, 1) + 0.1, mgp=c(2,1,0))}
if(!ranget0){boxplot(as.vector(bootobject$t[,indices])~factor(rep(1:length(indices),rep(nboot,length(indices)))),ylim=c(max(-5,min(as.vector(bootobject$t[,indices]))),min(5,max(as.vector(bootobject$t[,indices])))),xaxt="n",...)} else
{boxplot(as.vector(bootobject$t[,indices])~factor(rep(1:length(indices),rep(nboot,length(indices)))),ylim=c(max(-5, min(min(bootobject$t0[indices]),min(as.vector(bootobject$t[,indices])))), min(5, max(max(bootobject$t0[indices]),max(as.vector(bootobject$t[, indices]))))),xaxt="n",...)}
#if(prednames){axis(1, at = plotpos+.225, labels = rownames(bootobject$t0)[indices])} else {axis(1, at = plotpos+.225, labels = paste("x",(1:nr)[indices],sep=""))} 
if(xaxisticks)
  {
  if(prednames)
    {
    axis(1, at = plotpos, labels = rownames(bootobject$t0)[indices])
    }
    else
    {
    axis(1, at = plotpos, labels = paste("x",(1:nr)[indices],sep=""))
    }
  }
  else
  {
  if(prednames)
    {
    axis(1, at = plotpos, labels = rownames(bootobject$t0)[indices],lwd.ticks=0)
    }
    else
    {
    axis(1, at = plotpos, labels = paste("x",(1:nr)[indices],sep=""),lwd.ticks=0)
    }
  } 
abline(h=0,lty=2,col="blue",lwd=2)
points(plotpos,bootobject$t0[indices],col="red",pch=19)
if(articlestyle){par(oldpar)}
}
