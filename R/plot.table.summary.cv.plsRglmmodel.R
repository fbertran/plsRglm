plot.table.summary.cv.plsRglmmodel <- function(x,type=c("CVMC","CVQ2","CVPreChi2"), ...)
{
resCV=x[[type]]
mp<-barplot(resCV,col="lightblue")
text(mp, pmax(resCV/2,0.5), format(resCV/(sum(resCV)),digits = 2,nsmall=2), xpd = TRUE, col = "red")
}


