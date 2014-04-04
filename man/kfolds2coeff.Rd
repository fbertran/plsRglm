\name{kfolds2coeff}
\alias{kfolds2coeff}
\title{Extracts coefficients from k-fold cross validated partial least squares regression models}
\description{
This fonction extracts coefficients from k-fold cross validated partial least squares regression models
}
\usage{
kfolds2coeff(pls_kfolds)
}
\arguments{
  \item{pls_kfolds}{an object that is a k-fold cross validated partial least squares regression models either lm or glm}
}
\details{
This fonctions works for "lm" or "glm" plsR models.
}
\value{
  \item{coef.all}{matrix with the values of the coefficients for each leave one out step or \code{NULL} if another type of cross validation was used.}
}
\references{ ~put references to the literature/web site here ~ }
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Only for \code{NK=1} and leave one out CV}
\seealso{\code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}}, \code{\link{kfolds2Mclassed}} and \code{\link{kfolds2CVinfos_lm}} to extract and transform results from k-fold cross validation.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE)
bbb2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=4,keepcoeffs=TRUE)
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,2])
rm(list=c("XCornell","yCornell","Xpine","ypine","bbb","bbb2"))
}
\keyword{models}
\keyword{regression}
