\name{kfolds2CVinfos_lm}
\alias{kfolds2CVinfos_lm}
\title{Extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares models}
\description{
This function extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares models for both formula or classic specifications of the model.
}
\usage{
kfolds2CVinfos_lm(pls_kfolds, MClassed = FALSE,verbose=TRUE)
}
\arguments{
  \item{pls_kfolds}{an object computed using \code{\link{PLS_lm_kfoldcv}}}
  \item{MClassed}{should number of miss classed be computed}
  \item{verbose}{should infos be displayed ?}
}
\details{
The Mclassed option should only set to \code{TRUE} if the response is binary. 
}
\value{
  \item{list}{table of fit statistics for first group partition}
  \item{\dots}{\dots}
  \item{list}{table of fit statistics for last group partition}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{summary}} and \code{\link{cv.plsR}} instead.}
\seealso{\code{\link{kfolds2coeff}}, \code{\link{kfolds2Pressind}}, \code{\link{kfolds2Press}}, \code{\link{kfolds2Mclassedind}} and \code{\link{kfolds2Mclassed}} to extract and transforms results from k-fold cross-validation.}
\examples{\donttest{
data(Cornell)
summary(cv.plsR(Y~.,data=Cornell,nt=10,K=6,verbose=FALSE))


data(pine)
summary(cv.plsR(x11~.,data=pine,nt=10,NK=3,verbose=FALSE),verbose=FALSE)
data(pineNAX21)
summary(cv.plsR(x11~.,data=pineNAX21,nt=10,NK=3,
verbose=FALSE),verbose=FALSE)


data(aze_compl)
summary(cv.plsR(y~.,data=aze_compl,nt=10,K=8,NK=3,
verbose=FALSE),MClassed=TRUE,verbose=FALSE)
}
}
\keyword{models}
\keyword{regression}
