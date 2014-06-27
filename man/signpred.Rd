\name{signpred}
\alias{signpred}
\title{
Graphical assessment of the stability of selected variables
}
\description{
This fonctions plots, for each of the model, the 
}
\usage{
signpred(matbin, pred.lablength = max(sapply(rownames(matbin), nchar)), 
labsize = 1, plotsize = 12)
}
\arguments{
  \item{matbin}{Matrix with 0 or 1 entries. Each row per predictor and a column for every model. 0 means the predictor is not significant in the model and 1 that, on the contrary, it is significant.}
  \item{pred.lablength}{Maximum length of the predictors labels. Defaults to full label length.}
  \item{labsize}{Size of the predictors labels.}
  \item{plotsize}{Global size of the graph.}
}
\details{
This function is based on the \code{\link[bipartite]{visweb}} function from the bipartite package.
}
\value{A plot window.}
\references{Vazquez, P.D., Chacoff, N.,P. and Cagnolo, L. (2009) Evaluating multiple determinants of the structure of plant-animal mutualistic networks. \emph{Ecology}, 90:2039-2046.}
\author{Bernd Gruber with minor modifications from
\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}}
\seealso{See Also \code{\link{visweb}}}
\examples{
signpred(matrix(rbinom(160,1,.2),ncol=8,dimnames=list(as.character(1:20),as.character(1:8))))
}
\keyword{hplot}
