\name{infcrit.dof}
\alias{infcrit.dof}
\title{
Information criteria
}
\description{
This function computes information criteria for existing plsR model using Degrees of Freedom estimation.
}
\usage{
infcrit.dof(modplsR, naive = FALSE)
}
\arguments{
  \item{modplsR}{A plsR model i.e. an object returned by one of the functions \code{plsR}, \code{plsRmodel.default}, \code{plsRmodel.formula}, \code{PLS_lm} or \code{PLS_lm_formula}.}
  \item{naive}{A boolean.}
}
\details{
If \code{naive=FALSE} returns AIC, BIC and gmdl values for estimated and naive degrees of freedom. If \code{naive=TRUE} returns \code{NULL}.
}
\value{
\item{matrix}{AIC, BIC and gmdl values or \code{NULL}.}
}
\references{
M. Hansen, B. Yu, Model Selection and Minimum Descripion Length Principle, \emph{Journal of the American Statistical Association}, 96 (2001) 746-774.\cr 
N. Kraemer, M. Sugiyama, The Degrees of Freedom of Partial Least Squares Regression. Preprint (2010). \url{http://arxiv.org/abs/1002.4112}.\cr 
N. Kraemer, M.L. Braun, Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection, \emph{Proceedings of the 24th International Conference on Machine Learning}, Omni Press, (2007) 441-448. 
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{plsR.dof}} for degrees of freedom computation and \code{\link{infcrit.dof}} for computing information criteria directly from a previously fitted plsR model.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
infcrit.dof(modpls)
}
\keyword{models}
\keyword{regression}
\keyword{utilities}
