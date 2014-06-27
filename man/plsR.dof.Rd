\name{plsR.dof}
\alias{plsR.dof}
\title{
Computation of the Degrees of Freedom
}
\description{
This function computes the Degrees of Freedom using the Krylov representation of PLS and other quantities that are used to get information criteria values. For the time present, it only works with complete datasets.
}
\usage{
plsR.dof(modplsR, naive = FALSE)
}
\arguments{
  \item{modplsR}{A plsR model i.e. an object returned by one of the functions \code{plsR}, \code{plsRmodel.default}, \code{plsRmodel.formula}, \code{PLS_lm} or \code{PLS_lm_formula}.}
  \item{naive}{A boolean.}
}
\details{
If \code{naive=FALSE} returns values for estimated degrees of freedom and error dispersion. If \code{naive=TRUE} returns returns values for naive degrees of freedom and error dispersion.
The original code from Nicole Kraemer and Mikio L. Braun was unable to handle models with only one component.
}
\value{
  \item{DoF}{Degrees of Freedom}
  \item{sigmahat}{Estimates of dispersion}
  \item{Yhat}{Predicted values}
  \item{yhat}{Square Euclidean norms of the predicted values}
  \item{RSS}{Residual Sums of Squares}
}
\references{
N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of Partial Least Squares Regression. \emph{Journal of the American Statistical Association}, 106(494), 697-705.\cr
N. Kraemer, M. Sugiyama, M.L. Braun. (2009). Lanczos Approximations for the Speedup of Kernel Partial Least Squares Regression, \emph{Proceedings of the Twelfth International Conference on Artificial Intelligence and Statistics (AISTATS)}, 272-279.
}
\author{Nicole Kraemer, Mikio L. Braun with improvements from
\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{aic.dof}} and \code{\link{infcrit.dof}} for computing information criteria directly from a previously fitted plsR model.}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
plsR.dof(modpls) 
plsR.dof(modpls,naive=TRUE) 
}
\keyword{models}
\keyword{regression}
\keyword{utilities}
