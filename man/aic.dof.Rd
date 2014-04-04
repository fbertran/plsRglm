\name{aic.dof}
\alias{aic.dof}
\alias{bic.dof}
\alias{gmdl.dof}
\title{
Akaike and Bayesian Information Criteria and Generalized minimum description length 
}
\description{
This function computes the Akaike and Bayesian Information Criteria and the Generalized minimum description length.
}
\usage{
aic.dof(RSS, n, DoF, sigmahat)
bic.dof(RSS, n, DoF, sigmahat)
gmdl.dof(sigmahat, n, DoF, yhat)
}
\arguments{
  \item{RSS}{vector of residual sum of squares.}
  \item{n}{number of observations.}
  \item{DoF}{vector of Degrees of Freedom. The length of \code{DoF} is the same as the length of \code{RSS}.}
  \item{sigmahat}{Estimated model error. The length of \code{sigmahat} is the same as the length of \code{RSS}.}
  \item{yhat}{vector of squared norm of Yhat. The length of \code{yhat} is the same as the length of \code{sigmahat}.}
}
\details{
The gmdl criterion is defined as 
\deqn{gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)}{gmdl=n/2*log(S)+DoF/2*log(F)+1/2*log(n)}
with 
\deqn{S=\hat\sigma^2}{\sigma hat^2}
}
\value{
  \item{vector}{numerical values of the requested AIC, BIC or GMDL.}
}
\references{
M. Hansen, B. Yu, Model Selection and Minimum Descripion Length Principle, Journal of the American Statistical Association, 96 (2001) 746-774.\cr
N. Kraemer, M. Sugiyama, The Degrees of Freedom of Partial Least Squares Regression. Preprint (2010). \url{http://arxiv.org/abs/1002.4112}.\cr
N. Kraemer, M.L. Braun, Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection,
\emph{Proceedings of the 24th International Conference on Machine Learning}, Omni Press, (2007) 441-448. 
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
dof.object <- plsR.dof(modpls)
aic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
bic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
gmdl.dof(dof.object$sigmahat,modpls$nr,dof.object$DoF,dof.object$yhat)
naive.object <- plsR.dof(modpls,naive=TRUE)
aic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
bic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
gmdl.dof(naive.object$sigmahat,modpls$nr,naive.object$DoF,naive.object$yhat)
}
\keyword{models}
\keyword{regression}
\keyword{utilities}
