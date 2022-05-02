#' Akaike and Bayesian Information Criteria and Generalized minimum description
#' length
#' 
#' This function computes the Akaike and Bayesian Information Criteria and the
#' Generalized minimum description length.
#' 
#' The gmdl criterion is defined as
#' \deqn{gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)}{gmdl=n/2*log(S)+DoF/2*log(F)+1/2*log(n)}
#' with \deqn{S=\hat\sigma^2}{\sigma hat^2}
#' 
#' @aliases aic.dof
#' @param RSS vector of residual sum of squares.
#' @param n number of observations.
#' @param DoF vector of Degrees of Freedom. The length of \code{DoF} is the
#' same as the length of \code{RSS}.
#' @param sigmahat Estimated model error. The length of \code{sigmahat} is the
#' same as the length of \code{RSS}.
#' @param yhat vector of squared norm of Yhat. The length of \code{yhat} is the
#' same as the length of \code{sigmahat}.
#' @return \item{vector}{numerical values of the requested AIC, BIC or GMDL.}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{plsR.dof}} for degrees of freedom computation and
#' \code{\link{infcrit.dof}} for computing information criteria directly from a
#' previously fitted plsR model.
#' @references M. Hansen, B. Yu. (2001). Model Selection and Minimum Descripion
#' Length Principle, \emph{Journal of the American Statistical Association},
#' 96, 746-774.\cr N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of
#' Partial Least Squares Regression. \emph{Journal of the American Statistical
#' Association}, 106(494), 697-705.\cr N. Kraemer, M.L. Braun, Kernelizing PLS,
#' Degrees of Freedom, and Efficient Model Selection, \emph{Proceedings of the
#' 24th International Conference on Machine Learning}, Omni Press, (2007)
#' 441-448.
#' @keywords models regression utilities
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modpls <- plsR(yCornell,XCornell,4)
#' dof.object <- plsR.dof(modpls)
#' aic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
#' bic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
#' gmdl.dof(dof.object$sigmahat,modpls$nr,dof.object$DoF,dof.object$yhat)
#' naive.object <- plsR.dof(modpls,naive=TRUE)
#' aic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
#' bic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
#' gmdl.dof(naive.object$sigmahat,modpls$nr,naive.object$DoF,naive.object$yhat)
#' 
#' @export aic.dof
aic.dof <- function (RSS, n, DoF, sigmahat) 
{
    aic_temp <- RSS/n + 2 * (DoF/n) * sigmahat^2
    return(aic_temp)
}

