#' @keywords internal
#'
#' @name plsRglm-package
#' @references A short paper that sums up some of features of the package is available on \url{https://arxiv.org/}, Frédéric Bertrand and Myriam Maumy-Bertrand (2018), "plsRglm: Partial least squares linear and generalized linear regression for processing incomplete datasets by cross-validation and bootstrap techniques with R", *arxiv*, \url{https://arxiv.org/abs/1810.01005}, \url{https://github.com/fbertran/plsRglm/} et \url{https://fbertran.github.io/plsRglm/}
#' 
#' @examples
#' set.seed(314)
#' library(plsRglm)
#' data(Cornell)
#' cv.modpls<-cv.plsR(Y~.,data=Cornell,nt=6,K=6)
#' res.cv.modpls<-cvtable(summary(cv.modpls))
#' 
"_PACKAGE"

#' @importFrom graphics abline arrows axis barplot boxplot legend par plot points strwidth text
#' @importFrom stats AIC Gamma binomial coef contrasts family gaussian glm inverse.gaussian is.empty.model lm make.link model.frame model.matrix model.offset model.response model.weights na.exclude na.pass pnorm poisson predict rbinom residuals.glm weighted.mean
#' @importFrom MASS polr
#' @importFrom mvtnorm rmvnorm
#' @importFrom bipartite visweb
#' @importFrom car dataEllipse
#' @import boot 
NULL


