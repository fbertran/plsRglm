#' Information criteria
#' 
#' This function computes information criteria for existing plsR model using
#' Degrees of Freedom estimation.
#' 
#' If \code{naive=FALSE} returns AIC, BIC and gmdl values for estimated and
#' naive degrees of freedom. If \code{naive=TRUE} returns \code{NULL}.
#' 
#' @param modplsR A plsR model i.e. an object returned by one of the functions
#' \code{plsR}, \code{plsRmodel.default}, \code{plsRmodel.formula},
#' \code{PLS_lm} or \code{PLS_lm_formula}.
#' @param naive A boolean.
#' @return \item{matrix}{AIC, BIC and gmdl values or \code{NULL}.}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{plsR.dof}} for degrees of freedom computation and
#' \code{\link{infcrit.dof}} for computing information criteria directly from a
#' previously fitted plsR model.
#' @references M. Hansen, B. Yu. (2001). Model Selection and Minimum Descripion
#' Length Principle, \emph{Journal of the American Statistical Association},
#' 96, 746-774.\cr N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of
#' Partial Least Squares Regression. \emph{Journal of the American Statistical
#' Association}, 106(494), 697-705.\cr N. Kraemer, M. Sugiyama, M.L. Braun.
#' (2009). Lanczos Approximations for the Speedup of Kernel Partial Least
#' Squares Regression, \emph{Proceedings of the Twelfth International
#' Conference on Artificial Intelligence and Statistics (AISTATS)}, 272-279.
#' @keywords models regression utilities
#' @examples
#' 
#' data(Cornell)
#' XCornell<-Cornell[,1:7]
#' yCornell<-Cornell[,8]
#' modpls <- plsR(yCornell,XCornell,4)
#' infcrit.dof(modpls)
#' 
#' @export infcrit.dof
infcrit.dof <-
  function (modplsR, naive = FALSE) 
  {
    if (!(is.null(modplsR$weights) | identical(modplsR$weights, 
                                               rep(1L, modplsR$nr)) | identical(modplsR$weights, rep(1, 
                                                                                                     modplsR$nr)))) {
      naive = TRUE
    }
    if (modplsR$na.miss.X | modplsR$na.miss.Y) {
      naive = TRUE
    }
    if (!naive) {
      tempmodplsR_dof <- plsR.dof(modplsR, naive = FALSE)
      tempAIC.dof <- aic.dof(modplsR$RSS, modplsR$nr, tempmodplsR_dof$DoF, 
                             tempmodplsR_dof$sigmahat)
      tempBIC.dof <- bic.dof(modplsR$RSS, modplsR$nr, tempmodplsR_dof$DoF, 
                             tempmodplsR_dof$sigmahat)
      tempGMDL.dof <- gmdl.dof(tempmodplsR_dof$sigmahat, modplsR$nr, 
                               tempmodplsR_dof$DoF, tempmodplsR_dof$yhat)
      tempmodplsR_naive <- plsR.dof(modplsR, naive = TRUE)
      tempAIC.naive <- aic.dof(modplsR$RSS, modplsR$nr, tempmodplsR_naive$DoF, 
                               tempmodplsR_naive$sigmahat)
      tempBIC.naive <- bic.dof(modplsR$RSS, modplsR$nr, tempmodplsR_naive$DoF, 
                               tempmodplsR_naive$sigmahat)
      tempGMDL.naive <- gmdl.dof(tempmodplsR_naive$sigmahat, 
                                 modplsR$nr, tempmodplsR_naive$DoF, tempmodplsR_naive$yhat)
      InfCrit.dof <- t(rbind(tempmodplsR_dof$DoF, tempmodplsR_dof$sigmahat, 
                             tempAIC.dof, tempBIC.dof, tempGMDL.dof, tempmodplsR_naive$DoF, 
                             tempmodplsR_naive$sigmahat, tempAIC.naive, tempBIC.naive, 
                             tempGMDL.naive))
      dimnames(InfCrit.dof) <- list(paste("Nb_Comp_", 0:modplsR$computed_nt, 
                                          sep = ""), c("DoF.dof", "sigmahat.dof", "AIC.dof", 
                                                       "BIC.dof", "GMDL.dof", "DoF.naive", "sigmahat.naive", 
                                                       "AIC.naive", "BIC.naive", "GMDL.naive"))
    }
    else {
      tempmodplsR_naive <- plsR.dof(modplsR, naive = TRUE)
      tempAIC.naive <- aic.dof(modplsR$RSS, modplsR$nr, tempmodplsR_naive$DoF, 
                               tempmodplsR_naive$sigmahat)
      tempBIC.naive <- bic.dof(modplsR$RSS, modplsR$nr, tempmodplsR_naive$DoF, 
                               tempmodplsR_naive$sigmahat)
      tempGMDL.naive <- gmdl.dof(tempmodplsR_naive$sigmahat, 
                                 modplsR$nr, tempmodplsR_naive$DoF, tempmodplsR_naive$yhat)
      InfCrit.dof <- t(rbind(NA, NA, NA, NA, NA, tempmodplsR_naive$DoF, 
                             tempmodplsR_naive$sigmahat, tempAIC.naive, tempBIC.naive, 
                             tempGMDL.naive))
      dimnames(InfCrit.dof) <- list(paste("Nb_Comp_", 0:modplsR$computed_nt, 
                                          sep = ""), c("DoF.dof", "sigmahat.dof", "AIC.dof", 
                                                       "BIC.dof", "GMDL.dof", "DoF.naive", "sigmahat.naive", 
                                                       "AIC.naive", "BIC.naive", "GMDL.naive"))
      
    }
    return(InfCrit.dof)
  }
