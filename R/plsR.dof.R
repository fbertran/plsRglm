#' Computation of the Degrees of Freedom
#' 
#' This function computes the Degrees of Freedom using the Krylov
#' representation of PLS and other quantities that are used to get information
#' criteria values. For the time present, it only works with complete datasets.
#' 
#' If \code{naive=FALSE} returns values for estimated degrees of freedom and
#' error dispersion. If \code{naive=TRUE} returns returns values for naive
#' degrees of freedom and error dispersion. The original code from Nicole
#' Kraemer and Mikio L. Braun was unable to handle models with only one
#' component.
#' 
#' @param modplsR A plsR model i.e. an object returned by one of the functions
#' \code{plsR}, \code{plsRmodel.default}, \code{plsRmodel.formula},
#' \code{PLS_lm} or \code{PLS_lm_formula}.
#' @param naive A boolean.
#' @return \item{DoF}{Degrees of Freedom} \item{sigmahat}{Estimates of
#' dispersion} \item{Yhat}{Predicted values} \item{yhat}{Square Euclidean norms
#' of the predicted values} \item{RSS}{Residual Sums of Squares}
#' @author Nicole Kraemer, Mikio L. Braun with improvements from
#' Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{aic.dof}} and \code{\link{infcrit.dof}} for computing
#' information criteria directly from a previously fitted plsR model.
#' @references N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of
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
#' plsR.dof(modpls) 
#' plsR.dof(modpls,naive=TRUE) 
#' 
#' @export plsR.dof
plsR.dof <- function(modplsR,naive=FALSE){
temp.object <- vector("list",0)
dof.object <- vector("list",0)
temp.object$RSS <- as.vector(modplsR$RSS)
temp.object$Yhat <- matrix(mean(modplsR$dataY),nrow=modplsR$nr,ncol=1)
for(i in 1:modplsR$computed_nt){
temp.object$Yhat <- cbind(temp.object$Yhat,attr(modplsR$RepY, "scaled:center") + attr(modplsR$RepY, "scaled:scale") * modplsR$tt[,1:i,drop=FALSE] %*% modplsR$CoeffC[1:i])
}
colnames(temp.object$Yhat) <- paste("Nt_",0:modplsR$computed_nt,sep="")
if(!naive){
temp.object$TT <- sweep(modplsR$tt,MARGIN=2,FUN="/",sqrt(diag(crossprod(modplsR$tt))))

requireNamespace("plsdof")
pls.doftemp <- function (pls.object, n, y, K, m, DoF.max) 
{
    TT <- pls.object$TT
    Yhat <- pls.object$Yhat[, 2:(m + 1), drop=FALSE]
    TK = matrix(, m, m)
    KY <- plsdof::krylov(K, K %*% y, m)
    lambda <- eigen(K)$values
    tr.K <- vector(length = m)
    for (i in 1:m) {
        tr.K[i] <- sum(lambda^i)
    }
    BB = t(TT) %*% KY
    BB[row(BB) > col(BB)] = 0
    b <- t(TT) %*% y
    DoF = vector(length = m)
    Binv <- backsolve(BB, diag(m))
    tkt <- rep(0, m)
    ykv <- rep(0, m)
    KjT <- array(dim = c(m, n, m))
    dummy <- TT
    for (i in 1:m) {
        dummy <- K %*% dummy
        KjT[i, , ] <- dummy
    }
    trace.term = rep(0, m)
    for (i in 1:m) {
        Binvi <- Binv[1:i, 1:i, drop = FALSE]
        ci <- Binvi %*% b[1:i]
        Vi <- TT[, 1:i, drop = FALSE] %*% t(Binvi)
        trace.term[i] <- sum(ci * tr.K[1:i])
        ri <- y - Yhat[, i]
        for (j in 1:i) {
            KjTj = KjT[j, , ]
            if(is.null(dim(KjTj))){KjTj <- matrix(KjTj,ncol=1)}
            tkt[i] <- tkt[i] + ci[j] * plsdof::tr(t(TT[, 1:i, drop = FALSE]) %*% 
                KjTj[, 1:i, drop = FALSE])
            ri <- K %*% ri
            ykv[i] <- ykv[i] + sum(ri * Vi[, j])
        }
    }
    DoF <- trace.term + 1:m - tkt + ykv
    DoF[DoF > DoF.max] = DoF.max
    sigmahat = sqrt(pls.object$RSS[-1]/(n - DoF))
    return(list(DoF = DoF, sigmahat = sigmahat))
}

dof.object <- pls.doftemp(temp.object,n=modplsR$nr,y=modplsR$dataY,K=modplsR$ExpliX%*%t(modplsR$ExpliX),m=modplsR$computed_nt,DoF.max=min(modplsR$nr-1,modplsR$nc+1)-1)
temp.object$sigmahat <- c(sqrt(temp.object$RSS[1]/(modplsR$nr-1)), sqrt(temp.object$RSS[-1]/(modplsR$nr-dof.object$DoF)))
temp.object$DoF <- c(0, dof.object$DoF) + 1
} else
{
temp.object$DoF <- 1:(modplsR$computed_nt+1)
temp.object$sigmahat <- sqrt(temp.object$RSS/(modplsR$nr - temp.object$DoF))
}
temp.object$yhat <- vector("numeric",length=modplsR$computed_nt+1)
for(i in 1:(modplsR$computed_nt+1)){
temp.object$yhat[i] <- sum((temp.object$Yhat[,i])^2)
}
return(list(DoF=temp.object$DoF, sigmahat=temp.object$sigmahat, Yhat=temp.object$Yhat, yhat=temp.object$yhat, RSS=temp.object$RSS))
}
