\name{simul_data_UniYX}
\alias{simul_data_UniYX}
\title{Data generating function for univariate plsR models}
\description{
This function generates a single univariate response value \eqn{Y} and a vector of explanatory variables \eqn{(X_1,\ldots,X_{totdim})} drawn from a model with a given number of latent components.
}
\usage{
simul_data_UniYX(totdim, ncomp)
}
\arguments{
  \item{totdim}{Number of columns of the X vector (from \code{ncomp} to hardware limits)}
  \item{ncomp}{Number of latent components in the model (from 2 to 6)}
}
\details{
This function should be combined with the replicate function to give rise to a larger dataset. The algorithm used is a \R{} port of the one described in the article of Li which is a multivariate generalization of the algorithm of Naes and Martens.
}
\value{
  \item{vector}{\eqn{(Y,X_1,\ldots,X_{totdim})}}
}
\references{
T. Naes, H. Martens, Comparison of prediction methods for 
multicollinear data, Commun. Stat., Simul. 14 (1985) 545-576.
\url{http://dx.doi.org/10.1080/03610918508812458}
 
Baibing Li, Julian Morris, Elaine B. Martin,
Model selection for partial least squares regression, 
Chemometrics and Intelligent Laboratory Systems 64 (2002) 79-89.
\url{http://dx.doi.org/10.1016/S0169-7439(02)00051-5}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
%\note{}
\seealso{\code{\link{simul_data_YX}} and \code{\link{simul_data_complete}} for generating multivariate data}
\examples{
simul_data_UniYX(20,6)                          

\donttest{
dimX <- 6
Astar <- 2

simul_data_UniYX(dimX,Astar)
(dataAstar2 <- t(replicate(50,simul_data_UniYX(dimX,Astar))))
library(plspm)
resAstar2 <- plsreg1(dataAstar2[,2:7],dataAstar2[,1], cv=TRUE)
resAstar2$Q2
plsR(dataAstar2[,1],dataAstar2[,2:7],5)


dimX <- 6
Astar <- 3

simul_data_UniYX(dimX,Astar)
(dataAstar3 <- t(replicate(50,simul_data_UniYX(dimX,Astar))))
resAstar3 <- plsreg1(dataAstar3[,2:7],dataAstar3[,1],cv=TRUE)
resAstar3$Q2
plsR(dataAstar3[,1],dataAstar3[,2:7],5)


dimX <- 6
Astar <- 4

simul_data_UniYX(dimX,Astar)
(dataAstar4 <- t(replicate(50,simul_data_UniYX(dimX,Astar))))
resAstar4 <- plsreg1(dataAstar4[,2:7],dataAstar4[,1],cv=TRUE)
resAstar4$Q2
plsR(dataAstar4[,1],dataAstar4[,2:7],5)

rm(list=c("dimX","Astar"))
}
}\keyword{datagen}
\keyword{utilities}
