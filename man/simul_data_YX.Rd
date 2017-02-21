\name{simul_data_YX}
\alias{simul_data_YX}
\title{Data generating function for multivariate plsR models}
\description{
This function generates a single multivariate response value \eqn{\boldsymbol{Y}} and a vector of explinatory variables \eqn{(X_1,\ldots,X_{totdim})} drawn from a model with a given number of latent components.
}
\usage{
simul_data_YX(totdim, ncomp)
}
\arguments{
  \item{totdim}{Number of column of the X vector (from \code{ncomp} to hardware limits)}
  \item{ncomp}{Number of latent components in the model (from 2 to 6)}
}
\details{
This function should be combined with the replicate function to give rise to a larger dataset. The algorithm used is a \R{} port of the one described in the article of Li which is a multivariate generalization of the algorithm of Naes and Martens.
}
\value{\item{vector}{\eqn{(Y_1,\ldots,Y_r,X_1,\ldots,X_{totdim})}}}
\references{
T. Naes, H. Martens, Comparison of prediction methods for 
multicollinear data, Commun. Stat., Simul. 14 (1985) 545-576.\cr
\url{http://dx.doi.org/10.1080/03610918508812458}\cr
 Baibing Li, Julian Morris, Elaine B. Martin,
Model selection for partial least squares regression, 
Chemometrics and Intelligent Laboratory Systems 64 (2002) 79-89.\cr
\url{http://dx.doi.org/10.1016/S0169-7439(02)00051-5}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{
The value of \eqn{r} depends on the value of \code{ncomp} :
\tabular{cc}{
\code{ncomp} \tab \eqn{r} \cr
2 \tab 3 \cr
3 \tab 3 \cr
4 \tab 4 \cr
}
}
\seealso{\code{\link{simul_data_complete}} for highlighting the simulations parameters}
\examples{
simul_data_YX(20,6)                          

\donttest{
if(require(plsdepot)){
dimX <- 6
Astar <- 2
(dataAstar2 <- t(replicate(50,simul_data_YX(dimX,Astar))))
library(plsdepot)
resAstar2 <- plsreg2(dataAstar2[,4:9],dataAstar2[,1:3],comps=5)
resAstar2$Q2
resAstar2$Q2[,4]>0.0975

dimX <- 6
Astar <- 3
(dataAstar3 <- t(replicate(50,simul_data_YX(dimX,Astar))))
library(plsdepot)
resAstar3 <- plsreg2(dataAstar3[,4:9],dataAstar3[,1:3],comps=5)
resAstar3$Q2
resAstar3$Q2[,4]>0.0975

dimX <- 6
Astar <- 4
(dataAstar4 <- t(replicate(50,simul_data_YX(dimX,Astar))))
library(plsdepot)
resAstar4 <- plsreg2(dataAstar4[,5:10],dataAstar4[,1:4],comps=5)
resAstar4$Q2
resAstar4$Q2[,5]>0.0975

rm(list=c("dimX","Astar"))
}
}
}\keyword{datagen}
\keyword{utilities}
