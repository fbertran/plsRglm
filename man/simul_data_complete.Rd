\name{simul_data_complete}
\alias{simul_data_complete}
\title{Data generating detailed process for multivariate plsR models}
\description{
This function generates a single multivariate response value \eqn{\boldsymbol{Y}} and a vector of explinatory variables \eqn{(X_1,\ldots,X_{totdim})} drawn from a model with a given number of latent components.
}
\usage{
simul_data_complete(totdim, ncomp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{totdim}{Number of columns of the X vector (from \code{ncomp} to hardware limits)}
  \item{ncomp}{Number of latent components in the model (from 2 to 6)}
}
\details{
This function should be combined with the replicate function to give rise to a larger dataset. The algorithm used is a \R port of the one described in the article of Li which is a multivariate generalization of the algorithm of Naes and Martens.
}
\value{
  \item{simX}{Vector of explanatory variables}
  \item{HH}{Dimension of the response \eqn{\boldsymbol{Y}}}
  \item{eta}{See Li et al.}
  \item{r}{See Li et al.}
  \item{epsilon}{See Li et al.}
  \item{ksi}{See Li et al.}
  \item{f}{See Li et al.}
  \item{z}{See Li et al.}
  \item{Y}{See Li et al.}
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
\note{
The value of \eqn{r} depends on the value of \code{ncomp} :
\tabular{cc}{
\code{ncomp} \tab \eqn{r} \cr
2 \tab 3 \cr
3 \tab 3 \cr
4 \tab 4 \cr
}
}
\seealso{\code{\link{simul_data_YX}} for data simulation purpose}
\examples{
simul_data_complete(20,6)                          

dimX <- 6
Astar <- 2
simul_data_complete(dimX,Astar)


dimX <- 6
Astar <- 3
simul_data_complete(dimX,Astar)


dimX <- 6
Astar <- 4
simul_data_complete(dimX,Astar)

rm(list=c("dimX","Astar"))
}
\keyword{datagen}
\keyword{utilities}
