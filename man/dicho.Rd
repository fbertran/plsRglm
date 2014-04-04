\name{dicho}
\alias{dicho}
\title{Dichotomization}
\description{
This function takes a real value and converts it to 1 if it is positive and else to 0.
}
\usage{
dicho(val)
}
\arguments{
  \item{val}{A real value}%\code{val} 
}
%\details{
%  ~~ More details than the description above ~~
%}
\value{
0 or 1.
}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
%\note{ }
\seealso{\code{\link{ifelse}}}
\examples{
dimX <- 6
Astar <- 4
(dataAstar4 <- t(replicate(10,simul_data_YX(dimX,Astar))))

dicho(dataAstar4)

rm(list=c("dimX","Astar"))
}
\keyword{utilities}
