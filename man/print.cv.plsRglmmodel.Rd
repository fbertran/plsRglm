\name{print.cv.plsRglmmodel}
\alias{print.cv.plsRglmmodel}
\title{Print method for plsRglm models}
\description{
This function provides a print method for the class \code{"cv.plsRglmmodel"}
}
\usage{
\method{print}{cv.plsRglmmodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"cv.plsRglmmodel"}}
  \item{\dots}{not used}
}
%\details{}
\value{
\code{NULL}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
print(cv.plsRglm(dataY=yCornell,dataX=XCornell,nt=10,NK=1,
modele="pls-glm-family",family=gaussian(), verbose=FALSE))
rm(list=c("XCornell","yCornell","bbb"))
}
\keyword{methods}
\keyword{print}
