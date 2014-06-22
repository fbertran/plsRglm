\name{print.coef.plsRmodel}
\alias{print.coef.plsRmodel}
\title{Print method for plsR models}
\description{
This function provides a print method for the class \code{"coef.plsRmodel"}
}
\usage{
\method{print}{coef.plsRmodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"coef.plsRmodel"}}
  \item{\dots}{not used}
}
\value{
\code{NULL}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{r?gression}{r\'egression} PLS et de la \enc{r?gression}{r\'egression} logistique PLS : application aux \enc{donn?es}{donn\'ees} \enc{d'all?lotypage}{d'all\'elotypage}. \emph{Journal de la Soci?t? Fran?aise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}

\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
class(modpls)
print(coef(modpls))
rm(list=c("XCornell","yCornell","modpls"))
}
\keyword{methods}
\keyword{print}
