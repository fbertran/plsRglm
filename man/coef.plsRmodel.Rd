\name{coef.plsRmodel}
\alias{coef.plsRmodel}
\title{coef method for plsR models}
\description{
This function provides a coef method for the class \code{"plsRmodel"}
}
\usage{
\method{coef}{plsRmodel}(object,type=c("scaled","original"), \dots)
}
\arguments{
  \item{object}{an object of the class \code{"plsRmodel"}}
  \item{type}{if \code{scaled}, the coefficients of the predictors are given for the scaled predictors, if \code{original} the coefficients are to be used with the predictors on their original scale.}
  \item{\dots}{not used}
}
\value{
   \item{CoeffC}{Coefficients of the components.}
   \item{Std.Coeffs}{Coefficients of the scaled predictors.}
   \item{Coeffs}{Coefficients of the untransformed predictors (on their original scale).}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{r?gression}{r\'egression} PLS et de la \enc{r?gression}{r\'egression} logistique PLS : application aux \enc{donn?es}{donn\'ees} \enc{d'all?lotypage}{d'all\'elotypage}. \emph{Journal de la Soci?t? Fran?aise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{coef}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
class(modpls)
coef(modpls)
coef(modpls,type="scaled")
rm(list=c("XCornell","yCornell","modpls"))
}
\keyword{methods}
\keyword{coef}
