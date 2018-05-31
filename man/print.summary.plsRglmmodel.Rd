\name{print.summary.plsRglmmodel}
\alias{print.summary.plsRglmmodel}
\title{Print method for summaries of plsRglm models}
\description{
This function provides a print method for the class \code{"summary.plsRglmmodel"}
}
\usage{
\method{print}{summary.plsRglmmodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"summary.plsRglmmodel"}}
  \item{\dots}{not used}
}
%\details{}
\value{
  \item{language}{call of the model}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{r?gression}{r\'egression} PLS et de la \enc{r?gression}{r\'egression} logistique PLS : application aux \enc{donn?es}{donn\'ees} \enc{d'all?lotypage}{d'all\'elotypage}. \emph{Journal de la Soci?t? Fran?aise de Statistique}, 151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}

\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}} and \code{\link{summary}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
class(modplsglm)
print(summary(modplsglm))
rm(list=c("XCornell","yCornell","modplsglm"))
}
\keyword{methods}
\keyword{print}
