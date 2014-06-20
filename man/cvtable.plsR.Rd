\name{cvtable.plsR}
\alias{cvtable.plsR}
\title{Table method for summary of cross validated plsR models}
\description{
This function provides a table summary for the class \code{"summary.cv.plsRmodel"}
}
\usage{
cvtable.plsR(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"summary.cv.plsRmodel"}}
  \item{\dots}{further arguments to be passed to or from methods.}  
}
%\details{}
\value{
  \code{list}List of IC computed for each folds.
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{r?gression}{r\'egression} PLS et de la \enc{r?gression}{r\'egression} logistique PLS : application aux \enc{donn?es}{donn\'ees} \enc{d'all?lotypage}{d'all\'elotypage}. \emph{Journal de la Soci?t? Fran?aise de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{summary}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsR(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
cvtable.plsR(summary(bbb))
rm(list=c("XCornell","yCornell","bbb"))
}
\keyword{methods}
\keyword{print}
