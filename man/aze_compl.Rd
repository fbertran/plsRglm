\name{aze_compl}
\alias{aze_compl}
\docType{data}
\title{As aze without missing values}
\description{
This is a single imputation of the \code{\link{aze}} dataset which was collected on patients carrying a colon adenocarcinoma. It has 104 observations on 33 binary qualitative explanatory variables and one response variable \code{y} representing the cancer stage according to the to Astler-Coller classification (Astler and Coller, 1954). A microsattelite is a non-coding DNA sequence.
}
\usage{data(aze_compl)}
\format{
  A data frame with 104 observations on the following 34 variables.
  \describe{
    \item{\code{y}}{the response: a binary vector (Astler-Coller score).}
    \item{\code{D2S138}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D18S61}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D16S422}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D17S794}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D6S264}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D14S65}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D18S53}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D17S790}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D1S225}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D3S1282}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D9S179}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D5S430}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D8S283}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D11S916}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D2S159}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D16S408}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D5S346}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D10S191}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D13S173}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D6S275}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D15S127}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D1S305}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D4S394}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D20S107}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D1S197}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D1S207}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D10S192}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D3S1283}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D4S414}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D8S264}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D22S928}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{TP53}}{a binary vector that indicates whether this microsatellite is altered or not.}
    \item{\code{D9S171}}{a binary vector that indicates whether this microsatellite is altered or not.}
  }
}
\source{
Weber \emph{et al.} (2007). Allelotyping analyzes of synchronous primary and metastasis CIN colon cancers identified different subtypes. \emph{Int J Cancer}, 120(3), pages 524-32.
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frederic}{Fr\'ed\'eric} Bertrand (2010). Comparing the linear and the logistic PLS regression with qualitative predictors: application to allelotyping data. \emph{Journal de la Societe Francaise de Statistique}, 151(2), pages 1-18.
}
\examples{
data(aze_compl)
str(aze_compl)
}
\keyword{datasets}
