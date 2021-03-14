#' As aze without missing values
#' 
#' This is a single imputation of the \code{\link{aze}} dataset which was
#' collected on patients carrying a colon adenocarcinoma. It has 104
#' observations on 33 binary qualitative explanatory variables and one response
#' variable \code{y} representing the cancer stage according to the to
#' Astler-Coller classification (Astler and Coller, 1954). A microsattelite is
#' a non-coding DNA sequence.
#' 
#' 
#' @name aze_compl
#' @docType data
#' @format A data frame with 104 observations on the following 34 variables.
#' \describe{ \item{y}{the response: a binary vector (Astler-Coller
#' score).} \item{D2S138}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D18S61}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D16S422}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D17S794}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D6S264}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D14S65}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D18S53}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D17S790}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D1S225}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D3S1282}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D9S179}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D5S430}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D8S283}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D11S916}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D2S159}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D16S408}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D5S346}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D10S191}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D13S173}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D6S275}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D15S127}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D1S305}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D4S394}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D20S107}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D1S197}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D1S207}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D10S192}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D3S1283}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D4S414}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D8S264}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D22S928}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{TP53}{a binary vector that
#' indicates whether this microsatellite is altered or not.}
#' \item{D9S171}{a binary vector that indicates whether this
#' microsatellite is altered or not.} }
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Société Française de Statistique},
#' 151(2), pages 1-18.
#' @source Weber \emph{et al.} (2007). Allelotyping analyzes of synchronous
#' primary and metastasis CIN colon cancers identified different subtypes.
#' \emph{Int J Cancer}, 120(3), pages 524-32.
#' @keywords datasets
#' @examples
#' 
#' data(aze_compl)
#' str(aze_compl)
#' 
NULL





#' Microsatellites Dataset
#' 
#' This database was collected on patients carrying a colon adenocarcinoma. It
#' has 104 observations on 33 binary qualitative explanatory variables and one
#' response variable \code{y} representing the cancer stage according to the to
#' Astler-Coller classification (Astler and Coller, 1954). This dataset has
#' some missing data due to technical limits. A microsattelite is a non-coding
#' DNA sequence.
#' 
#' 
#' @name aze
#' @docType data
#' @format A data frame with 104 observations on the following 34 variables.
#' \describe{ \item{y}{the response: a binary vector (Astler-Coller
#' score).} \item{D2S138}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D18S61}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D16S422}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D17S794}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D6S264}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D14S65}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D18S53}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D17S790}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D1S225}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D3S1282}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D9S179}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D5S430}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D8S283}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D11S916}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D2S159}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D16S408}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D5S346}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D10S191}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D13S173}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D6S275}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D15S127}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D1S305}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D4S394}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D20S107}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D1S197}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D1S207}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D10S192}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D3S1283}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D4S414}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{D8S264}{a binary vector
#' that indicates whether this microsatellite is altered or not.}
#' \item{D22S928}{a binary vector that indicates whether this
#' microsatellite is altered or not.} \item{TP53}{a binary vector that
#' indicates whether this microsatellite is altered or not.}
#' \item{D9S171}{a binary vector that indicates whether this
#' microsatellite is altered or not.} }
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Société Française de Statistique},
#' 151(2), pages 1-18.
#' @source Weber \emph{et al.} (2007). Allelotyping analyzes of synchronous
#' primary and metastasis CIN colon cancers identified different subtypes.
#' \emph{Int J Cancer}, 120(3), pages 524-32.
#' @keywords datasets
#' @examples
#' 
#' data(aze)
#' str(aze)
#' 
NULL





#' Quality of wine dataset
#' 
#' Quality of Bordeaux wines (\code{Quality}) and four potentially predictive
#' variables (\code{Temperature}, \code{Sunshine}, \code{Heat} and
#' \code{Rain}).
#' 
#' 
#' @name bordeaux
#' @docType data
#' @format A data frame with 34 observations on the following 5 variables.
#' \describe{ \item{Temperature}{a numeric vector}
#' \item{Sunshine}{a numeric vector} \item{Heat}{a numeric
#' vector} \item{Rain}{a numeric vector} \item{Quality}{an
#' ordered factor with levels \code{1} < \code{2} < \code{3}} }
#' @references M. Tenenhaus. (2005). La regression logistique PLS. In J.-J.
#' Droesbeke, M. Lejeune, and G. Saporta, editors, Modeles statistiques pour
#' donnees qualitatives. Editions Technip, Paris.
#' @source P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS
#' generalised linear regression. \emph{Computational Statistics & Data
#' Analysis}, 48(1):17-46.
#' @keywords datasets
#' @examples
#' 
#' data(bordeaux)
#' str(bordeaux)
#' 
NULL





#' Quality of wine dataset
#' 
#' Quality of Bordeaux wines (\code{Quality}) and four potentially predictive
#' variables (\code{Temperature}, \code{Sunshine}, \code{Heat} and
#' \code{Rain}).
#' 
#' The value of x1 for the first observation was removed from the matrix of
#' predictors on purpose.
#' 
#' The \code{bordeauxNA} is a dataset with a missing value for testing purpose.
#' 
#' @name bordeauxNA
#' @docType data
#' @format A data frame with 34 observations on the following 5 variables.
#' \describe{ \item{Temperature}{a numeric vector}
#' \item{Sunshine}{a numeric vector} \item{Heat}{a numeric
#' vector} \item{Rain}{a numeric vector} \item{Quality}{an
#' ordered factor with levels \code{1} < \code{2} < \code{3}} }
#' @references M. Tenenhaus. (2005). La regression logistique PLS. In J.-J.
#' Droesbeke, M. Lejeune, and G. Saporta, editors, Modeles statistiques pour
#' donnees qualitatives. Editions Technip, Paris.
#' @source P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS
#' generalised linear regression. \emph{Computational Statistics & Data
#' Analysis}, 48(1):17-46.
#' @keywords datasets
#' @examples
#' 
#' data(bordeauxNA)
#' str(bordeauxNA)
#' 
NULL





#' Correlation matrix for simulating plsR datasets
#' 
#' A correlation matrix to simulate datasets
#' 
#' 
#' @name CorMat
#' @docType data
#' @format A data frame with 17 observations on the following 17 variables.
#' \describe{ \item{y}{a numeric vector} \item{x11}{a numeric
#' vector} \item{x12}{a numeric vector} \item{x13}{a numeric
#' vector} \item{x21}{a numeric vector} \item{x22}{a numeric
#' vector} \item{x31}{a numeric vector} \item{x32}{a numeric
#' vector} \item{x33}{a numeric vector} \item{x34}{a numeric
#' vector} \item{x41}{a numeric vector} \item{x42}{a numeric
#' vector} \item{x51}{a numeric vector} \item{x61}{a numeric
#' vector} \item{x62}{a numeric vector} \item{x63}{a numeric
#' vector} \item{x64}{a numeric vector} }
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @source Handmade.
#' @keywords datasets
#' @examples
#' 
#' data(CorMat)
#' str(CorMat)
#' 
NULL





#' Cornell dataset
#' 
#' The famous Cornell dataset. A mixture experiment on \code{X1}, \code{X2},
#' \code{X3}, \code{X4}, \code{X5}, \code{X6} and \code{X7} to analyse octane
#' degree (\code{Y}) in gazoline.
#' 
#' 
#' @name Cornell
#' @docType data
#' @format A data frame with 12 observations on the following 8 variables.
#' \describe{ \item{X1}{a numeric vector} \item{X2}{a numeric
#' vector} \item{X3}{a numeric vector} \item{X4}{a numeric
#' vector} \item{X5}{a numeric vector} \item{X6}{a numeric
#' vector} \item{X7}{a numeric vector} \item{Y}{response value:
#' a numeric vector} }
#' @references N. Kettaneh-Wold. Analysis of mixture data with partial least
#' squares. (1992). \emph{Chemometrics and Intelligent Laboratory Systems},
#' 14(1):57-69.
#' @source M. Tenenhaus. (1998). \emph{La regression PLS, Theorie et pratique}.
#' Editions Technip, Paris.
#' @keywords datasets
#' @examples
#' 
#' data(Cornell)
#' str(Cornell)
#' 
NULL





#' Fowlkes dataset
#' 
#' A classic dataset from Fowlkes.
#' 
#' 
#' @name fowlkes
#' @docType data
#' @format A data frame with 9949 observations on the following 13 variables.
#' \describe{ \item{Y}{binary response} \item{MA}{a numeric
#' vector} \item{MW}{a numeric vector} \item{NE}{a numeric
#' vector} \item{NW}{a numeric vector} \item{PA}{a numeric
#' vector} \item{SO}{a numeric vector} \item{SW}{a numeric
#' vector} \item{color}{a numeric vector} \item{age1}{a numeric
#' vector} \item{age2}{a numeric vector} \item{age3}{a numeric
#' vector} \item{sexe}{a numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(fowlkes)
#' str(fowlkes)
#' 
NULL





#' Complete Pine dataset
#' 
#' This is the complete caterpillar dataset from a 1973 study on pine_full
#' processionary caterpillars. It assesses the influence of some forest
#' settlement characteristics on the development of caterpillar colonies. The
#' response variable is the logarithmic transform of the average number of
#' nests of caterpillars per tree in an area of 500 square meters (\code{x11}).
#' There are k=10 potentially explanatory variables defined on n=55 areas.
#' 
#' These caterpillars got their names from their habit of moving over the
#' ground in incredibly long head-to-tail processions when leaving their nest
#' to create a new colony.
#' 
#' @name pine_full
#' @docType data
#' @format A data frame with 55 observations on the following 11 variables.
#' \describe{ \item{x1}{altitude (in meters)} \item{x2}{slope
#' (en degrees)} \item{x3}{number of pine_fulls in the area}
#' \item{x4}{height (in meters) of the tree sampled at the center of
#' the area} \item{x5}{diameter (in meters) of the tree sampled at the
#' center of the area} \item{x6}{index of the settlement density}
#' \item{x7}{orientation of the area (from 1 if southbound to 2
#' otherwise)} \item{x8}{height (in meters) of the dominant tree}
#' \item{x9}{number of vegetation strata} \item{x10}{mix
#' settlement index (from 1 if not mixed to 2 if mixed)}
#' \item{x11}{logarithmic transform of the average number of nests of
#' caterpillars per tree} }
#' @references J.-M. Marin, C. Robert. (2007). \emph{Bayesian Core: A Practical
#' Approach to Computational Bayesian Statistics}. Springer, New-York, pages
#' 48-49.
#' @source Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992),
#' \dQuote{La régression, nouveaux regards sur une ancienne
#' méthode statistique}, INRA,
#' \emph{Actualités Scientifiques et Agronomiques}, Masson,
#' Paris.
#' @keywords datasets
#' @examples
#' 
#' data(pine_full)
#' str(pine_full)
#' 
NULL





#' Complete Pine dataset
#' 
#' This is a supplementary dataset (used as a test set for the \code{pine}
#' dataset) that was extracted from a 1973 study on pine_sup processionary
#' caterpillars. It assesses the influence of some forest settlement
#' characteristics on the development of caterpillar colonies. The response
#' variable is the logarithmic transform of the average number of nests of
#' caterpillars per tree in an area of 500 square meters (\code{x11}). There
#' are k=10 potentially explanatory variables defined on n=22 areas.
#' 
#' These caterpillars got their names from their habit of moving over the
#' ground in incredibly long head-to-tail processions when leaving their nest
#' to create a new colony.\cr
#' 
#' The \code{pine_sup} dataset can be used as a test set to assess model
#' prediction error of a model trained on the \code{pine} dataset.
#' 
#' @name pine_sup
#' @docType data
#' @format A data frame with 22 observations on the following 11 variables.
#' \describe{ \item{x1}{altitude (in meters)} \item{x2}{slope
#' (en degrees)} \item{x3}{number of pine_sups in the area}
#' \item{x4}{height (in meters) of the tree sampled at the center of
#' the area} \item{x5}{diameter (in meters) of the tree sampled at the
#' center of the area} \item{x6}{index of the settlement density}
#' \item{x7}{orientation of the area (from 1 if southbound to 2
#' otherwise)} \item{x8}{height (in meters) of the dominant tree}
#' \item{x9}{number of vegetation strata} \item{x10}{mix
#' settlement index (from 1 if not mixed to 2 if mixed)}
#' \item{x11}{logarithmic transform of the average number of nests of
#' caterpillars per tree} }
#' @references J.-M. Marin, C. Robert. (2007). \emph{Bayesian Core: A Practical
#' Approach to Computational Bayesian Statistics}. Springer, New-York, pages
#' 48-49.
#' @source Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992),
#' \dQuote{La régression, nouveaux regards sur une ancienne
#' méthode statistique}, INRA,
#' \emph{Actualités Scientifiques et Agronomiques}, Masson,
#' Paris.
#' @keywords datasets
#' @examples
#' 
#' data(pine_sup)
#' str(pine_sup)
#' 
NULL





#' Pine dataset
#' 
#' The caterpillar dataset was extracted from a 1973 study on pine
#' processionary caterpillars. It assesses the influence of some forest
#' settlement characteristics on the development of caterpillar colonies. The
#' response variable is the logarithmic transform of the average number of
#' nests of caterpillars per tree in an area of 500 square meters (\code{x11}).
#' There are k=10 potentially explanatory variables defined on n=33 areas.
#' 
#' These caterpillars got their names from their habit of moving over the
#' ground in incredibly long head-to-tail processions when leaving their nest
#' to create a new colony.\cr
#' 
#' The \code{pine_sup} dataset can be used as a test set to assess model
#' prediction error of a model trained on the \code{pine} dataset.
#' 
#' @name pine
#' @docType data
#' @format A data frame with 33 observations on the following 11 variables.
#' \describe{ \item{x1}{altitude (in meters)} \item{x2}{slope
#' (en degrees)} \item{x3}{number of pines in the area}
#' \item{x4}{height (in meters) of the tree sampled at the center of
#' the area} \item{x5}{diameter (in meters) of the tree sampled at the
#' center of the area} \item{x6}{index of the settlement density}
#' \item{x7}{orientation of the area (from 1 if southbound to 2
#' otherwise)} \item{x8}{height (in meters) of the dominant tree}
#' \item{x9}{number of vegetation strata} \item{x10}{mix
#' settlement index (from 1 if not mixed to 2 if mixed)}
#' \item{x11}{logarithmic transform of the average number of nests of
#' caterpillars per tree} }
#' @references J.-M. Marin, C. Robert. (2007). \emph{Bayesian Core: A Practical
#' Approach to Computational Bayesian Statistics}. Springer, New-York, pages
#' 48-49.
#' @source Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992),
#' \dQuote{La régression, nouveaux regards sur une ancienne
#' méthode statistique}, INRA,
#' \emph{Actualités Scientifiques et Agronomiques}, Masson,
#' Paris.
#' @keywords datasets
#' @examples
#' 
#' data(pine)
#' str(pine)
#' 
NULL





#' Incomplete dataset from the pine caterpillars example
#' 
#' The caterpillar dataset was extracted from a 1973 study on pine
#' processionary caterpillars. It assesses the influence of some forest
#' settlement characteristics on the development of caterpillar colonies. There
#' are k=10 potentially explanatory variables defined on n=33 areas.\cr The
#' value of x2 for the first observation was removed from the matrix of
#' predictors on purpose.
#' 
#' These caterpillars got their names from their habit of moving over the
#' ground in incredibly long head-to-tail processions when leaving their nest
#' to create a new colony.\cr The \code{pineNAX21} is a dataset with a missing
#' value for testing purpose.
#' 
#' @name pineNAX21
#' @docType data
#' @format A data frame with 33 observations on the following 11 variables and
#' one missing value.  \describe{ \item{x1}{altitude (in meters)}
#' \item{x2}{slope (en degrees)} \item{x3}{number of pines in
#' the area} \item{x4}{height (in meters) of the tree sampled at the
#' center of the area} \item{x5}{diameter (in meters) of the tree
#' sampled at the center of the area} \item{x6}{index of the settlement
#' density} \item{x7}{orientation of the area (from 1 if southbound to
#' 2 otherwise)} \item{x8}{height (in meters) of the dominant tree}
#' \item{x9}{number of vegetation strata} \item{x10}{mix
#' settlement index (from 1 if not mixed to 2 if mixed)}
#' \item{x11}{logarithmic transform of the average number of nests of
#' caterpillars per tree} }
#' @source Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992).
#' \dQuote{La régression, nouveaux regards sur une ancienne
#' méthode statistique}, INRA,
#' \emph{Actualités Scientifiques et Agronomiques}, Masson,
#' Paris.
#' @keywords datasets
#' @examples
#' 
#' data(pineNAX21)
#' str(pineNAX21)
#' 
NULL





#' Incomplete dataset for the quality of wine dataset
#' 
#' Quality of Bordeaux wines (\code{Quality}) and four potentially predictive
#' variables (\code{Temperature}, \code{Sunshine}, \code{Heat} and
#' \code{Rain}).\cr The value of Temperature for the first observation was
#' remove from the matrix of predictors on purpose.
#' 
#' 
#' @name XbordeauxNA
#' @docType data
#' @format A data frame with 34 observations on the following 4 variables.
#' \describe{ \item{Temperature}{a numeric vector}
#' \item{Sunshine}{a numeric vector} \item{Heat}{a numeric
#' vector} \item{Rain}{a numeric vector} }
#' @references M. Tenenhaus. (2005). La regression logistique PLS. In J.-J.
#' Droesbeke, M. Lejeune, and G. Saporta, editors, Modeles statistiques pour
#' donnees qualitatives. Editions Technip, Paris.
#' @source P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS
#' generalised linear regression. \emph{Computational Statistics & Data
#' Analysis}, 48(1):17-46.
#' @keywords datasets
#' @examples
#' 
#' data(XbordeauxNA)
#' str(XbordeauxNA)
#' 
NULL





#' Incomplete dataset from the pine caterpillars example
#' 
#' The caterpillar dataset was extracted from a 1973 study on pine
#' processionary caterpillars. It assesses the influence of some forest
#' settlement characteristics on the development of caterpillar colonies. There
#' are k=10 potentially explanatory variables defined on n=33 areas.\cr The
#' value of x2 for the first observation was remove from the matrix of
#' predictors on purpose.
#' 
#' These caterpillars got their names from their habit of moving over the
#' ground in incredibly long head-to-tail processions when leaving their nest
#' to create a new colony.\cr The \code{XpineNAX21} is a dataset with a missing
#' value for testing purpose.
#' 
#' @name XpineNAX21
#' @docType data
#' @format A data frame with 33 observations on the following 10 variables and
#' one missing value.  \describe{ \item{x1}{altitude (in meters)}
#' \item{x2}{slope (en degrees)} \item{x3}{number of pines in
#' the area} \item{x4}{height (in meters) of the tree sampled at the
#' center of the area} \item{x5}{diameter (in meters) of the tree
#' sampled at the center of the area} \item{x6}{index of the settlement
#' density} \item{x7}{orientation of the area (from 1 if southbound to
#' 2 otherwise)} \item{x8}{height (in meters) of the dominant tree}
#' \item{x9}{number of vegetation strata} \item{x10}{mix
#' settlement index (from 1 if not mixed to 2 if mixed)} }
#' @source Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992).
#' \dQuote{La régression, nouveaux regards sur une ancienne
#' méthode statistique}, INRA,
#' \emph{Actualités Scientifiques et Agronomiques}, Masson,
#' Paris.
#' @keywords datasets
#' @examples
#' 
#' data(XpineNAX21)
#' str(XpineNAX21)
#' 
NULL
