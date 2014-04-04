\name{XpineNAX21}
\alias{XpineNAX21}
\docType{data}
\title{Missing data analysis for the pine dataset}
\description{
Le jeu de données comporte 33 observations sur la 
processionnaire du pin (variete de chenille). 
L'etude a ete realisee en 1973.
Les experimentateurs souhaitent connaitre l'influence
de caracteristiques forestieres sur le developpement
de la processionnaire du pin. 

L'unite qui represente l'observation est une parcelle
forestiere de 10 hectares d'un seul tenant. Cette parcelle
est consideree comme homogene par rapport aux variables
etudiees. Les valeurs de ces variables ont ete obtenues
comme des moyennes des valeurs mesurees sur des placettes
de 5 ares.

La valeur de x2 pour la première observation a été volontairement supprimée de la matrice des variables explicatives.
}
\usage{data(XpineNAX21)}
\format{
  A data frame with 33 observations on the following 10 variables.
  \describe{
    \item{\code{x1}}{altitude (en metre)}
    \item{\code{x2}}{pente (en degre)}
    \item{\code{x3}}{nombre de pins dans une placette de 5 ares}
    \item{\code{x4}}{hauteur de l'arbre echantillonne au centre de la placette (en metre)}
    \item{\code{x5}}{diametre de cet arbre (en metre)}
    \item{\code{x6}}{note de densite du peuplement}
    \item{\code{x7}}{orientation de la placette (1=sud, 2=nord)}
    \item{\code{x8}}{hauteur des arbres dominants (en metre)}
    \item{\code{x9}}{nombre de strates de vegetation}
    \item{\code{x10}}{melange du peuplement (1=pas melange, 2=melange)}
  }
}
\details{
La variable a expliquer est la variable x11 (log décimal du nombre de nids de
processionnaires par arbre d'une placette). C est pour des raisons
partiellement statistiques que  l'on analyse le log du nb de nids que 
le nb de nids lui-même.
}
\source{
R. Tomassone et al.
\enc{Actualités}{Actualit\'es} scientifiques et agronomiques
Masson, Paris (1992)
Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992), 
\dQuote{La \enc{régression}{r\'egression}, nouveaux regards sur une ancienne \enc{méthode}{m\'ethode} statistique}, 
INRA, \emph{\enc{Actualités}{Actualit\'es} Scientifiques et Agronomiques}, Masson, Paris, 188 p.
}
%\references{}
\examples{
data(XpineNAX21)
## maybe str(XpineNAX21) ; plot(XpineNAX21) ...
}
\keyword{datasets}
