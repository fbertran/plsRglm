\name{bootplsglm}
\alias{bootplsglm}
\title{Non-parametric Bootstrap for PLS generalized linear models}
\description{
Provides a wrapper for the bootstrap function \code{boot} from the \code{boot} R package.\cr Implements non-parametric bootstrap for PLS generalized linear models by case resampling.
}
\usage{
bootplsglm(object, typeboot="fmodel_np", R=250, statistic=coefs.plsRglmnp, 
sim="ordinary", stype="i", stabvalue=1e6,...)
}
\arguments{
  \item{object}{ ~~Explain \code{object} here~~ }
  \item{typeboot}{ ~~Explain \code{typeboot} here~~ }
  \item{R}{ ~~Explain \code{R} here~~ }
  \item{statistic}{ ~~Explain \code{statistic} here~~ }
  \item{sim}{ ~~Explain \code{sim} here~~ }
  \item{stype}{ ~~Explain \code{stype} here~~ }
  \item{stabvalue}{ ~~Explain \code{stabvalue} here~~ }
  \item{\dots}{ ~~Explain \code{\dots} here~~ }
}
\details{
  ~~ More details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{\code{\link[boot:boot]{boot}}}
\examples{\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

dataset <- cbind(y=yaze_compl,Xaze_compl)

library(boot)
# Lazraq-Cléroux PLS bootstrap Classic

aze_compl.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, 
nt=3, modele="pls-glm-logistic")
# The same
aze_compl.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, 
nt=3, modele="pls-glm-family",family=binomial)
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), 
sim="ordinary", stype="i", R=250)
# The same
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-family",
family=binomial), sim="ordinary", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)

boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))


plot(aze_compl.boot,index=2)
jack.after.boot(aze_compl.boot, index=2, useJ=TRUE, nt=3)
plot(aze_compl.boot, index=2,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), alpha=c(0.025, 0.975), 
sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), 
sim="ordinary", stype="i", index=1)

# PLS bootstrap balanced

aze_compl.boot <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", stype="i", R=250, 
nt=3, modele="pls-glm-logistic")
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), 
sim="balanced", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))



plot(aze_compl.boot)
jack.after.boot(aze_compl.boot, index=1, useJ=TRUE, nt=3)
plot(aze_compl.boot,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 250, 250), alpha=c(0.025, 0.975), 
sim="balanced", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), 
sim="balanced", stype="i", index=1)


# PLS permutation bootstrap

aze_compl.boot <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", 
R=250, nt=3, modele="pls-glm-logistic")
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), 
sim="permutation", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=6)


boxplots.bootpls(aze_compl.boot)


plot(aze_compl.boot)



data(aze)
Xaze<-aze[,2:34]
yaze<-aze$y

dataset <- cbind(y=yaze,Xaze)

library(boot)
aze.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, 
nt=3, modele="pls-glm-logistic")
aze.boot <- bootplsglm(plsRglm(yaze,Xaze,3,modele="pls-glm-logistic"), sim="ordinary", 
stype="i", R=250)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))



plot(aze.boot)
jack.after.boot(aze.boot, index=1, useJ=TRUE, nt=3)
plot(aze.boot,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), alpha=c(0.025, 
0.975), sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze.tilt.boot <- tilt.bootplsglm(plsRglm(yaze,Xaze,3, modele="pls-glm-logistic"), 
statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", 
stype="i", index=1)




data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
dataset <- cbind(y=ybordeaux,Xbordeaux)
options(contrasts = c("contr.treatment", "contr.poly"))
library(boot)
bordeaux.boot3<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), 
sim="permutation", stype="i", R=250)
#bordeaux.boot3 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", 
stype="i", R=250, nt=1, modele="pls-glm-polr")
boxplots.bootpls(bordeaux.boot3)
boxplots.bootpls(bordeaux.boot3,ranget0=TRUE)

bordeaux.boot2<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), 
sim="permutation", stype="i", R=250, strata=unclass(ybordeaux))
bordeaux.boot2 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", 
stype="i", R=250, nt=1, modele="pls-glm-polr", strata=unclass(ybordeaux))
boxplots.bootpls(bordeaux.boot2)


bordeaux.boot4<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), 
sim="balanced", stype="i", R=250)
#bordeaux.boot4 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", 
stype="i", R=250, nt=1, modele="pls-glm-polr")
#tempboot4 <- bordeaux.boot4
#tempboot4$t0=bordeaux.boot4$t0
#tempboot4$t=bordeaux.boot4$t[!is.na(bordeaux.boot4$t[,1]),]
#tempboot4$R=sum(!is.na(bordeaux.boot4$t[,1]))
#tempboot4$call$R<-sum(!is.na(bordeaux.boot4$t[,1]))
#tempboot4
boxplots.bootpls(bordeaux.boot4)
boot.ci(bordeaux.boot4)
confints.bootpls(bordeaux.boot4)
plots.confints.bootpls(confints.bootpls(bordeaux.boot4))

bordeaux.boot5<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), 
sim="balanced", stype="i", R=250, strata=unclass(ybordeaux))
#bordeaux.boot5 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", 
stype="i", R=250, nt=1, modele="pls-glm-polr", strata=unclass(ybordeaux))
#tempboot5 <- bordeaux.boot5
#tempboot5$t0=bordeaux.boot5$t0
#tempboot5$t=bordeaux.boot5$t[!is.na(bordeaux.boot5$t[,1]),]
#tempboot5$R=sum(!is.na(bordeaux.boot5$t[,1]))
#tempboot5$call$R<-sum(!is.na(bordeaux.boot5$t[,1]))
#tempboot5
boxplots.bootpls(bordeaux.boot5)
boot.ci(bordeaux.boot5)
confints.bootpls(bordeaux.boot5)
plots.confints.bootpls(confints.bootpls(bordeaux.boot5))


#install.packages("chemometrics")
library(chemometrics)
data(hyptis)
hyptis
yhyptis <- factor(hyptis$Group,ordered=TRUE)
Xhyptis <- as.data.frame(hyptis[,c(1:6)])
dataset <- cbind(y=yhyptis,Xhyptis)
options(contrasts = c("contr.treatment", "contr.poly"))
hyptis.boot3<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), 
sim="permutation", stype="i", R=250)
#library(boot)
#hyptis.boot3 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", 
stype="i", R=250, nt=3, modele="pls-glm-polr")
rownames(hyptis.boot3$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene",
"Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot3)
boxplots.bootpls(hyptis.boot3,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot3,ranget0=TRUE)
boxplots.bootpls(hyptis.boot3,ranget0=TRUE,xaxisticks=FALSE)

hyptis.boot2<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), 
sim="permutation", stype="i", R=250, strata=unclass(yhyptis))
#hyptis.boot2 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", 
stype="i", R=250, nt=3, modele="pls-glm-polr", strata=unclass(yhyptis))
rownames(hyptis.boot2$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene",
"Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot2)
boxplots.bootpls(hyptis.boot2,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot2,ranget0=TRUE)
boxplots.bootpls(hyptis.boot2,ranget0=TRUE,xaxisticks=FALSE)


hyptis.boot4<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), 
sim="balanced", stype="i", R=250)
#hyptis.boot4 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", 
stype="i", R=500, nt=3, modele="pls-glm-polr")
#tempboot4 <- hyptis.boot4
#tempboot4$t0=hyptis.boot4$t0
#tempboot4$t=hyptis.boot4$t[!is.na(hyptis.boot4$t[,1]),]
#tempboot4$R=sum(!is.na(hyptis.boot4$t[,1]))
#tempboot4$call$R<-sum(!is.na(hyptis.boot4$t[,1]))
#tempboot4
rownames(hyptis.boot4$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene",
"Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot4)
boxplots.bootpls(hyptis.boot4,xaxisticks=FALSE)
confints.bootpls(hyptis.boot4)
plots.confints.bootpls(confints.bootpls(hyptis.boot4),legendpos = "bottomleft",
xaxisticks=TRUE)
plots.confints.bootpls(confints.bootpls(hyptis.boot4),legendpos = "bottomleft",
xaxisticks=FALSE)


hyptis.boot5 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", 
stype="i", R=250, nt=3, modele="pls-glm-polr", strata=unclass(yhyptis))
#tempboot5 <- hyptis.boot5
#tempboot5$t0=hyptis.boot5$t0
#tempboot5$t=hyptis.boot5$t[!is.na(hyptis.boot5$t[,1]),]
#tempboot5$R=sum(!is.na(hyptis.boot5$t[,1]))
#tempboot5$call$R<-sum(!is.na(hyptis.boot5$t[,1]))
#tempboot5
rownames(hyptis.boot5$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene",
"Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot5)
boxplots.bootpls(hyptis.boot5,xaxisticks=FALSE)
confints.bootpls(hyptis.boot5)
plots.confints.bootpls(confints.bootpls(hyptis.boot5),legendpos = "bottomleft",
xaxisticks=TRUE)
plots.confints.bootpls(confints.bootpls(hyptis.boot5),legendpos = "bottomleft",
xaxisticks=FALSE)
}
}\keyword{models}
