## Test des performances sur le jeux de données bordeaux (Qualité pour glm, température pour le reste)

Unit: milliseconds

|expr|        min|         lq|       mean|     median|        uq|       max|
|----|-----------|-----------|-----------|-----------|----------|----------|
|polr| 122.559078| 125.493234| 150.261750| 147.028672| 157.50174| 465.13368|
|poisson|  46.882535|  48.626915|  56.368827|  49.498429|  64.71539| 107.12152|
|pls|   7.012132|   7.436828|   9.065606|   7.686204|   8.90115|  45.85319|

![](Test_naif_performances_plsrglm.png)


## Test des performances pour différents jeux de données(GLM logistic)

```R
microbenchmark(
+   single = plsRglm(yaze, Xaze,nt=10,modele="pls-glm-logistic",pvals.expli=T,verbose=F),
+   multi = plsRglmPar(yaze, Xaze,nt=10,modele="pls-glm-logistic",pvals.expli=T,verbose=F),
+   times=50)
```

### taille = 1 Aze

Unit: milliseconds

|   expr |      min|        lq|      mean|    median|        uq|       max| neval|
|--------|---------|----------|----------|----------|----------|----------|------|
| single | 755.9307|  786.7887|  853.3339|  865.6824|  900.4021| 1143.570 |    50|
|  multi |1178.6416| 1273.7029| 1367.7406| 1335.6206| 1459.3418| 2002.678 |    50|

![](glm_logistic_single_vs_multicore_1Aze.png)

### taille = 100 Aze

Unit: seconds

|   expr|      min|       lq|     mean|   median|       uq|      max| neval|
|-------|---------|---------|---------|---------|---------|---------|------|
| single| 25.57291| 25.64467| 25.85495| 25.72257| 25.93241| 26.96713|    20|
|  multi| 19.62857| 20.11816| 20.83233| 20.65258| 21.12395| 24.51206|    20|

![](glm_logistic_single_vs_multicore_100Aze.png)
