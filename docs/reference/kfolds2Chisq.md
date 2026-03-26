# Computes Predicted Chisquare for k-fold cross-validated partial least squares regression models.

This function computes Predicted Chisquare for k-fold cross validated
partial least squares regression models.

## Usage

``` r
kfolds2Chisq(pls_kfolds)
```

## Arguments

- pls_kfolds:

  a k-fold cross validated partial least squares regression glm model

## Value

- list:

  Total Predicted Chisquare vs number of components for the first group
  partition

- list():

  ...

- list:

  Total Predicted Chisquare vs number of components for the last group
  partition

## Note

Use
[`cv.plsRglm`](https://fbertran.github.io/plsRglm/reference/cv.plsRglm.md)
to create k-fold cross validated partial least squares regression glm
models.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.md),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md),
[`kfolds2Chisqind`](https://fbertran.github.io/plsRglm/reference/kfolds2Chisqind.md),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.md)
and
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.md)
to extract and transforms results from k-fold cross validation.

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=16,verbose=FALSE)
bbb2 <- cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=5,verbose=FALSE)
kfolds2Chisq(bbb)
#> [[1]]
#> [1] 55.70774 24.52966 20.84377
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] 62.06902 21.17915 19.65666
#> 
rm(list=c("XCornell","yCornell","bbb","bbb2"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsRglm(object=ypine,dataX=Xpine,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
bbb2 <- cv.plsRglm(object=ypine,dataX=Xpine,nt=10,modele="pls-glm-gaussian",K=10,verbose=FALSE)
kfolds2Chisq(bbb)
#> [[1]]
#> [1] 14.18541 12.66755 11.29117 12.16633
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#>  [1] 16.01037 15.04104 12.82465 13.46816 13.22914 14.40349 14.61850 14.95845
#>  [9] 14.86388 14.91391
#> 
                  
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=10,modele="pls",K=10,verbose=FALSE)
kfolds2Press(bbbNA)
#> [[1]]
#> [1] 13.67592 14.91544 14.81330 12.40331 12.06289 12.11392 12.83221 21.45890
#> [9] 26.94021
#> 
kfolds2Chisq(bbbNA)
#> [[1]]
#> [1] 13.67592 14.91544 14.81330 12.40331 12.06289 12.11392 12.83221 21.45890
#> [9] 26.94021
#> 
bbbNA2 <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
bbbNA3 <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian",K=10,
verbose=FALSE)
kfolds2Chisq(bbbNA2)
#> [[1]]
#> [1] 13.20970 31.65945 28.43378 22.55204
#> 
kfolds2Chisq(bbbNA3)
#> [[1]]
#> [1] 14.32418 14.87616 11.89853 12.62791 14.28605 13.89974 21.94249 12.59579
#> [9] 19.89630
#> 
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2","bbbNA3"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",
family="binomial",verbose=FALSE))
#> [[1]]
#> [1]  186.4604  268.9234 1574.2969 6408.2692
#> 
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic",
verbose=FALSE))
#> [[1]]
#> [1]     247.6534     580.9059   24264.5656 5442717.3419
#> 
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-family",
family=binomial(),K=10,verbose=FALSE))
#> [[1]]
#>  [1] 1.986115e+02 4.625618e+02 1.210752e+04 1.416175e+05 5.670987e+05
#>  [6] 9.457270e+05 3.882446e+06 1.077786e+07 7.300684e+07 1.367845e+08
#> 
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-logistic",
K=10,verbose=FALSE))
#> [[1]]
#>  [1]    193.2747    396.3554   4373.4478  23786.8329  46570.8844  63133.7421
#>  [7]  96978.2077 179158.3425 180320.0528 210506.6641
#> 
rm(list=c("Xaze_compl","yaze_compl"))
# }
```
