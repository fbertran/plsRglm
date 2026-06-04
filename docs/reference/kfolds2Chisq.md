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
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/47/>

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
#> [1] NA NA NA
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA
#> 
rm(list=c("XCornell","yCornell","bbb","bbb2"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsRglm(object=ypine,dataX=Xpine,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
bbb2 <- cv.plsRglm(object=ypine,dataX=Xpine,nt=10,modele="pls-glm-gaussian",K=10,verbose=FALSE)
kfolds2Chisq(bbb)
#> [[1]]
#> [1] NA NA NA NA
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
                  
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=10,modele="pls",K=10,verbose=FALSE)
kfolds2Press(bbbNA)
#> [[1]]
#> [1] 14.57505 14.63146 13.91856 13.41291 15.14871 15.27842 47.43333 59.65936
#> [9] 15.49328
#> 
kfolds2Chisq(bbbNA)
#> [[1]]
#> [1] 14.57505 14.63146 13.91856 13.41291 15.14871 15.27842 47.43333 59.65936
#> [9] 15.49328
#> 
bbbNA2 <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
bbbNA3 <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian",K=10,
verbose=FALSE)
kfolds2Chisq(bbbNA2)
#> [[1]]
#> [1] NA NA NA NA
#> 
kfolds2Chisq(bbbNA3)
#> [[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2","bbbNA3"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",
family="binomial",verbose=FALSE))
#> [[1]]
#> [1] NA NA NA NA
#> 
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic",
verbose=FALSE))
#> [[1]]
#> [1] NA NA NA NA
#> 
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-family",
family=binomial(),K=10,verbose=FALSE))
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
kfolds2Chisq(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-logistic",
K=10,verbose=FALSE))
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
rm(list=c("Xaze_compl","yaze_compl"))
# }
```
