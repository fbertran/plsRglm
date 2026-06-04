# Computes individual Predicted Chisquare for k-fold cross validated partial least squares regression models.

This function computes individual Predicted Chisquare for k-fold cross
validated partial least squares regression models.

## Usage

``` r
kfolds2Chisqind(pls_kfolds)
```

## Arguments

- pls_kfolds:

  a k-fold cross validated partial least squares regression glm model

## Value

- list:

  Individual PChisq vs number of components for the first group
  partition

- list():

  ...

- list:

  Individual PChisq vs number of components for the last group partition

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
[`kfolds2Chisq`](https://fbertran.github.io/plsRglm/reference/kfolds2Chisq.md),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.md)
and
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.md)
to extract and transforms results from k-fold cross-validation.

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
kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 
#> [[1]][[3]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 
#> [[1]][[7]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 
#> [[1]][[8]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 
#> [[1]][[9]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 
#> [[1]][[10]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 
#> [[1]][[11]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 
#> [[1]][[12]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 
#> 
kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA
#> 
#> 
rm(list=c("XCornell","yCornell","bbb","bbb2"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsRglm(object=ypine,dataX=Xpine,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
bbb2 <- cv.plsRglm(object=ypine,dataX=Xpine,nt=10,modele="pls-glm-gaussian",K=10,verbose=FALSE)
kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA
#> 
#> 
kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
                  
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=10,modele="pls",K=10,verbose=FALSE)
kfolds2Pressind(bbbNA)
#> [[1]]
#> [[1]][[1]]
#> [1] 1.1240258 1.3442425 0.9058616 0.6864513 0.8559083 0.8182619 0.9949588
#> [8] 1.1207491 4.9255270
#> 
#> [[1]][[2]]
#> [1]  1.9847940  1.5194752  1.0427757  0.5824120 15.4231786  8.6156547  3.3106528
#> [8]  0.5240196  0.6727968
#> 
#> [[1]][[3]]
#> [1] 0.5368438 0.6478287 0.9479271 0.4078544 0.3486778 0.3690438 0.3368690
#> [8] 0.2056592 0.3311661
#> 
#> [[1]][[4]]
#> [1] 4.618568 4.230738 3.758332 3.998854 3.872621 3.883942 3.963446 3.969494
#> [9] 3.856248
#> 
#> [[1]][[5]]
#> [1] 0.4877758 0.6947223 1.1021406 0.4197469 0.3579208 0.3987580 0.3964675
#> [8] 0.3839523 0.3101484
#> 
#> [[1]][[6]]
#> [1] 0.9626907 0.5602840 0.7669114 0.3614088 0.3769860 0.4006152 0.4589790
#> [8] 0.4835146 0.4363950
#> 
#> [[1]][[7]]
#> [1] 1.0253261 0.6714788 0.5522761 0.5964822 0.6313382 0.8195699 0.6601231
#> [8] 0.7192523 0.9289249
#> 
#> [[1]][[8]]
#> [1] 0.0884500 0.8203197 0.6541690 0.9896632 0.6982577 0.6266222 0.6496231
#> [8] 0.6208857 0.7021905
#> 
#> [[1]][[9]]
#> [1] 1.417490 1.414474 1.315603 1.296645 1.261910 1.275216 1.260917 1.271807
#> [9] 1.527300
#> 
#> [[1]][[10]]
#> [1] 0.7348488 0.6609862 0.5138980 1.1468310 1.0621610 0.9921056 1.1127299
#> [8] 1.2428777 1.2177092
#> 
#> 
kfolds2Chisqind(bbbNA)
#> [[1]]
#> [[1]][[1]]
#> [1] 1.1240258 1.3442425 0.9058616 0.6864513 0.8559083 0.8182619 0.9949588
#> [8] 1.1207491 4.9255270
#> 
#> [[1]][[2]]
#> [1]  1.9847940  1.5194752  1.0427757  0.5824120 15.4231786  8.6156547  3.3106528
#> [8]  0.5240196  0.6727968
#> 
#> [[1]][[3]]
#> [1] 0.5368438 0.6478287 0.9479271 0.4078544 0.3486778 0.3690438 0.3368690
#> [8] 0.2056592 0.3311661
#> 
#> [[1]][[4]]
#> [1] 4.618568 4.230738 3.758332 3.998854 3.872621 3.883942 3.963446 3.969494
#> [9] 3.856248
#> 
#> [[1]][[5]]
#> [1] 0.4877758 0.6947223 1.1021406 0.4197469 0.3579208 0.3987580 0.3964675
#> [8] 0.3839523 0.3101484
#> 
#> [[1]][[6]]
#> [1] 0.9626907 0.5602840 0.7669114 0.3614088 0.3769860 0.4006152 0.4589790
#> [8] 0.4835146 0.4363950
#> 
#> [[1]][[7]]
#> [1] 1.0253261 0.6714788 0.5522761 0.5964822 0.6313382 0.8195699 0.6601231
#> [8] 0.7192523 0.9289249
#> 
#> [[1]][[8]]
#> [1] 0.0884500 0.8203197 0.6541690 0.9896632 0.6982577 0.6266222 0.6496231
#> [8] 0.6208857 0.7021905
#> 
#> [[1]][[9]]
#> [1] 1.417490 1.414474 1.315603 1.296645 1.261910 1.275216 1.260917 1.271807
#> [9] 1.527300
#> 
#> [[1]][[10]]
#> [1] 0.7348488 0.6609862 0.5138980 1.1468310 1.0621610 0.9921056 1.1127299
#> [8] 1.2428777 1.2177092
#> 
#> 
bbbNA2 <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=4,modele="pls-glm-gaussian",verbose=FALSE)
bbbNA3 <- cv.plsRglm(object=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian",
K=10,verbose=FALSE)
kfolds2Chisqind(bbbNA2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA
#> 
#> 
kfolds2Chisqind(bbbNA3)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> 
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2","bbbNA3"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
kfolds2Chisqind(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",
family=binomial(),verbose=FALSE))
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA
#> 
#> 
kfolds2Chisqind(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic",
verbose=FALSE))
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA
#> 
#> 
kfolds2Chisqind(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-family",
family=binomial(),K=10,verbose=FALSE))
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisqind(cv.plsRglm(object=yaze_compl,dataX=Xaze_compl,nt=10,
modele="pls-glm-logistic",K=10,verbose=FALSE))
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
rm(list=c("Xaze_compl","yaze_compl"))
# }
```
