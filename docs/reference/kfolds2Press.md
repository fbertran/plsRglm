# Computes PRESS for k-fold cross validated partial least squares regression models.

This function computes PRESS for k-fold cross validated partial least
squares regression models.

## Usage

``` r
kfolds2Press(pls_kfolds)
```

## Arguments

- pls_kfolds:

  a k-fold cross validated partial least squares regression model

## Value

- list:

  Press vs number of components for the first group partition

- list():

  ...

- list:

  Press vs number of components for the last group partition

## Note

Use [`cv.plsR`](https://fbertran.github.io/plsRglm/reference/cv.plsR.md)
to create k-fold cross validated partial least squares regression
models.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/47/>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md),
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
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
kfolds2Press(cv.plsR(object=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),
nt=6,K=12,NK=1,verbose=FALSE))
#> [[1]]
#> [1] 55.70774 41.43274 20.27397 21.24240 24.51801
#> 
kfolds2Press(cv.plsR(object=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),
nt=6,K=6,NK=1,verbose=FALSE))
#> [[1]]
#> [1] 55.56262 40.00977 19.05723 19.30291 22.60499
#> 
rm(list=c("XCornell","yCornell"))

# \donttest{
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
kfolds2Press(cv.plsR(object=ypine,dataX=Xpine,nt=10,NK=1,verbose=FALSE))
#> [[1]]
#>  [1] 12.17362 13.32722 11.89863 10.76002 10.78439 10.87965 11.09345 13.43135
#>  [9] 14.32411 14.54937
#> 
kfolds2Press(cv.plsR(object=ypine,dataX=Xpine,nt=10,NK=2,verbose=FALSE))
#> [[1]]
#>  [1] 12.728449 11.554432  9.794321  9.631880  9.640802  9.904351 10.462499
#>  [8] 10.344699 10.385719 10.604705
#> 
#> [[2]]
#>  [1] 13.82029 15.47818 15.22997 13.84117 14.29024 14.38570 15.81367 16.75103
#>  [9] 17.70783 18.81556
#> 

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
kfolds2Press(cv.plsR(object=ypine,dataX=XpineNAX21,nt=10,NK=1,verbose=FALSE))
#> [[1]]
#> [1] 13.94508 14.72641 14.90373 15.47681 15.55398 15.75101 16.21811 18.29196
#> [9] 19.58054
#> 
kfolds2Press(cv.plsR(object=ypine,dataX=XpineNAX21,nt=10,NK=2,verbose=FALSE))
#> [[1]]
#> [1] 14.44350 14.29197 15.49908 16.12835 24.10455 18.42922 19.68966 25.39844
#> [9] 24.77465
#> 
#> [[2]]
#> [1] 15.55066 17.75194 16.72691 14.86151 15.45137 15.88554 16.09298 22.61966
#> [9] 22.42832
#> 
rm(list=c("Xpine","XpineNAX21","ypine"))
# }
```
