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
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

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
#> [1] 51.17644 30.35008 19.73907 18.50431 29.77615
#> 
rm(list=c("XCornell","yCornell"))

# \donttest{
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
kfolds2Press(cv.plsR(object=ypine,dataX=Xpine,nt=10,NK=1,verbose=FALSE))
#> [[1]]
#>  [1] 11.97615 12.00141 10.69039 10.30474 10.70328 10.57269 11.61333 11.70210
#>  [9] 11.98101 12.10232
#> 
kfolds2Press(cv.plsR(object=ypine,dataX=Xpine,nt=10,NK=2,verbose=FALSE))
#> [[1]]
#>  [1] 13.31403 13.81107 12.48677 13.62458 17.49581 18.48085 18.01662 17.96863
#>  [9] 17.94488 18.08595
#> 
#> [[2]]
#>  [1] 12.234426 12.271449 12.117853 10.314004  9.485939  9.493277  9.443128
#>  [8]  9.439081 10.075595 10.187946
#> 

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
kfolds2Press(cv.plsR(object=ypine,dataX=XpineNAX21,nt=10,NK=1,verbose=FALSE))
#> [[1]]
#> [1] 13.27636 13.56210 13.74395 12.95239 11.95871 12.45488 14.71568 18.22976
#> [9] 17.77810
#> 
kfolds2Press(cv.plsR(object=ypine,dataX=XpineNAX21,nt=10,NK=2,verbose=FALSE))
#> [[1]]
#> [1] 13.24267 16.33994 15.91731 14.25737 13.78751 14.49856 17.29316 17.10006
#> [9] 16.57428
#> 
#> [[2]]
#> [1] 12.22707 15.25976 17.52732 16.28843 14.28885 12.85276 13.40780 22.87170
#> [9] 42.66878
#> 
rm(list=c("Xpine","XpineNAX21","ypine"))
# }
```
