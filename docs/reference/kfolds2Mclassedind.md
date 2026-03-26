# Number of missclassified individuals per group for k-fold cross validated partial least squares regression models.

This function indicates the number of missclassified individuals per
group for k-fold cross validated partial least squares regression
models.

## Usage

``` r
kfolds2Mclassedind(pls_kfolds)
```

## Arguments

- pls_kfolds:

  a k-fold cross validated partial least squares regression model used
  on binary data

## Value

- list:

  Number of missclassified individuals per group vs number of components
  for the first group partition

- list():

  ...

- list:

  Number of missclassified individuals per group vs number of components
  for the last group partition

## Note

Use [`cv.plsR`](https://fbertran.github.io/plsRglm/reference/cv.plsR.md)
or
[`cv.plsRglm`](https://fbertran.github.io/plsRglm/reference/cv.plsRglm.md)
to create k-fold cross validated partial least squares regression models
or generalized linear ones.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.md),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md)
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
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
kfolds2Mclassedind(cv.plsR(object=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=1,verbose=FALSE))
#> [[1]]
#> [[1]][[1]]
#>  [1] 5 6 6 7 7 8 8 8 8 8
#> 
#> [[1]][[2]]
#>  [1] 6 5 5 5 5 5 4 3 4 4
#> 
#> [[1]][[3]]
#>  [1] 4 4 4 5 5 5 5 5 5 5
#> 
#> [[1]][[4]]
#>  [1] 4 5 5 5 5 5 5 5 5 5
#> 
#> [[1]][[5]]
#>  [1] 6 5 4 4 4 4 4 4 4 4
#> 
#> [[1]][[6]]
#>  [1] 5 6 4 4 4 4 4 5 5 5
#> 
#> [[1]][[7]]
#>  [1] 7 7 5 5 5 5 5 5 5 5
#> 
#> [[1]][[8]]
#>  [1] 8 8 8 8 8 8 8 8 8 8
#> 
#> 
kfolds2Mclassedind(cv.plsR(object=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=2,verbose=FALSE))
#> [[1]]
#> [[1]][[1]]
#>  [1] 5 6 3 4 7 7 6 6 6 6
#> 
#> [[1]][[2]]
#>  [1] 6 7 7 7 7 6 6 6 6 6
#> 
#> [[1]][[3]]
#>  [1] 8 8 7 7 7 7 7 7 7 6
#> 
#> [[1]][[4]]
#>  [1] 4 5 4 5 4 4 4 4 4 4
#> 
#> [[1]][[5]]
#>  [1] 5 5 7 6 5 5 6 6 6 6
#> 
#> [[1]][[6]]
#>  [1] 5 3 4 3 3 3 4 4 4 4
#> 
#> [[1]][[7]]
#>  [1] 5 5 5 5 5 5 5 5 6 6
#> 
#> [[1]][[8]]
#>  [1] 6 8 8 8 8 8 7 7 7 7
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>  [1] 7 6 5 5 5 5 5 5 5 5
#> 
#> [[2]][[2]]
#>  [1] 1 2 3 5 5 5 5 5 5 5
#> 
#> [[2]][[3]]
#>  [1] 6 6 7 6 6 7 8 7 7 7
#> 
#> [[2]][[4]]
#>  [1] 5 4 5 5 7 7 7 7 7 7
#> 
#> [[2]][[5]]
#>  [1] 6 8 6 6 6 6 6 6 7 7
#> 
#> [[2]][[6]]
#>  [1] 8 8 7 7 7 7 7 7 7 7
#> 
#> [[2]][[7]]
#>  [1] 5 6 6 6 5 5 5 5 5 5
#> 
#> [[2]][[8]]
#>  [1] 6 7 7 6 7 6 5 5 5 5
#> 
#> 
rm(list=c("Xaze_compl","yaze_compl"))
# }
```
