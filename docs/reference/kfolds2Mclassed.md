# Number of missclassified individuals for k-fold cross validated partial least squares regression models.

This function indicates the total number of missclassified individuals
for k-fold cross validated partial least squares regression models.

## Usage

``` r
kfolds2Mclassed(pls_kfolds)
```

## Arguments

- pls_kfolds:

  a k-fold cross validated partial least squares regression model used
  on binary data

## Value

- list:

  Total number of missclassified individuals vs number of components for
  the first group partition

- list():

  ...

- list:

  Total number of missclassified individuals vs number of components for
  the last group partition

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
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md)
and
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.md)
to extract and transforms results from k-fold cross validation.

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
kfolds2Mclassed(cv.plsR(object=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=1,verbose=FALSE))
#> [[1]]
#>  [1] 45 46 42 45 44 45 45 44 44 43
#> 
kfolds2Mclassed(cv.plsR(object=yaze_compl,dataX=Xaze_compl,nt=10,K=8,NK=2,verbose=FALSE))
#> [[1]]
#>  [1] 41 45 46 47 46 46 44 44 44 44
#> 
#> [[2]]
#>  [1] 44 43 44 44 46 46 46 47 48 49
#> 
rm(list=c("Xaze_compl","yaze_compl"))
# }
```
