# Plot method for table of summary of cross validated plsRglm models

This function provides a table method for the class
`"summary.cv.plsRglmmodel"`

## Usage

``` r
# S3 method for class 'table.summary.cv.plsRglmmodel'
plot(x, type = c("CVMC", "CVQ2Chi2", "CVPreChi2"), ...)
```

## Arguments

- x:

  an object of the class `"table.summary.cv.plsRglmmodel"`

- type:

  the type of cross validation criterion to plot.

- ...:

  further arguments to be passed to or from methods.

## Value

`NULL`

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`summary`](https://rdrr.io/r/base/summary.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,
modele="pls-glm-family",family=gaussian(), verbose=FALSE)
plot(cvtable(summary(bbb,verbose=FALSE)),type="CVQ2Chi2")
#> 
#> CV Q2Chi2 criterion:
#> 0 1 2 
#> 0 0 1 
#> 
#> CV PreChi2 criterion:
#> 1 2 
#> 0 1 

rm(list=c("bbb"))

# \donttest{
data(Cornell)
plot(cvtable(summary(cv.plsRglm(Y~.,data=Cornell,nt=10,NK=100,
modele="pls-glm-family",family=gaussian(), verbose=FALSE),
verbose=FALSE)),type="CVQ2Chi2")
#> 
#> CV Q2Chi2 criterion:
#>  0  1  2 
#>  0 30 70 
#> 
#> CV PreChi2 criterion:
#>  1  2  3  4  5 
#>  0 20 55 22  3 

# }
```
