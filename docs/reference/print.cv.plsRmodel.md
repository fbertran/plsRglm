# Print method for plsR models

This function provides a print method for the class `"cv.plsRmodel"`

## Usage

``` r
# S3 method for class 'cv.plsRmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"cv.plsRmodel"`

- ...:

  not used

## Value

`NULL`

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`print`](https://rdrr.io/r/base/print.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
print(cv.plsR(object=yCornell,dataX=XCornell,nt=10,K=6, verbose=FALSE))
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 6
rm(list=c("XCornell","yCornell","bbb"))
#> Warning: object 'bbb' not found
```
