# Print method for plsRglm models

This function provides a print method for the class `"cv.plsRglmmodel"`

## Usage

``` r
# S3 method for class 'cv.plsRglmmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"cv.plsRglmmodel"`

- ...:

  not used

## Value

`NULL`

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparaison de la régression PLS et de la régression logistique PLS :
application aux données d'allélotypage. *Journal de la Société Française
de Statistique*, 151(2), pages 1-18.
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
print(cv.plsRglm(object=yCornell,dataX=XCornell,nt=10,NK=1,
modele="pls-glm-family",family=gaussian(), verbose=FALSE))
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 5
rm(list=c("XCornell","yCornell","bbb"))
#> Warning: object 'bbb' not found
```
