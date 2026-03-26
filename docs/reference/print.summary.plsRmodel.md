# Print method for summaries of plsR models

This function provides a print method for the class
`"summary.plsRmodel"`

## Usage

``` r
# S3 method for class 'summary.plsRmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"summary.plsRmodel"`

- ...:

  not used

## Value

- language:

  call of the model

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparaison de la régression PLS et de la régression logistique PLS :
application aux données d'allélotypage. *Journal de la Société Française
de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`print`](https://rdrr.io/r/base/print.html) and
[`summary`](https://rdrr.io/r/base/summary.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
class(modpls)
#> [1] "plsRglmmodel"
print(summary(modpls))
#> Call:
#> plsRglmmodel.default(dataX = XCornell, nt = 3, modele = "pls", 
#>     dataY = yCornell)
rm(list=c("XCornell","yCornell","modpls"))
```
