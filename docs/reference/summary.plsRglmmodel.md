# Summary method for plsRglm models

This function provides a summary method for the class `"plsRglmmodel"`

## Usage

``` r
# S3 method for class 'plsRglmmodel'
summary(object, ...)
```

## Arguments

- object:

  an object of the class `"plsRglmmodel"`

- ...:

  further arguments to be passed to or from methods.

## Value

- call :

  function call of plsRglmmodel

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
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
class(modplsglm)
#> [1] "plsRglmmodel"
summary(modplsglm)
#> Call:
#> plsRglmmodel.default(dataX = XCornell, nt = 3, modele = "pls-glm-gaussian", 
#>     dataY = yCornell)
rm(list=c("XCornell","yCornell","modplsglm"))
```
