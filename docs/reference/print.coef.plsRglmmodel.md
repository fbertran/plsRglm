# Print method for plsRglm models

This function provides a print method for the class
`"coef.plsRglmmodel"`

## Usage

``` r
# S3 method for class 'coef.plsRglmmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"coef.plsRglmmodel"`

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
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-family",family=gaussian())
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
print(coef(modplsglm))
#> Coefficients of the components
#> Coeff_Comp_Reg 1 Coeff_Comp_Reg 2 Coeff_Comp_Reg 3 
#>        3.1434906        2.8745210        0.9766999 
#> Coefficients of the predictors (original scale)
#>                 [,1]
#> Intercept  87.652763
#> X1         -5.930456
#> X2         -2.069198
#> X3         -9.607722
#> X4         -4.994568
#> X5          2.603934
#> X6         14.721801
#> X7        -20.912671
rm(list=c("XCornell","yCornell","modplsglm"))
```
