# Print method for plsR models

This function provides a print method for the class `"coef.plsRmodel"`

## Usage

``` r
# S3 method for class 'coef.plsRmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"coef.plsRmodel"`

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
print(coef(modpls))
#> Coefficients of the components
#> Coeff_Comp_Reg 1 Coeff_Comp_Reg 2 Coeff_Comp_Reg 3 
#>        0.4820365        0.2731127        0.1030689 
#> Coefficients of the predictors (original scale)
#>                 [,1]
#> Intercept  92.675989
#> X1         -9.828318
#> X2         -6.960181
#> X3        -16.666239
#> X4         -8.421802
#> X5         -4.388934
#> X6         10.161304
#> X7        -34.528959
rm(list=c("XCornell","yCornell","modpls"))
```
