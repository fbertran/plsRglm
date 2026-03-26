# coef method for plsR models

This function provides a coef method for the class `"plsRglmmodel"`

## Usage

``` r
# S3 method for class 'plsRglmmodel'
coef(object, type = c("scaled", "original"), ...)
```

## Arguments

- object:

  an object of the class `"plsRglmmodel"`

- type:

  if `scaled`, the coefficients of the predictors are given for the
  scaled predictors, if `original` the coefficients are to be used with
  the predictors on their original scale.

- ...:

  not used

## Value

An object of class `coef.plsRglmmodel`.  

- CoeffC:

  Coefficients of the components.

- Std.Coeffs:

  Coefficients of the scaled predictors in the regression function.

- Coeffs:

  Coefficients of the untransformed predictors (on their original
  scale).

## See also

[`coef`](https://rdrr.io/r/stats/coef.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsRglm(yCornell,XCornell,3,modele="pls-glm-family",family=gaussian())
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
class(modpls)
#> [1] "plsRglmmodel"
coef(modpls)
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
coef(modpls,type="scaled")
#> Coefficients of the components
#> Coeff_Comp_Reg 1 Coeff_Comp_Reg 2 Coeff_Comp_Reg 3 
#>        3.1434906        2.8745210        0.9766999 
#> Coefficients of the predictors (scaled scale)
#>                 [,1]
#> Intercept 88.5833333
#> X1        -0.5473211
#> X2        -0.4045974
#> X3        -0.5171213
#> X4        -1.1338145
#> X5         0.1486891
#> X6         4.3120090
#> X7        -0.5663178
rm(list=c("XCornell","yCornell","modpls"))
```
