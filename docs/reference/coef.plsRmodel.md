# coef method for plsR models

This function provides a coef method for the class `"plsRmodel"`

## Usage

``` r
# S3 method for class 'plsRmodel'
coef(object, type = c("scaled", "original"), ...)
```

## Arguments

- object:

  an object of the class `"plsRmodel"`

- type:

  if `scaled`, the coefficients of the predictors are given for the
  scaled predictors, if `original` the coefficients are to be used with
  the predictors on their original scale.

- ...:

  not used

## Value

An object of class `coef.plsRmodel`.  

- CoeffC:

  Coefficients of the components.

- Std.Coeffs:

  Coefficients of the scaled predictors.

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
coef(modpls)
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
coef(modpls,type="scaled")
#> Coefficients of the components
#> Coeff_Comp_Reg 1 Coeff_Comp_Reg 2 Coeff_Comp_Reg 3 
#>        0.4820365        0.2731127        0.1030689 
#> Coefficients of the predictors (scaled scale)
#>                  [,1]
#> Intercept  0.00000000
#> X1        -0.13909167
#> X2        -0.20869374
#> X3        -0.13755531
#> X4        -0.29316826
#> X5        -0.03843049
#> X6         0.45638984
#> X7        -0.14338442
rm(list=c("XCornell","yCornell","modpls"))
```
