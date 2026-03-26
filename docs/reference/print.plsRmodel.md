# Print method for plsR models

This function provides a print method for the class `"plsRmodel"`

## Usage

``` r
# S3 method for class 'plsRmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"plsRmodel"`

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
print(modpls)
#> Number of required components:
#> [1] 3
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                 [,1]
#> Intercept  92.675989
#> X1         -9.828318
#> X2         -6.960181
#> X3        -16.666239
#> X4         -8.421802
#> X5         -4.388934
#> X6         10.161304
#> X7        -34.528959
#> Information criteria and Fit statistics:
#>                AIC      RSS_Y      R2_Y R2_residY RSS_residY    AIC.std
#> Nb_Comp_0 82.01205 467.796667        NA        NA 11.0000000  37.010388
#> Nb_Comp_1 53.15173  35.742486 0.9235940 0.9235940  0.8404663   8.150064
#> Nb_Comp_2 41.08283  11.066606 0.9763431 0.9763431  0.2602256  -3.918831
#> Nb_Comp_3 32.06411   4.418081 0.9905556 0.9905556  0.1038889 -12.937550
#>            DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof DoF.naive
#> Nb_Comp_0 1.000000    6.5212706 46.0708838 47.7893514 27.59461         1
#> Nb_Comp_1 2.740749    1.8665281  4.5699686  4.9558156 21.34020         2
#> Nb_Comp_2 5.085967    1.1825195  2.1075461  2.3949331 27.40202         3
#> Nb_Comp_3 5.121086    0.7488308  0.8467795  0.9628191 24.40842         4
#>           sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3      0.7431421  0.7363469  0.8256118   19.01033
rm(list=c("XCornell","yCornell","modpls"))
```
