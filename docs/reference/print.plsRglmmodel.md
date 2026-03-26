# Print method for plsRglm models

This function provides a print method for the class `"plsRglmmodel"`

## Usage

``` r
# S3 method for class 'plsRglmmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"plsRglmmodel"`

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
print(modplsglm)
#> Number of required components:
#> [1] 3
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                 [,1]
#> Intercept  87.652763
#> X1         -5.930456
#> X2         -2.069198
#> X3         -9.607722
#> X4         -4.994568
#> X5          2.603934
#> X6         14.721801
#> X7        -20.912671
#> Information criteria and Fit statistics:
#>                AIC      BIC Chi2_Pearson_Y      RSS_Y      R2_Y R2_residY
#> Nb_Comp_0 82.01205 82.98186     467.796667 467.796667        NA        NA
#> Nb_Comp_1 53.15173 54.60645      35.742486  35.742486 0.9235940 0.9235940
#> Nb_Comp_2 31.46903 33.40866       4.966831   4.966831 0.9893825 0.9893825
#> Nb_Comp_3 31.54404 33.96857       4.230693   4.230693 0.9909561 0.9909561
#>           RSS_residY
#> Nb_Comp_0 467.796667
#> Nb_Comp_1  35.742486
#> Nb_Comp_2   4.966831
#> Nb_Comp_3   4.230693
#> Model with all the required components:
#> 
#> Call:  glm(formula = YwotNA ~ ., family = family, data = tttrain)
#> 
#> Coefficients:
#> (Intercept)         tt.1         tt.2         tt.3  
#>     88.5833       3.1435       2.8745       0.9767  
#> 
#> Degrees of Freedom: 11 Total (i.e. Null);  8 Residual
#> Null Deviance:       467.8 
#> Residual Deviance: 4.231     AIC: 31.54
rm(list=c("XCornell","yCornell","modplsglm"))
```
