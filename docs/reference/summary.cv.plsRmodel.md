# Summary method for plsR models

This function provides a summary method for the class `"cv.plsRmodel"`

## Usage

``` r
# S3 method for class 'cv.plsRmodel'
summary(object, ...)
```

## Arguments

- object:

  an object of the class `"cv.plsRmodel"`

- ...:

  further arguments to be passed to or from methods.

## Value

An object of class `"summary.cv.plsRglmmodel"`.

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
summary(cv.plsR(Y~.,data=Cornell,nt=10,K=6, verbose=FALSE), verbose=FALSE)
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y        Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA          NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.8775326  0.0975  0.87753262 57.28983  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8757271  0.0975 -0.01474319 36.26944  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.7730217  0.0975 -0.82644991 20.21260   4.418081 0.9905556
#> Nb_Comp_4 33.76477  0.2356285  0.0975 -2.36759728 14.87832   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -2.7185700  0.0975 -3.86487254 20.96388   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975          NA       NA   3.496074 0.9925265
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000001    0.7633342  0.9711319  1.1359499 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
#> attr(,"computed_nt")
#> [1] 6
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"
```
