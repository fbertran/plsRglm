# Bootstrap confidence intervals

This function is a wrapper for
[`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html) to derive
bootstrap-based confidence intervals from a `"boot"` object.

## Usage

``` r
confints.bootpls(bootobject, indices = NULL, typeBCa = TRUE)
```

## Arguments

- bootobject:

  an object of class `"boot"`

- indices:

  the indices of the predictor for which CIs should be calculated.
  Defaults to `NULL`: all the predictors will be used.

- typeBCa:

  shall BCa bootstrap based CI derived ? Defaults to `TRUE`. This is a
  safety option since sometimes computing BCa bootstrap based CI fails
  whereas the other types of CI can still be derived.

## Value

Matrix with the limits of bootstrap based CI for all (defaults) or only
the selected predictors (`indices` option). The limits are given in that
order: Normal Lower then Upper Limit, Basic Lower then Upper Limit,
Percentile Lower then Upper Limit, BCa Lower then Upper Limit.

## See also

See also
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md) and
[`bootplsglm`](https://fbertran.github.io/plsRglm/reference/bootplsglm.md).

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data(Cornell)

#Lazraq-Cleroux PLS (Y,X) bootstrap
set.seed(250)
modpls <- plsR(Y~.,data=Cornell,3)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
Cornell.bootYX <- bootpls(modpls, R=250, verbose=FALSE)
confints.bootpls(Cornell.bootYX,2:8)
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#>                                                                         
#> X1 -0.2305299 -0.03654653 -0.2146155 -0.01243502 -0.2657483 -0.063567788
#> X2 -0.3824730 -0.12056633 -0.4240731 -0.16474400 -0.2526435  0.006685662
#> X3 -0.2262325 -0.03807142 -0.2115428 -0.01464437 -0.2604663 -0.063567788
#> X4 -0.4336032 -0.19861671 -0.4793055 -0.22524165 -0.3610949 -0.107030999
#> X5 -0.2895056  0.13307318 -0.3083408  0.07915147 -0.1560125  0.231479782
#> X6  0.3197348  0.65767612  0.3256605  0.67125328  0.2415264  0.587119147
#> X7 -0.2387634 -0.03963758 -0.2590735 -0.03271142 -0.2540574 -0.027695351
#>                          
#> X1 -0.2867282 -0.07494113
#> X2 -0.2795110 -0.11744873
#> X3 -0.2795040 -0.07955903
#> X4 -0.4109452 -0.17018880
#> X5 -0.1803183  0.17569760
#> X6  0.3172633  0.64752609
#> X7 -0.2222602  0.03146667
#> attr(,"typeBCa")
#> [1] TRUE
confints.bootpls(Cornell.bootYX,2:8,typeBCa=FALSE)
#>                                                                         
#> X1 -0.2305299 -0.03654653 -0.2146155 -0.01243502 -0.2657483 -0.063567788
#> X2 -0.3824730 -0.12056633 -0.4240731 -0.16474400 -0.2526435  0.006685662
#> X3 -0.2262325 -0.03807142 -0.2115428 -0.01464437 -0.2604663 -0.063567788
#> X4 -0.4336032 -0.19861671 -0.4793055 -0.22524165 -0.3610949 -0.107030999
#> X5 -0.2895056  0.13307318 -0.3083408  0.07915147 -0.1560125  0.231479782
#> X6  0.3197348  0.65767612  0.3256605  0.67125328  0.2415264  0.587119147
#> X7 -0.2387634 -0.03963758 -0.2590735 -0.03271142 -0.2540574 -0.027695351
#> attr(,"typeBCa")
#> [1] FALSE
# }
```
