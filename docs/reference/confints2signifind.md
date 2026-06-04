# Detect significant predictors from bootstrap confidence intervals

This function converts the matrix returned by
[`confints.bootpls`](https://fbertran.github.io/plsRglm/reference/confints.bootpls.md)
into a logical indicator of significance. A predictor is marked as
significant when the selected confidence interval lies entirely on one
side of the reference value.

## Usage

``` r
confints2signifind(ic_bootobject, typeIC, threshold = 0)
```

## Arguments

- ic_bootobject:

  a matrix returned by
  [`confints.bootpls`](https://fbertran.github.io/plsRglm/reference/confints.bootpls.md).

- typeIC:

  type of confidence interval to use. Defaults to `"BCa"` when BCa
  limits are available and to `"Percentile"` otherwise.

- threshold:

  reference value to test against. Defaults to `0`.

## Value

Named logical vector. `TRUE` means that the selected confidence interval
excludes `threshold`; `FALSE` means that it contains `threshold`; `NA`
means that the interval could not be evaluated.

## See also

[`confints.bootpls`](https://fbertran.github.io/plsRglm/reference/confints.bootpls.md)
and
[`signpred`](https://fbertran.github.io/plsRglm/reference/signpred.md).

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data(Cornell)
set.seed(250)
modpls <- plsR(Y ~ ., data = Cornell, 3)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
Cornell.bootYT <- bootpls(modpls, typeboot = "fmodel_np", R = 250, verbose = FALSE)
temp.ci <- confints.bootpls(Cornell.bootYT, indices = 2:8, typeBCa = FALSE)
confints2signifind(temp.ci)
#>    X1    X2    X3    X4    X5    X6    X7 
#>  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
# }
```
