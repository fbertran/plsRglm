# Summary method for experimental multivariate PLS2 CV models

Summarizes repeated k-fold cross-validation results from
[`cv.plsRmulti`](https://fbertran.github.io/plsRglm/reference/cv.plsRmulti.md).

## Usage

``` r
# S3 method for class 'cv.plsRmultiModel'
summary(object, verbose = TRUE, ...)
```

## Arguments

- object:

  An object of class `"cv.plsRmultiModel"`.

- verbose:

  Should progress information be displayed?

- ...:

  Further arguments passed to methods.

## Value

A list of per-partition summary matrices with the same aggregate columns
used by
[`summary.cv.plsRmodel`](https://fbertran.github.io/plsRglm/reference/summary.cv.plsRmodel.md)
for `Q2`, `PRESS`, and `RSS`, plus response-specific `PRESS`, `RSS`,
`Q2`, and `R2` columns.

## Details

The returned object inherits from `"summary.cv.plsRmodel"` so that
[`cvtable`](https://fbertran.github.io/plsRglm/reference/cvtable.md) and
the existing plot method can be reused for the aggregated multivariate
criteria.

## See also

[`cv.plsRmulti`](https://fbertran.github.io/plsRglm/reference/cv.plsRmulti.md),
[`cvtable`](https://fbertran.github.io/plsRglm/reference/cvtable.md)

## Examples

``` r
set.seed(123)
X <- matrix(rnorm(60 * 4), ncol = 4)
Y <- cbind(
  y1 = X[, 1] - 0.5 * X[, 2] + rnorm(60, sd = 0.1),
  y2 = 0.3 * X[, 2] + X[, 3] + rnorm(60, sd = 0.1)
)

cv_fit <- cv.plsRmulti(Y, X, nt = 2, K = 3, NK = 1, verbose = FALSE)
summary(cv_fit, verbose = FALSE)
#> [[1]]
#>           AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0  NA          NA      NA          NA        NA 118.000000        NA
#> Nb_Comp_1  NA -0.08370946  0.0975 -0.08370946 59.406927  54.818131 0.5354396
#> Nb_Comp_2  NA -0.76796595  0.0975 -0.63140215  4.509567   2.764228 0.9765743
#>           AIC.std   PRESS_y1     RSS_y1        Q2_y1     R2_y1  PRESS_y2
#> Nb_Comp_0      NA         NA 59.0000000           NA        NA        NA
#> Nb_Comp_1      NA 31.6161004 31.8142050  0.006226923 0.4607762 27.790827
#> Nb_Comp_2      NA  0.7826158  0.5521728 -0.417338684 0.9906411  3.726952
#>              RSS_y2      Q2_y2     R2_y2
#> Nb_Comp_0 59.000000         NA        NA
#> Nb_Comp_1 23.003926 -0.2080906 0.6101030
#> Nb_Comp_2  2.212055 -0.6848366 0.9625075
#> attr(,"computed_nt")
#> [1] 2
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmultiModel" "summary.cv.plsRmodel"     
```
