# Experimental multivariate-response PLS2 models

`plsRmulti()` implements an experimental complete-case linear PLS2 fit
for multivariate numeric responses. It is intentionally separate from
[`plsR`](https://fbertran.github.io/plsRglm/reference/plsR.md) so the
current PLS1 API remains unchanged.

## Usage

``` r
plsRmulti(object, ...)

# Default S3 method
plsRmultiModel(
  object,
  dataX,
  nt = 2,
  limQ2set = 0.0975,
  dataPredictY,
  modele = "pls",
  family = NULL,
  typeVC = "none",
  EstimXNA = FALSE,
  scaleX = TRUE,
  scaleY = NULL,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights,
  sparse = FALSE,
  sparseStop = FALSE,
  naive = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'formula'
plsRmultiModel(
  object,
  data,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  typeVC = "none",
  EstimXNA = FALSE,
  scaleX = TRUE,
  scaleY = NULL,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights,
  subset,
  contrasts = NULL,
  sparse = FALSE,
  sparseStop = FALSE,
  naive = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  For the default method, a numeric multivariate response matrix or data
  frame with at least two columns. For the formula method, a formula of
  the form `cbind(y1, y2, ...) ~ .`.

- ...:

  Not used. Extra arguments are rejected in this experimental release.

- dataX:

  Numeric predictor matrix or data frame.

- nt:

  Number of components to extract.

- limQ2set:

  Kept for interface compatibility. Not supported by `plsRmulti`.

- dataPredictY:

  Kept for interface compatibility. Not supported by `plsRmulti`; fit
  first and then use [`predict`](https://rdrr.io/r/stats/predict.html).

- modele:

  Only `"pls"` is supported.

- family:

  Not supported in this experimental release.

- typeVC:

  Only `"none"` is supported.

- EstimXNA:

  Not supported in this experimental release.

- scaleX:

  Should predictors be scaled?

- scaleY:

  Should responses be scaled? Defaults to `TRUE`.

- pvals.expli:

  Not supported in this experimental release.

- alpha.pvals.expli:

  Not supported in this experimental release.

- MClassed:

  Not supported in this experimental release.

- tol_Xi:

  Tolerance used for degeneracy checks during component extraction.

- weights:

  Not supported in this experimental release.

- sparse:

  Not supported in this experimental release.

- sparseStop:

  Not supported in this experimental release.

- naive:

  Not supported in this experimental release.

- verbose:

  Should informational messages be displayed?

- data:

  An optional data frame for the formula method.

- subset:

  An optional subset for the formula method.

- contrasts:

  Optional contrasts for the formula method.

## Value

An object of class `"plsRmultiModel"` with multivariate analogues of the
linear `plsR` outputs, including the extracted scores `tt`, X loadings
`pp`, response score coefficients `CoeffC`, coefficient matrix `Coeffs`,
intercept vector `CoeffConstante`, scaled response matrix `RepY`, and
fitted response matrices `YChapeau`, `Std.ValsPredictY`, and
`ValsPredictY`.

## Details

This experimental release supports complete-case linear PLS2 fitting,
prediction, repeated k-fold cross-validation via
[`cv.plsRmulti`](https://fbertran.github.io/plsRglm/reference/cv.plsRmulti.md),
and bootstrap resampling via
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md). It
still does not support missing values, weights, sparse extraction,
classification diagnostics, or GLM families.

## See also

[`predict.plsRmultiModel`](https://fbertran.github.io/plsRglm/reference/predict.plsRmultiModel.md),
[`cv.plsRmulti`](https://fbertran.github.io/plsRglm/reference/cv.plsRmulti.md),
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md),
[`plsR`](https://fbertran.github.io/plsRglm/reference/plsR.md)

## Examples

``` r
set.seed(123)
X <- matrix(rnorm(60 * 4), ncol = 4)
Y <- cbind(
  y1 = X[, 1] - 0.5 * X[, 2] + rnorm(60, sd = 0.1),
  y2 = 0.3 * X[, 2] + X[, 3] + rnorm(60, sd = 0.1)
)

fit <- plsRmulti(Y, X, nt = 2, verbose = FALSE)
fit
#> Experimental multivariate PLS2 model
#> Number of responses:
#> [1] 2
#> Number of required components:
#> [1] 2
#> Number of successfully computed components:
#> [1] 2
#> Coefficient matrix:
#>               y1         y2
#> X.1  1.001797150 0.01301818
#> X.2 -0.483516920 0.32994042
#> X.3 -0.029327449 0.94052217
#> X.4 -0.004080383 0.18787469
head(predict(fit))
#>            y1         y2
#> 1 -0.72491996  0.0330750
#> 2  0.05416003 -0.8182979
#> 3  1.75767102 -0.6124818
#> 4  0.59341225 -0.7343838
#> 5  0.61387418  1.3420938
#> 6  1.61054712 -0.5235085
```
