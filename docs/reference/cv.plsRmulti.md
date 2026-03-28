# Experimental k-fold cross-validation for multivariate-response PLS2

`cv.plsRmulti()` performs repeated k-fold cross-validation for the
experimental complete-case linear
[`plsRmulti`](https://fbertran.github.io/plsRglm/reference/plsRmulti.md)
workflow.

## Usage

``` r
cv.plsRmulti(object, ...)

# Default S3 method
cv.plsRmultiModel(
  object,
  dataX,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  K = 5,
  NK = 1,
  grouplist = NULL,
  random = TRUE,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepfolds = FALSE,
  keepdataY = TRUE,
  keepMclassed = FALSE,
  EstimXNA = FALSE,
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
cv.plsRmultiModel(
  object,
  data = NULL,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  K = 5,
  NK = 1,
  grouplist = NULL,
  random = TRUE,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepfolds = FALSE,
  keepdataY = TRUE,
  keepMclassed = FALSE,
  EstimXNA = FALSE,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights = NULL,
  subset = NULL,
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

  Number of components to extract in each fold fit.

- limQ2set:

  Threshold used by
  [`cvtable`](https://fbertran.github.io/plsRglm/reference/cvtable.md)
  for the aggregated `Q2` criterion.

- modele:

  Only `"pls"` is supported.

- family:

  Not supported in this experimental release.

- K:

  Number of groups for each partition.

- NK:

  Number of repeated partitions.

- grouplist:

  Optional user-supplied partitions.

- random:

  Should the folds be generated randomly?

- scaleX:

  Should predictors be scaled?

- scaleY:

  Should responses be scaled? Defaults to `TRUE`.

- keepcoeffs:

  Should standardized coefficient vectors be stored for each fold fit?

- keepfolds:

  Should training indices be stored for each fold fit?

- keepdataY:

  Kept for interface compatibility. Observed fold responses are stored
  so that summaries can be computed.

- keepMclassed:

  Not supported in this experimental release.

- EstimXNA:

  Not supported in this experimental release.

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

An object of class `"cv.plsRmultiModel"` with repeated fold predictions,
observed fold responses, optional coefficient vectors and fold indices,
and the reference full-data `"plsRmultiModel"` fit used for aggregated
summary metrics.

## Details

Only the linear multivariate-response PLS2 mode is supported here.
Missing values, weights, sparse extraction options, classification
diagnostics, and GLM families remain out of scope for this experimental
API.

## See also

[`plsRmulti`](https://fbertran.github.io/plsRglm/reference/plsRmulti.md),
[`summary.cv.plsRmultiModel`](https://fbertran.github.io/plsRglm/reference/summary.cv.plsRmultiModel.md),
[`cvtable`](https://fbertran.github.io/plsRglm/reference/cvtable.md),
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md)

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
