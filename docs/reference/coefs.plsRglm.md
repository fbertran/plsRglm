# Coefficients for bootstrap computations of PLSGLR models

A function passed to `boot` to perform bootstrap.

## Usage

``` r
coefs.plsRglm(
  dataset,
  ind,
  nt,
  modele,
  family = NULL,
  fit_backend = "stats",
  maxcoefvalues,
  ifbootfail,
  verbose
)
```

## Arguments

- dataset:

  dataset to resample

- ind:

  indices for resampling

- nt:

  number of components to use

- modele:

  type of modele to use, see
  [plsRglm](https://fbertran.github.io/plsRglm/reference/plsRglm.md)

- family:

  glm family to use, see
  [plsRglm](https://fbertran.github.io/plsRglm/reference/plsRglm.md)

- fit_backend:

  backend used for repeated non-ordinal score-space GLM fits. Use
  `"stats"` or `"fastglm"`.

- maxcoefvalues:

  maximum values allowed for the estimates of the coefficients to
  discard those coming from singular bootstrap samples

- ifbootfail:

  value to return if the estimation fails on a bootstrap sample

- verbose:

  should info messages be displayed ?

## Value

estimates on a bootstrap sample or `ifbootfail` value if the bootstrap
computation fails.

## See also

See also
[`bootplsglm`](https://fbertran.github.io/plsRglm/reference/bootplsglm.md).

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)

# (Y,X) bootstrap of a PLSGLR model
# statistic=coefs.plsRglm is the default for (Y,X) bootstrap of a PLSGLR models.
set.seed(250)
modplsglm <- plsRglm(Y~.,data=Cornell,1,modele="pls-glm-family",family=gaussian)
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
Cornell.bootYX <- bootplsglm(modplsglm, R=250, typeboot="plsmodel", 
statistic=coefs.plsRglm, verbose=FALSE)
```
