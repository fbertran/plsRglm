# Coefficients for permutation bootstrap computations of PLSGLR models

A function passed to `boot` to perform bootstrap.

## Usage

``` r
permcoefs.plsRglmnp(
  dataRepYtt,
  ind,
  nt,
  modele,
  family = NULL,
  maxcoefvalues,
  wwetoile,
  ifbootfail
)
```

## Arguments

- dataRepYtt:

  components' coordinates to bootstrap

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

- maxcoefvalues:

  maximum values allowed for the estimates of the coefficients to
  discard those coming from singular bootstrap samples

- wwetoile:

  values of the Wstar matrix in the original fit

- ifbootfail:

  value to return if the estimation fails on a bootstrap sample

## Value

estimates on a bootstrap sample or `ifbootfail` value if the bootstrap
computation fails.

## Note

\~~some notes\~~

## See also

See also
[`bootplsglm`](https://fbertran.github.io/plsRglm/reference/bootplsglm.md)

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
Cornell.bootYT <- bootplsglm(modplsglm, R=250, statistic=permcoefs.plsRglmnp, verbose=FALSE)
```
