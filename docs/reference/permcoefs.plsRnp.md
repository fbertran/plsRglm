# Coefficients computation for permutation bootstrap

A function passed to `boot` to perform bootstrap.

## Usage

``` r
permcoefs.plsRnp(
  dataRepYtt,
  ind,
  nt,
  modele,
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

## See also

See also
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS (Y,X) bootstrap
# statistic=coefs.plsR is the default for (Y,X) resampling of PLSR models.
set.seed(250)
modpls <- plsR(yCornell,XCornell,1)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
Cornell.bootYT <- bootpls(modpls, R=250, typeboot="fmodel_np", sim="permutation",
statistic=permcoefs.plsRnp, verbose=FALSE)
```
