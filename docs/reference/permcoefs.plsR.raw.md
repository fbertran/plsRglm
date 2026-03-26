# Raw coefficients for permutation bootstrap computations of PLSR models

A function passed to `boot` to perform bootstrap.

## Usage

``` r
permcoefs.plsR.raw(
  dataset,
  ind,
  nt,
  modele,
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
  [plsR](https://fbertran.github.io/plsRglm/reference/plsR.md)

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
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md).

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
set.seed(250)
modpls <- permcoefs.plsR.raw(Cornell[,-8],1:nrow(Cornell),nt=3,
maxcoefvalues=1e5,ifbootfail=rep(0,3),verbose=FALSE)
```
