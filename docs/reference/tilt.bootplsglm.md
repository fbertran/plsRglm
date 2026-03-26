# Non-parametric tilted bootstrap for PLS generalized linear regression models

Provides a wrapper for the bootstrap function `tilt.boot` from the
`boot` R package.  
Implements non-parametric tilted bootstrap for PLS generalized linear
regression models by case resampling : the `tilt.boot` function will run
an initial bootstrap with equal resampling probabilities (if required)
and will use the output of the initial run to find resampling
probabilities which put the value of the statistic at required values.
It then runs an importance resampling bootstrap using the calculated
probabilities as the resampling distribution.

## Usage

``` r
tilt.bootplsglm(
  object,
  typeboot = "fmodel_np",
  statistic = coefs.plsRglm,
  R = c(499, 250, 250),
  alpha = c(0.025, 0.975),
  sim = "ordinary",
  stype = "i",
  index = 1,
  stabvalue = 1e+06,
  ...
)
```

## Arguments

- object:

  An object of class `plsRbetamodel` to bootstrap

- typeboot:

  The type of bootstrap. Either (Y,X) boostrap (`typeboot="plsmodel"`)
  or (Y,T) bootstrap (`typeboot="fmodel_np"`). Defaults to (Y,T)
  resampling.

- statistic:

  A function which when applied to data returns a vector containing the
  statistic(s) of interest. `statistic` must take at least two
  arguments. The first argument passed will always be the original data.
  The second will be a vector of indices, frequencies or weights which
  define the bootstrap sample. Further, if predictions are required,
  then a third argument is required which would be a vector of the
  random indices used to generate the bootstrap predictions. Any further
  arguments can be passed to statistic through the `...` argument.

- R:

  The number of bootstrap replicates. Usually this will be a single
  positive integer. For importance resampling, some resamples may use
  one set of weights and others use a different set of weights. In this
  case `R` would be a vector of integers where each component gives the
  number of resamples from each of the rows of weights.

- alpha:

  The alpha level to which tilting is required. This parameter is
  ignored if `R[1]` is 0 or if `theta` is supplied, otherwise it is used
  to find the values of `theta` as quantiles of the initial uniform
  bootstrap. In this case `R[1]` should be large enough that
  `min(c(alpha, 1-alpha))*R[1] > 5`, if this is not the case then a
  warning is generated to the effect that the `theta` are extreme values
  and so the tilted output may be unreliable.

- sim:

  A character string indicating the type of simulation required.
  Possible values are `"ordinary"` (the default), `"balanced"`,
  `"permutation"`, or `"antithetic"`.

- stype:

  A character string indicating what the second argument of `statistic`
  represents. Possible values of stype are `"i"` (indices - the
  default), `"f"` (frequencies), or `"w"` (weights).

- index:

  The index of the statistic of interest in the output from `statistic`.
  By default the first element of the output of `statistic` is used.

- stabvalue:

  Upper bound for the absolute value of the coefficients.

- ...:

  ny further arguments can be passed to `statistic`.

## Value

An object of class "boot".

## See also

[`tilt.boot`](https://rdrr.io/pkg/boot/man/boot.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

dataset <- cbind(y=yaze_compl,Xaze_compl)

# Lazraq-Cleroux PLS bootstrap Classic

aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic", family=NULL), statistic=coefs.plsRglm, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
boxplots.bootpls(aze_compl.tilt.boot,1:2)


aze_compl.tilt.boot2 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic"), statistic=coefs.plsRglm, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
boxplots.bootpls(aze_compl.tilt.boot2,1:2)


aze_compl.tilt.boot3 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-family", family=binomial), statistic=coefs.plsRglm, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
boxplots.bootpls(aze_compl.tilt.boot3,1:2)



# PLS bootstrap balanced

aze_compl.tilt.boot4 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, 
modele="pls-glm-logistic"), statistic=coefs.plsRglm, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
boxplots.bootpls(aze_compl.tilt.boot4,1:2)

# }
```
