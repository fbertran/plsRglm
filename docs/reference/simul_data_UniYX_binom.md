# Data generating function for univariate binomial plsR models

This function generates a single univariate binomial response value
\\Y\\ and a vector of explanatory variables \\(X_1,\ldots,X\_{totdim})\\
drawn from a model with a given number of latent components.

## Usage

``` r
simul_data_UniYX_binom(totdim, ncomp, link = "logit", offset = 0)
```

## Arguments

- totdim:

  Number of columns of the X vector (from `ncomp` to hardware limits)

- ncomp:

  Number of latent components in the model (from 2 to 6)

- link:

  Character specification of the link function in the mean model (mu).
  Currently, "`logit`", "`probit`", "`cloglog`", "`cauchit`", "`log`",
  "`loglog`" are supported. Alternatively, an object of class "link-glm"
  can be supplied.

- offset:

  Offset on the linear scale

## Value

- vector:

  \\(Y,X_1,\ldots,X\_{totdim})\\

## Details

This function should be combined with the replicate function to give
rise to a larger dataset. The algorithm used is a modification of a port
of the one described in the article of Li which is a multivariate
generalization of the algorithm of Naes and Martens.

## References

T. Naes, H. Martens, Comparison of prediction methods for multicollinear
data, Commun. Stat., Simul. 14 (1985) 545-576.  
Morris, Elaine B. Martin, Model selection for partial least squares
regression, Chemometrics and Intelligent Laboratory Systems 64 (2002),
79-89,
[doi:10.1016/S0169-7439(02)00051-5](https://doi.org/10.1016/S0169-7439%2802%2900051-5)
.

## See also

[`simul_data_UniYX`](https://fbertran.github.io/plsRglm/reference/simul_data_UniYX.md)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
layout(matrix(1:6,nrow=2))
# logit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4)))[,1])
# probit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="probit")))[,1])
# cloglog link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="cloglog")))[,1])
# cauchit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="cauchit")))[,1])
# loglog link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="loglog")))[,1])
# log link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="log")))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced

layout(1)


layout(matrix(1:6,nrow=2))
# logit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,offset=5)))[,1])
# probit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="probit",offset=5)))[,1])
# cloglog link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="cloglog",offset=5)))[,1])
# cauchit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="cauchit",offset=5)))[,1])
# loglog link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="loglog",offset=5)))[,1])
# log link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="log",offset=5)))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced

layout(1)


layout(matrix(1:6,nrow=2))
# logit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,offset=-5)))[,1])
# probit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="probit",offset=-5)))[,1])
# cloglog link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="cloglog",offset=-5)))[,1])
# cauchit link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="cauchit",offset=-5)))[,1])
# loglog link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="loglog",offset=-5)))[,1])
# log link
hist(t(replicate(100,simul_data_UniYX_binom(4,4,link="log",offset=-5)))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced

layout(1)
# }
```
