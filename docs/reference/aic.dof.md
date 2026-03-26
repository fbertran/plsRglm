# Akaike and Bayesian Information Criteria and Generalized minimum description length

This function computes the Akaike and Bayesian Information Criteria and
the Generalized minimum description length.

## Usage

``` r
aic.dof(RSS, n, DoF, sigmahat)

bic.dof(RSS, n, DoF, sigmahat)

gmdl.dof(sigmahat, n, DoF, yhat)
```

## Arguments

- RSS:

  vector of residual sum of squares.

- n:

  number of observations.

- DoF:

  vector of Degrees of Freedom. The length of `DoF` is the same as the
  length of `RSS`.

- sigmahat:

  Estimated model error. The length of `sigmahat` is the same as the
  length of `RSS`.

- yhat:

  vector of squared norm of Yhat. The length of `yhat` is the same as
  the length of `sigmahat`.

## Value

- vector:

  numerical values of the requested AIC, BIC or GMDL.

## Details

The gmdl criterion is defined as
\$\$gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)\$\$
with \$\$S=\hat\sigma^2\$\$

## References

M. Hansen, B. Yu. (2001). Model Selection and Minimum Descripion Length
Principle, *Journal of the American Statistical Association*, 96,
746-774.  
N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of Partial Least
Squares Regression. *Journal of the American Statistical Association*,
106(494), 697-705.  
N. Kraemer, M.L. Braun, Kernelizing PLS, Degrees of Freedom, and
Efficient Model Selection, *Proceedings of the 24th International
Conference on Machine Learning*, Omni Press, (2007) 441-448.

## See also

[`plsR.dof`](https://fbertran.github.io/plsRglm/reference/plsR.dof.md)
for degrees of freedom computation and
[`infcrit.dof`](https://fbertran.github.io/plsRglm/reference/infcrit.dof.md)
for computing information criteria directly from a previously fitted
plsR model.

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
dof.object <- plsR.dof(modpls)
aic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
#>          [,1]     [,2]     [,3]      [,4]      [,5]
#> [1,] 46.07088 4.569969 2.107546 0.8467795 0.8232505
bic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
#>          [,1]     [,2]     [,3]      [,4]      [,5]
#> [1,] 47.78935 4.955816 2.394933 0.9628191 0.9357846
gmdl.dof(dof.object$sigmahat,modpls$nr,dof.object$DoF,dof.object$yhat)
#> [1] 27.59461 21.34020 27.40202 24.40842 24.23105
naive.object <- plsR.dof(modpls,naive=TRUE)
aic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
#>          [,1]     [,2]     [,3]      [,4]      [,5]
#> [1,] 46.07088 4.169957 1.537029 0.7363469 0.8721072
bic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
#>          [,1]    [,2]     [,3]      [,4]      [,5]
#> [1,] 47.78935 4.45882 1.686092 0.8256118 0.9964867
gmdl.dof(naive.object$sigmahat,modpls$nr,naive.object$DoF,naive.object$yhat)
#> [1] 27.59461 18.37545 17.71117 19.01033 24.16510
```
