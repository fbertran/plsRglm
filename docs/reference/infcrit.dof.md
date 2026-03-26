# Information criteria

This function computes information criteria for existing plsR model
using Degrees of Freedom estimation.

## Usage

``` r
infcrit.dof(modplsR, naive = FALSE)
```

## Arguments

- modplsR:

  A plsR model i.e. an object returned by one of the functions `plsR`,
  `plsRmodel.default`, `plsRmodel.formula`, `PLS_lm` or
  `PLS_lm_formula`.

- naive:

  A boolean.

## Value

- matrix:

  AIC, BIC and gmdl values or `NULL`.

## Details

If `naive=FALSE` returns AIC, BIC and gmdl values for estimated and
naive degrees of freedom. If `naive=TRUE` returns `NULL`.

## References

M. Hansen, B. Yu. (2001). Model Selection and Minimum Descripion Length
Principle, *Journal of the American Statistical Association*, 96,
746-774.  
N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of Partial Least
Squares Regression. *Journal of the American Statistical Association*,
106(494), 697-705.  
N. Kraemer, M. Sugiyama, M.L. Braun. (2009). Lanczos Approximations for
the Speedup of Kernel Partial Least Squares Regression, *Proceedings of
the Twelfth International Conference on Artificial Intelligence and
Statistics (AISTATS)*, 272-279.

## See also

[`plsR.dof`](https://fbertran.github.io/plsRglm/reference/plsR.dof.md)
for degrees of freedom computation and `infcrit.dof` for computing
information criteria directly from a previously fitted plsR model.

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
infcrit.dof(modpls)
#>            DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof DoF.naive
#> Nb_Comp_0 1.000000    6.5212706 46.0708838 47.7893514 27.59461         1
#> Nb_Comp_1 2.740749    1.8665281  4.5699686  4.9558156 21.34020         2
#> Nb_Comp_2 5.085967    1.1825195  2.1075461  2.3949331 27.40202         3
#> Nb_Comp_3 5.121086    0.7488308  0.8467795  0.9628191 24.40842         4
#> Nb_Comp_4 5.103312    0.7387162  0.8232505  0.9357846 24.23105         5
#>           sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4      0.7846050  0.8721072  0.9964867   24.16510
```
