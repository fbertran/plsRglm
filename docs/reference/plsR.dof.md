# Computation of the Degrees of Freedom

This function computes the Degrees of Freedom using the Krylov
representation of PLS and other quantities that are used to get
information criteria values. For the time present, it only works with
complete datasets.

## Usage

``` r
# S3 method for class 'dof'
plsR(modplsR, naive = FALSE)
```

## Arguments

- modplsR:

  A plsR model i.e. an object returned by one of the functions `plsR`,
  `plsRmodel.default`, `plsRmodel.formula`, `PLS_lm` or
  `PLS_lm_formula`.

- naive:

  A boolean.

## Value

- DoF:

  Degrees of Freedom

- sigmahat:

  Estimates of dispersion

- Yhat:

  Predicted values

- yhat:

  Square Euclidean norms of the predicted values

- RSS:

  Residual Sums of Squares

## Details

If `naive=FALSE` returns values for estimated degrees of freedom and
error dispersion. If `naive=TRUE` returns returns values for naive
degrees of freedom and error dispersion. The original code from Nicole
Kraemer and Mikio L. Braun was unable to handle models with only one
component.

## References

N. Kraemer, M. Sugiyama. (2011). The Degrees of Freedom of Partial Least
Squares Regression. *Journal of the American Statistical Association*,
106(494), 697-705.  
N. Kraemer, M. Sugiyama, M.L. Braun. (2009). Lanczos Approximations for
the Speedup of Kernel Partial Least Squares Regression, *Proceedings of
the Twelfth International Conference on Artificial Intelligence and
Statistics (AISTATS)*, 272-279.

## See also

[`aic.dof`](https://fbertran.github.io/plsRglm/reference/aic.dof.md) and
[`infcrit.dof`](https://fbertran.github.io/plsRglm/reference/infcrit.dof.md)
for computing information criteria directly from a previously fitted
plsR model.

## Author

Nicole Kraemer, Mikio L. Braun with improvements from Frédéric
Bertrand  
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
plsR.dof(modpls) 
#> $DoF
#> [1] 1.000000 2.740749 5.085967 5.121086 5.103312
#> 
#> $sigmahat
#> [1] 6.5212706 1.8665281 1.1825195 0.7488308 0.7387162
#> 
#> $Yhat
#>        Nt_0     Nt_1     Nt_2     Nt_3     Nt_4
#> 1  88.58333 95.03164 96.49529 97.55864 97.46485
#> 2  88.58333 96.36241 97.51798 97.59151 97.72308
#> 3  88.58333 95.91107 97.56157 97.44534 97.54025
#> 4  88.58333 94.98711 92.14504 91.80793 91.88728
#> 5  88.58333 88.36924 87.98130 85.99452 85.97024
#> 6  88.58333 93.65634 91.12235 91.77507 91.62905
#> 7  88.58333 81.65428 81.33705 81.49670 81.50559
#> 8  88.58333 82.31664 82.49582 82.57542 82.61982
#> 9  88.58333 82.02100 82.28564 82.52399 82.53440
#> 10 88.58333 82.56111 83.12821 83.26027 83.30569
#> 11 88.58333 82.05655 81.23614 81.92924 81.97612
#> 12 88.58333 88.07261 89.69361 89.04136 88.84362
#> 
#> $yhat
#> [1] 94164.08 94596.14 94620.81 94627.46 94627.57
#> 
#> $RSS
#> [1] 467.796667  35.742486  11.066606   4.418081   4.309235
#> 
plsR.dof(modpls,naive=TRUE) 
#> $DoF
#> [1] 1 2 3 4 5
#> 
#> $sigmahat
#> [1] 6.5212706 1.8905683 1.1088836 0.7431421 0.7846050
#> 
#> $Yhat
#>        Nt_0     Nt_1     Nt_2     Nt_3     Nt_4
#> 1  88.58333 95.03164 96.49529 97.55864 97.46485
#> 2  88.58333 96.36241 97.51798 97.59151 97.72308
#> 3  88.58333 95.91107 97.56157 97.44534 97.54025
#> 4  88.58333 94.98711 92.14504 91.80793 91.88728
#> 5  88.58333 88.36924 87.98130 85.99452 85.97024
#> 6  88.58333 93.65634 91.12235 91.77507 91.62905
#> 7  88.58333 81.65428 81.33705 81.49670 81.50559
#> 8  88.58333 82.31664 82.49582 82.57542 82.61982
#> 9  88.58333 82.02100 82.28564 82.52399 82.53440
#> 10 88.58333 82.56111 83.12821 83.26027 83.30569
#> 11 88.58333 82.05655 81.23614 81.92924 81.97612
#> 12 88.58333 88.07261 89.69361 89.04136 88.84362
#> 
#> $yhat
#> [1] 94164.08 94596.14 94620.81 94627.46 94627.57
#> 
#> $RSS
#> [1] 467.796667  35.742486  11.066606   4.418081   4.309235
#> 
```
