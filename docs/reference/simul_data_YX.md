# Data generating function for multivariate plsR models

This function generates a single multivariate response value
\\\boldsymbol{Y}\\ and a vector of explinatory variables
\\(X_1,\ldots,X\_{totdim})\\ drawn from a model with a given number of
latent components.

## Usage

``` r
simul_data_YX(totdim, ncomp)
```

## Arguments

- totdim:

  Number of column of the X vector (from `ncomp` to hardware limits)

- ncomp:

  Number of latent components in the model (from 2 to 6)

## Value

- vector:

  \\(Y_1,\ldots,Y_r,X_1,\ldots,X\_{totdim})\\

## Details

This function should be combined with the replicate function to give
rise to a larger dataset. The algorithm used is a port of the one
described in the article of Li which is a multivariate generalization of
the algorithm of Naes and Martens.

## Note

The value of \\r\\ depends on the value of `ncomp` :

|         |       |
|---------|-------|
| `ncomp` | \\r\\ |
| 2       | 3     |
| 3       | 3     |
| 4       | 4     |

## References

T. Naes, H. Martens, Comparison of prediction methods for multicollinear
data, Commun. Stat., Simul. 14 (1985) 545-576.  
Morris, Elaine B. Martin, Model selection for partial least squares
regression, Chemometrics and Intelligent Laboratory Systems 64 (2002)
79-89,
[doi:10.1016/S0169-7439(02)00051-5](https://doi.org/10.1016/S0169-7439%2802%2900051-5)
.

## See also

[`simul_data_complete`](https://fbertran.github.io/plsRglm/reference/simul_data_complete.md)
for highlighting the simulations parameters

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
simul_data_YX(20,6)                          
#>        Y 1        Y 2        Y 3        Y 4         X1         X2         X3 
#> -4.0719461 -4.0370330 -4.0459331 -4.0552833 -2.0681875 -2.6831433 -2.6892432 
#>         X4         X5         X6         X7         X8         X9        X10 
#> -1.2244089 -0.5400524 -1.1837768 -1.1693010  0.2837632 -0.2179875 -0.6139394 
#>        X11        X12        X13        X14        X15        X16        X17 
#> -1.8399970 -1.9919613 -2.4186501 -1.3666176 -0.2062973 -0.6271765 -1.8297457 
#>        X18        X19        X20 
#> -1.9994416 -2.4101170 -1.3788371 

# \donttest{
if(require(plsdepot)){
dimX <- 6
Astar <- 2
(dataAstar2 <- t(replicate(50,simul_data_YX(dimX,Astar))))
library(plsdepot)
resAstar2 <- plsreg2(dataAstar2[,4:9],dataAstar2[,1:3],comps=5)
resAstar2$Q2
resAstar2$Q2[,4]>0.0975

dimX <- 6
Astar <- 3
(dataAstar3 <- t(replicate(50,simul_data_YX(dimX,Astar))))
library(plsdepot)
resAstar3 <- plsreg2(dataAstar3[,4:9],dataAstar3[,1:3],comps=5)
resAstar3$Q2
resAstar3$Q2[,4]>0.0975

dimX <- 6
Astar <- 4
(dataAstar4 <- t(replicate(50,simul_data_YX(dimX,Astar))))
library(plsdepot)
resAstar4 <- plsreg2(dataAstar4[,5:10],dataAstar4[,1:4],comps=5)
resAstar4$Q2
resAstar4$Q2[,5]>0.0975

rm(list=c("dimX","Astar"))
}
#> Loading required package: plsdepot
#> 
#> Attaching package: ‘plsdepot’
#> The following object is masked from ‘package:chemometrics’:
#> 
#>     nipals
# }
```
