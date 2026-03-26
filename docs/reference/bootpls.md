# Non-parametric Bootstrap for PLS models

Provides a wrapper for the bootstrap function `boot` from the `boot` R
package.  
Implements non-parametric bootstraps for PLS Regression models by either
(Y,X) or (Y,T) resampling.

## Usage

``` r
bootpls(
  object,
  typeboot = "plsmodel",
  R = 250,
  statistic = NULL,
  sim = "ordinary",
  stype = "i",
  stabvalue = 1e+06,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object of class `plsRmodel` to bootstrap

- typeboot:

  The type of bootstrap. Either (Y,X) boostrap (`typeboot="plsmodel"`)
  or (Y,T) bootstrap (`typeboot="fmodel_np"`). Defaults to (Y,X)
  resampling.

- R:

  The number of bootstrap replicates. Usually this will be a single
  positive integer. For importance resampling, some resamples may use
  one set of weights and others use a different set of weights. In this
  case `R` would be a vector of integers where each component gives the
  number of resamples from each of the rows of weights.

- statistic:

  A function which when applied to data returns a vector containing the
  statistic(s) of interest. `statistic` must take at least two
  arguments. The first argument passed will always be the original data.
  The second will be a vector of indices, frequencies or weights which
  define the bootstrap sample. Further, if predictions are required,
  then a third argument is required which would be a vector of the
  random indices used to generate the bootstrap predictions. Any further
  arguments can be passed to statistic through the `...` argument.

- sim:

  A character string indicating the type of simulation required.
  Possible values are `"ordinary"` (the default), `"balanced"`,
  `"permutation"`, or `"antithetic"`.

- stype:

  A character string indicating what the second argument of `statistic`
  represents. Possible values of stype are `"i"` (indices - the
  default), `"f"` (frequencies), or `"w"` (weights).

- stabvalue:

  A value to hard threshold bootstrap estimates computed from atypical
  resamplings. Especially useful for Generalized Linear Models.

- verbose:

  should info messages be displayed ?

- ...:

  Other named arguments for `statistic` which are passed unchanged each
  time it is called. Any such arguments to `statistic` should follow the
  arguments which `statistic` is required to have for the simulation.
  Beware of partial matching to arguments of `boot` listed above.

## Value

An object of class `"boot"`. See the Value part of the help of the
function [`boot`](https://rdrr.io/pkg/boot/man/boot.html).

## Details

More details on bootstrap techniques are available in the help of the
[`boot`](https://rdrr.io/pkg/boot/man/boot.html) function.

## References

A. Lazraq, R. Cleroux, and J.-P. Gauchi. (2003). Selecting both latent
and explanatory variables in the PLS1 regression model. *Chemometrics
and Intelligent Laboratory Systems*, 66(2):117-126.  
P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS generalised
linear regression. *Computational Statistics & Data Analysis*,
48(1):17-46.  
A. C. Davison and D. V. Hinkley. (1997). *Bootstrap Methods and Their
Applications*. Cambridge University Press, Cambridge.

## See also

[`boot`](https://rdrr.io/pkg/boot/man/boot.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS ordinary bootstrap
set.seed(250)
modpls <- plsR(yCornell,XCornell,3)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 

#(Y,X) resampling
Cornell.bootYX <- bootpls(modpls, R=250, verbose=FALSE)

#(Y,T) resampling
Cornell.bootYT <- bootpls(modpls, typeboot="fmodel_np", R=250, verbose=FALSE)

# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.bootYX,indices=2:8)

# Confidence intervals plotting
confints.bootpls(Cornell.bootYX,indices=2:8)
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#>                                                                         
#> X1 -0.2305299 -0.03654653 -0.2146155 -0.01243502 -0.2657483 -0.063567788
#> X2 -0.3824730 -0.12056633 -0.4240731 -0.16474400 -0.2526435  0.006685662
#> X3 -0.2262325 -0.03807142 -0.2115428 -0.01464437 -0.2604663 -0.063567788
#> X4 -0.4336032 -0.19861671 -0.4793055 -0.22524165 -0.3610949 -0.107030999
#> X5 -0.2895056  0.13307318 -0.3083408  0.07915147 -0.1560125  0.231479782
#> X6  0.3197348  0.65767612  0.3256605  0.67125328  0.2415264  0.587119147
#> X7 -0.2387634 -0.03963758 -0.2590735 -0.03271142 -0.2540574 -0.027695351
#>                          
#> X1 -0.2867282 -0.07494113
#> X2 -0.2795110 -0.11744873
#> X3 -0.2795040 -0.07955903
#> X4 -0.4109452 -0.17018880
#> X5 -0.1803183  0.17569760
#> X6  0.3172633  0.64752609
#> X7 -0.2222602  0.03146667
#> attr(,"typeBCa")
#> [1] TRUE
plots.confints.bootpls(confints.bootpls(Cornell.bootYX,indices=2:8))
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.bootYX$t[,-1])~factor(rep(1:7,rep(250,7))), 
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.bootYX$t0[-1],col="red",pch=19)



# \donttest{
library(boot)
boot.ci(Cornell.bootYX, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
#> BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#> Based on 250 bootstrap replicates
#> 
#> CALL : 
#> boot.ci(boot.out = Cornell.bootYX, conf = c(0.9, 0.95), type = c("norm", 
#>     "basic", "perc", "bca"), index = 2)
#> 
#> Intervals : 
#> Level      Normal              Basic         
#> 90%   (-0.2149, -0.0521 )   (-0.1971, -0.0401 )   
#> 95%   (-0.2305, -0.0365 )   (-0.2146, -0.0124 )  
#> 
#> Level     Percentile            BCa          
#> 90%   (-0.2380, -0.0811 )   (-0.2573, -0.0857 )   
#> 95%   (-0.2657, -0.0636 )   (-0.2867, -0.0749 )  
#> Calculations and Intervals on Original Scale
#> Some basic intervals may be unstable
#> Some percentile intervals may be unstable
#> Some BCa intervals may be unstable
plot(Cornell.bootYX,index=2)

jack.after.boot(Cornell.bootYX, index=2, useJ=TRUE, nt=3)

plot(Cornell.bootYX,index=2,jack=TRUE)


car::dataEllipse(Cornell.bootYX$t[,2], Cornell.bootYX$t[,3], cex=.3, 
levels=c(.5, .95, .99), robust=TRUE)

rm(Cornell.bootYX)


# PLS balanced bootstrap

set.seed(225)
Cornell.bootYX <- bootpls(modpls, sim="balanced", R=250, verbose=FALSE)
boot.array(Cornell.bootYX, indices=TRUE)
#>        [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#>   [1,]    5    4    9    2    4   11    5    3    1     5     9     6
#>   [2,]    1   10    4    1    8    3    7    2   11     8     5    11
#>   [3,]    2    6    1    8    4    1   10    5    1    10     5     7
#>   [4,]    7    9    2   12    5    7    6    9    7     7     6     8
#>   [5,]    5   11    3    8    4    4    6    2    7     1     5     5
#>   [6,]    2    6    4   10    9   11    7    8    5    12    12     4
#>   [7,]    6    6   11    2    7    7    9    2    4     8     2    12
#>   [8,]    5    9    9   11   10   10    7    2    3    10     2     6
#>   [9,]    6    2    3    8    2    9    8   10    6     9    11     7
#>  [10,]    5    8    4    4    2   12    4    1   10    12     8    12
#>  [11,]    2    1    4    6   11   10    9   10    7     6     4     6
#>  [12,]    1    1    1    1   10    2    9    1   12     7    10     1
#>  [13,]    1    2   11    9   10    4    1   10    1     7    12     3
#>  [14,]   11    3    1    7    5   12    6    6    9     8     1     7
#>  [15,]    3   11    8    7    1    9    6    1    4     1     4     5
#>  [16,]   10   12   12    8    9    6    3    5    1     2     8    12
#>  [17,]    1    8    7    8    1    3   10    3    7     9     8     9
#>  [18,]    5    1    1    1    5    4    4   11    8     7     3     9
#>  [19,]    2    5    3    2   12    7    4    7   10    11     4    10
#>  [20,]    2    4    2   11   10   10   10    5    7     7     4    12
#>  [21,]    7   10    8    1   12    4   11    7    9    12     2     7
#>  [22,]    5    5   11    6    8    6   12    9    7     9     4     5
#>  [23,]    8   12    6   12    9    8    1    5    1     2    11     4
#>  [24,]    2    2    3    6   10    1    9   11   11     8     2     2
#>  [25,]    4    7   11    7    7    7    1    6    1     5     3    11
#>  [26,]    4   11    8    8    1    2    6    9   10    11     6    10
#>  [27,]    3    4   12    8   11    8    1   10    9     1     3    12
#>  [28,]    6    8    6    6    4    7   10    9    3     7    10     8
#>  [29,]    5    1    2   12   10    9    4    9    8     6    11     8
#>  [30,]   11   12    4    8    8   10    3   10    8     2     2     7
#>  [31,]    6   11    3    7    5   11    3   12    6     5     4    10
#>  [32,]    8    7    8    4    1    7    3    7    1    11     4     6
#>  [33,]    2   11    5    3    2    8   11    9    2     2     6     9
#>  [34,]    6    1    6    4    8    6    6    3    8     2     3     5
#>  [35,]    6    4   10    5   11    7   10   11    4     5     7     4
#>  [36,]    4    8    1   10    4    4    6    7    6     7     2    10
#>  [37,]    9    1    6    4   10    6    8   11   12     2     4     7
#>  [38,]   11    9    6    7    6    3    3   11    7     6     3    11
#>  [39,]    7    4    1    7    3    3    9   12    2     3    10     3
#>  [40,]    6    2    5    1    9    8   10   12    2     3     5     1
#>  [41,]    8   12    1    9    9    8    6    6   12    12     2     9
#>  [42,]    1    3    7   11   10   10   11   11    9     1     2    11
#>  [43,]   11    2   10    9    6    9    1    9    7    10     1     5
#>  [44,]   11    2   12    7   11    1   12    6    9    10     6     5
#>  [45,]    5    4    9    4    9   10    3    6    9    12     9    12
#>  [46,]    8    2    1   11    6    8   10   11   11    12     9     4
#>  [47,]    6    9    5    1    1    6    9    7    1    10     8     5
#>  [48,]    6    5    3    3   12   10   11   12   11     4    11    11
#>  [49,]    1    1    6    1    6   10    5    6   12     6     4     4
#>  [50,]    2    2    1   12   11    2   12    7    5     5     3     9
#>  [51,]    8    8    9    5   10    3    6   12    3     1     5     6
#>  [52,]    4   12    9    9    5   10    9    6    8     8     5    11
#>  [53,]    3    5    1    4    7   12    8    9   12    10    10     6
#>  [54,]    3    5   11    3   11   12    7   11   10     4     9     2
#>  [55,]    5    7    7    6    2    1    9   11    2     9     7    11
#>  [56,]    8    4    2    7    8    1    6    6   11     3     8     2
#>  [57,]   12    5    9    2    6   10    1   12   11     7     5    12
#>  [58,]    1   12    6   12    7    5    2    6    9     1     6    12
#>  [59,]    7    7    6    2    5    3    9    8    5    10     5     5
#>  [60,]    7    5    2    5    4    7    2    8    4     8     6     2
#>  [61,]    8    8    8   12   10    1   10    8   12    11     1     5
#>  [62,]    1    9    9    4    3    4    7    1    1    11     7     6
#>  [63,]    4    7    5    3    6    7    8    8   10    11    10     8
#>  [64,]    2   10   11    6    6    9   11    6   11     3     8     8
#>  [65,]   12    5    6   10   10   12    1    9    6     5     6     7
#>  [66,]   11    9    2   12    5    2    9   11    5     2     6     7
#>  [67,]   10    1    8   12    6   12    9   10   12     7     4     2
#>  [68,]    7    5    8    3    6   12   12    8    3     6    12     1
#>  [69,]    4    4    3    5    4    2    5   12    8     4     1    10
#>  [70,]   10   11   11    8    8    6    1    9    7     8    11    11
#>  [71,]    3    3    1    1   11    1    1    2    4     6     2     8
#>  [72,]   12   12    8    8    8   12    2   10    3    11     5    11
#>  [73,]    7    9    6    5    1   11    8    5    2    11    12     3
#>  [74,]    6    9    1   10    9    3    1    3    1     4     5     4
#>  [75,]    5    5    8    2    8    2    8    4   11     5     3     4
#>  [76,]    8    2    4    9    8   12   10    6    3     4    11     9
#>  [77,]    5    7    9   10    8    9    4    7   12    12    12     8
#>  [78,]    6    9    9    3    7    3    4   12    4     7     3    11
#>  [79,]   12    2   12    3    5   12    8   12    3     8     8    12
#>  [80,]   12    4    7    6    8    5    7    3    8     9    11     3
#>  [81,]    4    5    6    8   10   12   12    7    6    11     1     4
#>  [82,]    8    3   10   11    5    5    4   12   10     2     8     4
#>  [83,]    9    2    2    9   12    9   12   11    5    11     4     8
#>  [84,]    4    7    9    9    1    2   12    7    9     3    11    12
#>  [85,]    9    3   12    2   11    3    3    4    4     7     8     2
#>  [86,]    6    2    2    7   12    3    8   10    2     5    10    10
#>  [87,]    7    5    6    4    8    4    8    8   12     5     1     7
#>  [88,]    5    4   11    4    2    4    4    9    3     1     8    12
#>  [89,]    1    7    5   10    1   11   12    5    8    10     3     3
#>  [90,]    9    2    6    8    4    8    1    1    9     6     2     3
#>  [91,]   12    5    4   12    7    4    8    3    3    11     5     1
#>  [92,]    5   12   10    4    3   11    9   11    3    12     9     1
#>  [93,]    2    9    4    2    7    8    5    1    6    11     4     5
#>  [94,]   10    1    2    7    7    2    9    2    5    12     6    11
#>  [95,]    3    1    9    4    6   11    2   12   12     4     6     4
#>  [96,]    8    9    3    1    1   11    9    3   12    12    11     5
#>  [97,]    9    1    1    5   11    2   11   11    6    11     7     1
#>  [98,]    8   11    4    1    3    5    9    1   10     8     8     2
#>  [99,]    8    6    3    4    7    1    6    6    7     7     2     9
#> [100,]   10    9   10    2   11    7    1    5    9     5     4     6
#> [101,]    5    6    2    3    6    8   10    1   11     4     2     5
#> [102,]    7    6    1    9    1    1   10    8    9    12    11    11
#> [103,]    2    1    7    3   10    5    3    6    8    10     8     6
#> [104,]    1    9    5    4    4   11   10    2    5     2     2     2
#> [105,]    2    1    2    8    7   12    7   10    9     3     9     5
#> [106,]    1    1   12    8    3    6    8    6    3     3     9    10
#> [107,]    3    5   10   10   11    5    8   10   11     1    12     7
#> [108,]    4   12    7   10    1   10   11    1   10    12     3     8
#> [109,]   11   10    3    7   12    4   12   11    6     7     4    10
#> [110,]   11    4    5    2   12    8    7    7   12     8     3     2
#> [111,]    5   11    3   12    9    3    5    6    6    11    10     4
#> [112,]   11    6    6    5    4    1    6   11   11    12    10     9
#> [113,]    9   12    3    9    9    8    3   10   12     5    11    10
#> [114,]   12    8    6    3    2    4    1    1   12     7     7    12
#> [115,]    9    4    3    9   10    7    4    5    2     9     9     8
#> [116,]    4    7    3    4    6    5    4    6   10     5     3    12
#> [117,]    1    8    2   10   11   10   11    1   12    12     1     2
#> [118,]    6    7    3    7    7    8    9   11   12     9    12     5
#> [119,]    3    9    6    5    8    8    2    8   10     4     1     2
#> [120,]    5    9   11    2   10   10    2   10   12     4    11     6
#> [121,]    9    6    9    4    8    2   12   11   12     3     2     4
#> [122,]    9   11    7   11    2    4   12    4    2     4     9     9
#> [123,]    3    1    4    5    9    7   10    3    7     9     7    10
#> [124,]    5    6    9    3    3    4   11   11    8     1     6     6
#> [125,]    6    8    9    7    2    4    7   10    3     5     2     7
#> [126,]    3    7    2    2   11    7   12    6   10    10     7     1
#> [127,]   10    5    8    3   10    1   11    5   10     4     3    11
#> [128,]    9    6    3    4    2    2   10    2    1     2    11     1
#> [129,]    4    9    1    4    3    9   10    3   12     9     2     6
#> [130,]   11   10    9   10    3    6    8   11   12     6     3     5
#> [131,]    4    3   11    9    7   10    5    8    1    12     9     4
#> [132,]   11    1    9    6   10   10    6    5    2     2    11     4
#> [133,]   10    7    9    9   10    6    6    3    8     3     5     6
#> [134,]    3    8    6   10   12    3    4   11    7    10     4    11
#> [135,]    4    1    6    5    2    8   10    5    9    11     2    11
#> [136,]   11    4   11    2   11    8   10    7    7     9     2    12
#> [137,]    5    2   10    4   12    6    7    5    1     1     9     4
#> [138,]    7    7    1    8    3    4    3    2    6     5     7     1
#> [139,]   12    2    8    3    2    2    7   11   11     4     6    10
#> [140,]    3    3    7    6    8    9    4    6    8     4     7     1
#> [141,]    4    2    3   12    9   10   10   10    5     6    12    10
#> [142,]    4    9    8    8    3    3    4   10   10     1     7     7
#> [143,]   12    3   10    9    2    8    8    6   10     7     9    10
#> [144,]    7    3    4    5   11    2   12   12    7    12     6    10
#> [145,]    1    3    6    9   12   12   10   12    8    11     7     4
#> [146,]    7    3   11    3   12    4   10    4   10     2     5     6
#> [147,]   12    3    8    3    5   11    2    6    7     5     3    10
#> [148,]    3   12    6    4    7   12    6   11    4     2     5     3
#> [149,]    8    5   12    2   10    2    7    4    2     8    11     2
#> [150,]   12    6    7    6   12    4   12    9    3    12     4     3
#> [151,]   11    5    5    4    7   10    2   12    9    10     6     2
#> [152,]    9    3    8    2    7    7    1    6    9     3     8     3
#> [153,]    4    2    7    1    5    3    6    6    7     8    11     4
#> [154,]    1    5    7    5    8   12    4    6    6     9     6     9
#> [155,]    4    8    4    1    8    2    6    4    6    10     9     6
#> [156,]    4   11    3   10    3    2    5    2   12    11     1     5
#> [157,]   12    7    1   11    2    6    3    8    3     8    11     7
#> [158,]    2    2    3   11   11   11    3    6    7     9    11     6
#> [159,]    8   11    2   10    1   11   12   11    3     4     2     5
#> [160,]    4   10    1    5   10    1    6   11    1     6     5     2
#> [161,]   11    9    3   11    6   10    2   11    6     7     9     8
#> [162,]   10   11    1    7    2    6   12    8    6     3    10    11
#> [163,]   11   10    7    1    2    3    4    9    8     6     7     7
#> [164,]    9    2   10    8    1    8    3   10    5     4     9     5
#> [165,]    5    1    4    6   10    7    5    6    2     3     2     2
#> [166,]    8   12   10    3    7    7    3    3   11    12     9    12
#> [167,]    4    9    4   11    8    9    4    6    8     4     3    10
#> [168,]    9    5    3    5    2   10   12    7   12    10     4    11
#> [169,]    2    3    2   10    5    8    7    5   12     7     2    10
#> [170,]   10    6   10    8    6   12    6    3    9     6     6     3
#> [171,]    6    5    6    7    7   12    1    9    4     4    12     3
#> [172,]    5   12    5   12    1    8    5    7    7     1     8     8
#> [173,]    9    7    2    2    1    9    8    5    3     5     5     2
#> [174,]   12    2   10    4   10    5   10    7    5    10     1     6
#> [175,]    1    5    1   11   11   12   10    1    3    12     5     1
#> [176,]    4    8   12    4    9    7    4   12    8     4     1     9
#> [177,]    1   10    2    4   11    4    9    4    2     5     9     8
#> [178,]   10    6    6    2    4    6    5    2    1    10     9     8
#> [179,]    9    8    7    7    8    3    1    5   11     9    11     5
#> [180,]    1    9   11    6    2   10    1   11   11     6    12     3
#> [181,]    9    9   10    6    9    6    6    7    2    12     3     5
#> [182,]   12    2    2    9    7    3    7    4   11    10     3     8
#> [183,]    8    2    9    8   12    1    8    9   10     2    10     5
#> [184,]   11    5   12   10   12    9    1   10    2    12    12     5
#> [185,]   10    1   11    4    8    9    5    6    1     8    10     1
#> [186,]    4   12    5    9   12    1    5   12    9     8     4     9
#> [187,]    3    9    9    9    6    5    7    5    8     4     6    10
#> [188,]    1    5    1    1    2    6   11    8   10     1     5     6
#> [189,]    4    6    5    6    4   12    9    8    2     1    12     9
#> [190,]    3    5    2    4   11   10   12    6   10    11     3    12
#> [191,]    7    9    5    4    8   12    3    4   12     3     8    12
#> [192,]    5    3    7    5    2   11    8    9   11     6     6     7
#> [193,]    3   10   10    8    9    5   10    3    6     9     7     4
#> [194,]    8    7    4    2    7    3    9    3    4     9     2     2
#> [195,]    7   12    8    5    4   10    7    1   10     7     1     7
#> [196,]    4   10    3   12   10   12    7    7    9     8     6     2
#> [197,]    9    2    1    7   11    5    2   10    6     5    12     8
#> [198,]    6    2   11   12    2    8   11    1    9     1     8    11
#> [199,]    4    3    1    3   12    2   11    7    6     4     5     7
#> [200,]    9    2    5    5    9    7    2    2    7     3     5     4
#> [201,]    2    1   10   10    8    5    1    2    6    11     9     4
#> [202,]    5    1   12   11    1    8    5    3    4     9     9     2
#> [203,]    6   12    2    6    8    4    3   12    2    12     1     3
#> [204,]   12    8    3    1   10    4    3    3    8     3     3    10
#> [205,]    1    6   12    5   10    2    7    6    3     9     2     2
#> [206,]    7    8    5    7   12    3    4    6    7     3    11    12
#> [207,]    1   10   12    3   12    3   12    4    3    10     9     6
#> [208,]    2   12    9   10   11    7    8   10    7     9     4     7
#> [209,]    6   11    4    7    9   11   11    1    8    11     8     7
#> [210,]    8    1    9    3    7    7    4    4    3     1     8     1
#> [211,]    5    3    1    9   12    6   10    4    5    12     3    11
#> [212,]    2    2    9    8    7    7   10    9    1     5     4     7
#> [213,]   10    1    4    5    1    9    8    6    9     8     3    11
#> [214,]    4    2    6    3   11    3    5    9    4     1     8    12
#> [215,]    4   10    5    4    9    3    8   12    5     1     6     4
#> [216,]    2    3    5    3    1   11    7    5    1    12     2     9
#> [217,]    3   12   11    7   12    9    6   11   10     6     4     8
#> [218,]    8    6    3    2    9    5    8    9    9     2     7     3
#> [219,]    1    7    6   11    3    4    1   11   12    10     2     9
#> [220,]   11   11    3    8   10    9    8    9    9    10     8     6
#> [221,]    7    8    1    7    4    8   10    5    7     2     3     4
#> [222,]   12    2    6    8   11    3    2    2    3     4    10     3
#> [223,]   10    1    4   11    4   12   10    2    9    11    10     3
#> [224,]    3    8    4    5    8    4    1    8   11     9    12     3
#> [225,]    1   12    7    5   12    5    1    8    2    11     7    11
#> [226,]    9    1   11   12   11    1    3   10    1     4     8     8
#> [227,]    1   11    5    4    1   10    5   12    3    10     5    12
#> [228,]   12    1   10   12    9    5   11    4    2     6    11     5
#> [229,]   10    6    2    1   12   11    8   11   11     9     2     5
#> [230,]    5    4    2   11    1   10   11    7    6     3     2     9
#> [231,]    1    1    8    2    6    7    1   12   12    10     4     3
#> [232,]    3    8   10    3    8   10    3    5   10     7     6     6
#> [233,]    7    5    6    2    2    2   12    6    3     4     2     1
#> [234,]    4    3   10   11    8    7    7    6    5     4     4     1
#> [235,]    6    4    2    6    3   11    6    5    4     5    11     7
#> [236,]    5    4    5   11    2    9    7    5   12     5     7    12
#> [237,]    9    4    3   12    3   12    6    1   10     3    11     1
#> [238,]    7    3    5   11    8    5   11    1    5     7     1    10
#> [239,]    1    7    6    3    7    2    7    1    5     2     4     5
#> [240,]   10   11    1    2   12    8   12    4    7     5     1     3
#> [241,]    8    5    9   11   11    1    6    3    4    11     2     2
#> [242,]    5   11   10   10   12   10   12    4    1     9     4     1
#> [243,]   12    3    7    5    5   10    3    6    8    12     1    10
#> [244,]    5    1    6    2   10    9    8    5    9    12     5     9
#> [245,]   10    4    3    6   12   11   11   12    1    12     2     6
#> [246,]    1    3    2   11    7    9    1    5    5     8    12     5
#> [247,]   10    4    8    9   12    7   10   10    1     9     8     2
#> [248,]    5   12    1    7    5   12   10    8    1     7     8     7
#> [249,]    9    7    3    8    8    5    9    7   10     7     1    11
#> [250,]    7    7    6   11   11    2    9    7    9     3     7    10

# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.bootYX,indices=2:8)

# Confidence intervals plotting
confints.bootpls(Cornell.bootYX,indices=2:8)
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#>                                                                         
#> X1 -0.2192145 -0.03502591 -0.2119416 -0.02617619 -0.2520071 -0.066241725
#> X2 -0.3828012 -0.11801756 -0.4193349 -0.15286851 -0.2645190  0.001947385
#> X3 -0.2154695 -0.03582192 -0.2088689 -0.02820433 -0.2469063 -0.066241725
#> X4 -0.4379789 -0.21238307 -0.4660525 -0.22425891 -0.3620776 -0.120283989
#> X5 -0.2967356  0.13546366 -0.2921161  0.09312255 -0.1699835  0.215255066
#> X6  0.3276202  0.67729268  0.3397828  0.68421162  0.2285681  0.572996919
#> X7 -0.2394465 -0.03754892 -0.2166742 -0.03032279 -0.2564461 -0.070094620
#>                          
#> X1 -0.2459524 -0.04831319
#> X2 -0.2918411 -0.11466891
#> X3 -0.2394376 -0.04850420
#> X4 -0.4645601 -0.20487969
#> X5 -0.2196585  0.18011275
#> X6  0.3325666  0.69760910
#> X7 -0.2561586 -0.06959096
#> attr(,"typeBCa")
#> [1] TRUE
plots.confints.bootpls(confints.bootpls(Cornell.bootYX,indices=2:8))
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints
#> Warning: extreme order statistics used as endpoints

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.bootYX$t[,-1])~factor(rep(1:7,rep(250,7))), 
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.bootYX$t0[-1],col="red",pch=19)



library(boot)
boot.ci(Cornell.bootYX, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), 
index=2, verbose=FALSE)
#> BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#> Based on 250 bootstrap replicates
#> 
#> CALL : 
#> boot.ci(boot.out = Cornell.bootYX, conf = c(0.9, 0.95), type = c("norm", 
#>     "basic", "perc", "bca"), index = 2, verbose = FALSE)
#> 
#> Intervals : 
#> Level      Normal              Basic         
#> 90%   (-0.2044, -0.0498 )   (-0.1967, -0.0406 )   
#> 95%   (-0.2192, -0.0350 )   (-0.2119, -0.0262 )  
#> 
#> Level     Percentile            BCa          
#> 90%   (-0.2376, -0.0815 )   (-0.2216, -0.0718 )   
#> 95%   (-0.2520, -0.0662 )   (-0.2460, -0.0483 )  
#> Calculations and Intervals on Original Scale
#> Some basic intervals may be unstable
#> Some percentile intervals may be unstable
#> Some BCa intervals may be unstable
plot(Cornell.bootYX,index=2)

jack.after.boot(Cornell.bootYX, index=2, useJ=TRUE, nt=3)

plot(Cornell.bootYX,index=2,jack=TRUE)


rm(Cornell.bootYX)

# PLS permutation bootstrap

set.seed(500)
Cornell.bootYX <- bootpls(modpls, sim="permutation", R=1000, verbose=FALSE)
boot.array(Cornell.bootYX, indices=TRUE)
#>         [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#>    [1,]    7   11    5    2    1   12    4    9   10     6     3     8
#>    [2,]   10    6   12    9    8    2   11    4    7     1     5     3
#>    [3,]    6   12    1    5    2   10   11    3    4     9     8     7
#>    [4,]    5    7    9    4   11    2    6    8   12     1    10     3
#>    [5,]    5    9    7    6   12   10    4   11    8     2     1     3
#>    [6,]    8   12    7   10    1    4    3    2    5     6     9    11
#>    [7,]    1    2    5    6    4    7   12   10    8    11     9     3
#>    [8,]    9    1   12    7    2   11    8    3    5     6    10     4
#>    [9,]    2   10    5    6   11   12    8    1    4     3     7     9
#>   [10,]   10    1    9    7   12    3    4   11    5     8     2     6
#>   [11,]    8    3    6   12   10    7    1    9   11     4     2     5
#>   [12,]    9    2    4    1    3   10    6   12   11     7     5     8
#>   [13,]    5   10    1    3    2   11   12    4    7     8     6     9
#>   [14,]    7    4    3    2    9   12   10   11    8     5     1     6
#>   [15,]    9    5   11    3    4    6    8    2    1     7    10    12
#>   [16,]    3    9    1   12    4    5    6    7   11     2     8    10
#>   [17,]    1    3    8    2    5    6   11   10    4     7    12     9
#>   [18,]    5   12    2   10    4   11    9    3    1     6     7     8
#>   [19,]   12    1    5    4    3    2    6   10    9     7    11     8
#>   [20,]   12   11    6    8    2    5    1    3    7     9     4    10
#>   [21,]    3    4   12    9    7    8   11    6    1     2     5    10
#>   [22,]    4    2    1    5    6   10    7   11    3     8     9    12
#>   [23,]    2    9    5    1   12    3    8   10   11     6     4     7
#>   [24,]    9    5    7    4   12    2    6    1   10    11     3     8
#>   [25,]    5    8    4    1    6   12   10    9   11     2     3     7
#>   [26,]    5    7    8   12    3    2    4   10    6     9    11     1
#>   [27,]    5    8    1   11    6    9    3    2   12     4    10     7
#>   [28,]    3    8   12    4   11   10    2    5    7     6     1     9
#>   [29,]    9    5    7    1    2    4    6   12   10    11     3     8
#>   [30,]    7    3   11    6    2    1    8   12    9     5    10     4
#>   [31,]    9    1    6    2    8    3   12   11    4     5     7    10
#>   [32,]    3    8    7    9    4    6    2   10    1     5    12    11
#>   [33,]    3    4    6    8   12    1    7   10   11     2     9     5
#>   [34,]    9   10   12    8    3    4    6    5    2     1    11     7
#>   [35,]    7    8    6    3   12    4    9    5    2    10    11     1
#>   [36,]    5    9    6   11    7    4   12    8    3     1     2    10
#>   [37,]    9    6    8    3   12   10    4    5   11     1     7     2
#>   [38,]    7    9   12   10    5    6    2    3    1     8    11     4
#>   [39,]    6    1    2   10   11    4    7   12    9     8     5     3
#>   [40,]    6    1    5    9    3    7   10   11    4     2     8    12
#>   [41,]   11    1    4    8    5    9   10   12    7     6     3     2
#>   [42,]    2    3    1   12    8    4    9    6   11     5     7    10
#>   [43,]    4    3    6    5    8   11    2   12   10     7     9     1
#>   [44,]    4    7   11    3   10    5    2    1    6     8     9    12
#>   [45,]    2    1   10    9    3    6    5   12    4     8     7    11
#>   [46,]   12    8   10    7    3   11    4    1    9     6     5     2
#>   [47,]   11    3    7    5    8    4    9    2    1     6    12    10
#>   [48,]   11    5    1    3    2    8    7    4   12     6    10     9
#>   [49,]    2    3    4    8    7    1   11    5   10     9     6    12
#>   [50,]   12   10    5    2    6    3    1    7   11     8     9     4
#>   [51,]    1   10    4    6    5    7    3    8    9     2    12    11
#>   [52,]    6    9    7    2   11    8   12    3    4     1    10     5
#>   [53,]   11    9    8    4    6   12    5    1    7     2    10     3
#>   [54,]   12    9    8    3    4    1   11    6    7     2    10     5
#>   [55,]    4    7    6    3   10   12    2    9    1    11     8     5
#>   [56,]   10    8    4    5    9    7   12    6   11     2     1     3
#>   [57,]    8    5   12    7    9    2    6    1   11     3    10     4
#>   [58,]    9    1    6    5   12   11    7    3   10     4     8     2
#>   [59,]    3    7    4    9    6    2    1   12    5    10     8    11
#>   [60,]   12    7   11    4    8    9   10    6    5     2     1     3
#>   [61,]    9    1    3    7    4    5    8    6   10    11    12     2
#>   [62,]    8    5    2   11   10    1    7    4    6     3     9    12
#>   [63,]    1    7   10    3    4    6   12    8    9     2     5    11
#>   [64,]    6   11    8    7    3    4    5    1    2     9    10    12
#>   [65,]   10    5    4    3    7    6    8   12    1    11     9     2
#>   [66,]   11    7    6   12    3    2    9   10    4     8     1     5
#>   [67,]    5    7    9    2    6    3   12   11    1     8    10     4
#>   [68,]   11    3    8    6   10    4    7    2    1     9    12     5
#>   [69,]    5    7    4   12    2    3   11   10    9     1     8     6
#>   [70,]    1    4    3   12   10    6    5   11    2     8     9     7
#>   [71,]   11    2   12    7    8    1    6   10    4     3     9     5
#>   [72,]    4    5    9    8    7   10   12    2    3    11     6     1
#>   [73,]    8   10   11    5    7    4    9    1    2     6    12     3
#>   [74,]    9    3    5    2    8    4    1    6   12    10    11     7
#>   [75,]   12    2    8    9    7    4    3   10    6     1     5    11
#>   [76,]    4    6    3    1   10    2    8   12    5     9     7    11
#>   [77,]    1    5   11    4    2    6   10    3   12     8     9     7
#>   [78,]    2    3    8   10    7   12    1    6    5     9    11     4
#>   [79,]   11    4    8   10    6   12    5    7    3     1     9     2
#>   [80,]    4    8    3    6    9    1   11    7    2    12     5    10
#>   [81,]    6    9    7   10    3   11    1    8    5     4    12     2
#>   [82,]    4    9   11    1   12    5    7    8    3    10     2     6
#>   [83,]    9    3    1   11    4    6    8    5   12    10     7     2
#>   [84,]    1    2   11    5    9    4    6    8   10     7     3    12
#>   [85,]   12    7    9    5   10    2    4    1   11     8     6     3
#>   [86,]    7   12    3    5    2   11    1    4    6     8    10     9
#>   [87,]    8    1   11   12    2    4    9    6    5    10     7     3
#>   [88,]   10    2    5   12    6    1    3   11    4     9     7     8
#>   [89,]   12    1    4    2    8    5    9    7   11     6    10     3
#>   [90,]    1   12   11    9    6    5    4    2    8    10     7     3
#>   [91,]   11    2    1    7   12    8    3    4   10     5     6     9
#>   [92,]    5    1    8    2    4    6   12   10    9     3     7    11
#>   [93,]    2    3    1    5   10   12    4    9    7     8     6    11
#>   [94,]    6    1    2    5   12   10   11    8    4     7     9     3
#>   [95,]    7    6   12    5   10    4    3    2    1    11     8     9
#>   [96,]    2    3    9   10    7    1    5    6    4    12     8    11
#>   [97,]    8   12    3    4   10    9    2    1   11     7     5     6
#>   [98,]    8    6    2    1    7   10   11    5    9     3    12     4
#>   [99,]   10    8    1    9   11    5    3    7   12     2     6     4
#>  [100,]    4    7    1   10    6   11    3   12    2     5     8     9
#>  [101,]    6    7    9    2   10    8    3   12    1     4     5    11
#>  [102,]    9   12    6    7    3    4   11    1    5     2     8    10
#>  [103,]    2    3   12    1   11    6    5    7   10     8     9     4
#>  [104,]    6    8    7    1   11    4    2   12    3     5    10     9
#>  [105,]   11    9    4    3    8    6    7   10   12     2     1     5
#>  [106,]    5    6   11    2   12    8   10    9    7     4     3     1
#>  [107,]    3    8    4    1   12    2    9    6   11     7     5    10
#>  [108,]    2    8    1   10   12    4    9   11    5     7     6     3
#>  [109,]    7    4    6    2   11    9    1   10    3     8    12     5
#>  [110,]    2   11    1   12    7    4    8    3   10     9     6     5
#>  [111,]    5    7   12    6    2    1    9    3    4     8    10    11
#>  [112,]   10   12    4    9    3    8    1    6   11     5     2     7
#>  [113,]    5    9    7    1   10   11    8   12    3     6     4     2
#>  [114,]    4    2    3    6   10   11    7    8    1    12     5     9
#>  [115,]   12    1    9   11    3    8    6    4    5     7     2    10
#>  [116,]    9    6    4    3   11    5    2    1   12     8    10     7
#>  [117,]    8    9    5    7    1   10   11    3    6     2    12     4
#>  [118,]    3   12    9    2   10    1    5    8    7     4    11     6
#>  [119,]   10    1   12    3    5    8    6    4    9     2    11     7
#>  [120,]    1   11   10    2   12    3    7    5    4     8     9     6
#>  [121,]    8    5    3    7    2    9    4   12   10    11     6     1
#>  [122,]    7    9   10    5   12    6    3   11    8     4     1     2
#>  [123,]    2    8    5    1    3   11   10    7    9     4    12     6
#>  [124,]    8    9    2   11    5    7    6    4    3     1    10    12
#>  [125,]    8   10    1    2    5    7    6    4   12     9    11     3
#>  [126,]    3    5    6    7    2    1   10   12    9    11     4     8
#>  [127,]    4   10    3    1    7   12    5    2    8     9    11     6
#>  [128,]    7    3   10    2    5    8    4   12    9     1    11     6
#>  [129,]    7   10   11    3    1    4    9    2    6    12     8     5
#>  [130,]   12    6    1    2    3    7   11    5    8    10     9     4
#>  [131,]    8    1   11    9    5    7    2   10    3    12     6     4
#>  [132,]    9    5   11    8    3    6   10    2    1    12     7     4
#>  [133,]    8    5    9    7   11    1   10   12    6     2     3     4
#>  [134,]   11    6    9   10    3    8    5    7    2     1     4    12
#>  [135,]    8   11   12    9    3    2   10    1    7     5     4     6
#>  [136,]    9    2   10    4   12    8    6    1    5     3     7    11
#>  [137,]   10    7    9    8    4    3    1   11    6    12     2     5
#>  [138,]   12    2    7   10    5   11    8    6    1     4     3     9
#>  [139,]   11   10    8    3    9    4    2    5    6     7    12     1
#>  [140,]   11    9    4    5    3    6    1    8   12     2     7    10
#>  [141,]    2    6   12    1    3    7    4    5    8     9    10    11
#>  [142,]   12    1    3    8   10   11    5    4    9     7     2     6
#>  [143,]    6    8    1   12    3    7    2    9    5    10     4    11
#>  [144,]   11    1    8   12    5    2    3    9    6    10     7     4
#>  [145,]   10    6    7    9    1    5   11   12    4     8     2     3
#>  [146,]    6    5    8    9    1    4    3   11    7    10     2    12
#>  [147,]    6    2    1    7    8    9   10    5   11     4     3    12
#>  [148,]   10   12    9    3    5    2    8    6    4    11     7     1
#>  [149,]   10    3    4    7    9    8   11    2    6     5     1    12
#>  [150,]    1    9    2    4   10   12   11    8    6     7     5     3
#>  [151,]    4    1    6   11   10    3    8    5    9    12     2     7
#>  [152,]    8    3    2    7    6    4   10    5    1    12    11     9
#>  [153,]    1    2   11   12    6    8    3    4    5    10     9     7
#>  [154,]    1    9    6    5    3    4   10   11    8     2     7    12
#>  [155,]    9    7    5    3    2   12    8    4   11     1    10     6
#>  [156,]   11    8    3    2   10    6    5   12    7     4     1     9
#>  [157,]    4    2   10    5   11    3    7   12    8     6     9     1
#>  [158,]    2   12    4    6   11    1   10    9    8     5     3     7
#>  [159,]    7    8    5    9    4    3    6    1   11    12    10     2
#>  [160,]    2    1    8   11    5    9    7    4   12    10     6     3
#>  [161,]    7    1   12    6    8    9   10    5   11     4     2     3
#>  [162,]   10   11    6    7    2    5    1    3   12     4     9     8
#>  [163,]   12    8    6    1    7    2   11    3    5     4    10     9
#>  [164,]    6    2   12    8    1    5   10    3    7     4     9    11
#>  [165,]    1    3   11   12   10    9    5    4    2     6     8     7
#>  [166,]   12    5    1    2    8    4    9    3    7    10     6    11
#>  [167,]    1    8    5   11    3   12    7    9    2     4    10     6
#>  [168,]   10   11    3    4   12    5    2    7    9     6     1     8
#>  [169,]   11    9    5    1    4   10   12    3    6     2     7     8
#>  [170,]    8   10   12    4    3    2    1    7   11     9     5     6
#>  [171,]    2    5    3    8   10    7    9    4   12    11     6     1
#>  [172,]    8    6    5   11    9   10    1    2    7    12     3     4
#>  [173,]    4   10   12    8    2    9    7    1    5    11     6     3
#>  [174,]   11   10    6    5    3    2    4    7   12     9     1     8
#>  [175,]   12    3   10    5    1   11    8    9    4     7     6     2
#>  [176,]    9    6    7    8   11    3    5   10    2    12     4     1
#>  [177,]    5    8    2    3    6   11   10    9    1     7     4    12
#>  [178,]    5   12    8    3    9    4    1   10   11     6     2     7
#>  [179,]    3   12    6   10    9    8    1   11    2     4     7     5
#>  [180,]    4   12    1    3    8   10   11    5    7     2     9     6
#>  [181,]   11    4    5    6    3    2    1    8   12     9    10     7
#>  [182,]    7   12    3    1    6    2    5   11    8    10     9     4
#>  [183,]    1    4    9    3    5   10   11    7   12     8     6     2
#>  [184,]   12    9    3    5   11    6    8    4    1     7    10     2
#>  [185,]    7    4    6    8    9    5    1   10   11    12     2     3
#>  [186,]    2    9    8   12    1   11    5    3   10     4     7     6
#>  [187,]    2    4   11    5   10    8    7    1    3     9     6    12
#>  [188,]    4    2    9    7    1    8   11    3    5    12     6    10
#>  [189,]    6   11    1    8   10    7    5   12    3     2     9     4
#>  [190,]   12    7   10   11    8    3    4    6    5     9     1     2
#>  [191,]    5    4   11    2   10    1   12    3    9     7     6     8
#>  [192,]    1    8    2    7   10    5    6   11    4     3     9    12
#>  [193,]   10    2    7    5    6   12   11    3    4     9     1     8
#>  [194,]    4    8    5   10   11    6    3   12    7     1     9     2
#>  [195,]    3   12    9   11    2    8   10    6    5     7     1     4
#>  [196,]    8    4   12    2    3    9   10    7    5     6     1    11
#>  [197,]   10    3    5    2    9    4    7   12    6     1     8    11
#>  [198,]    1    6    2   10    9    8    7    4    5    11     3    12
#>  [199,]    3    1   10   11    7    6    2   12    8     5     4     9
#>  [200,]   11    4    3    5    9   10   12    1    6     7     8     2
#>  [201,]    4    8    6   11    2    3    1    9    7    12    10     5
#>  [202,]    7   11    2    5   12    3   10    6    1     4     9     8
#>  [203,]    4    8    6    5   11    9    3    1    2    12    10     7
#>  [204,]    5   10    7    3    9    2    4    8   11     1    12     6
#>  [205,]    8    9   10    5    4   12   11    2    3     6     1     7
#>  [206,]    3   10    9    5    6    8    4    1   11     2    12     7
#>  [207,]    2    7    5    3    9   12    6    1   11    10     8     4
#>  [208,]    1    4   12    2    9    7    8    6    3     5    10    11
#>  [209,]   10    9    6    7   12    1    8    2    3     5    11     4
#>  [210,]    1    9    8   12   11    2    7    6    4    10     3     5
#>  [211,]    1   10    5    7    4    8    2    3    6    12     9    11
#>  [212,]    3    1    4    7   12    2   11   10    6     5     8     9
#>  [213,]    3    9    2   10   12    5    6    8   11     1     4     7
#>  [214,]    4    2    3   12   10   11    1    9    8     6     7     5
#>  [215,]    2    5    3    4    1    6   10    8    9     7    12    11
#>  [216,]    1    2    4   10    9    6    3    5    7    12    11     8
#>  [217,]   11    7    6    1    2    4    9    3    8    10    12     5
#>  [218,]    5    2    3   11    6    1    8   12   10     7     4     9
#>  [219,]    8    2   11    3    6    4   10   12    5     9     7     1
#>  [220,]   12    1   10    2   11    4    7    8    6     9     5     3
#>  [221,]    3    4    9    2    5   11   12    1   10     8     7     6
#>  [222,]    9    5    6    4    2    7    8   11   10     3    12     1
#>  [223,]   11    2    5    9    1    8    3   10    6     7     4    12
#>  [224,]    4    7    9    5   10   12    2    8    1    11     3     6
#>  [225,]   10    3   12    1   11    2    8    7    9     5     4     6
#>  [226,]    7    3    9    1   10    6    5    2    8     4    11    12
#>  [227,]    8    3    2    5    4   11    1    9   12     6    10     7
#>  [228,]    8    3    5    9    4   11    2   12    1     6    10     7
#>  [229,]   10    4    5    2    8    1    3   12    6    11     9     7
#>  [230,]    2   11    6    1    4   12    3    8   10     7     9     5
#>  [231,]    2    9    4    8    1    3    6   12    5     7    10    11
#>  [232,]   12    2    3    7    1   10    8    4    6     5    11     9
#>  [233,]    5   11    6    4    7   12   10    1    8     9     2     3
#>  [234,]   10    1    6    4    3    8   11    7    5     2    12     9
#>  [235,]    6    3    7    8   12    2    5    4    1    10     9    11
#>  [236,]   11    5    4    6    7    1    9    2   12     8     3    10
#>  [237,]    6    5    4    2    7   12    8    1   11     9     3    10
#>  [238,]   10    4    6   12    7    1    3   11    5     8     9     2
#>  [239,]   12    2    7    3    8    9    6    5    4    11    10     1
#>  [240,]    2    3    9    5    4   10   12    6    1     8     7    11
#>  [241,]    7    1    4    8   11    6    5    9   12     3    10     2
#>  [242,]    1   11   10    4    5   12    2    3    8     7     6     9
#>  [243,]    3    6    7    2    1    9   12    5   11    10     8     4
#>  [244,]    8    2    4    9    5    3    6   12   11    10     7     1
#>  [245,]    9   10    1    5    8    4   12    6    3     2     7    11
#>  [246,]   11   12    3   10    6    9    1    2    5     4     8     7
#>  [247,]    9    1    6    8   11    7    2    5   12    10     4     3
#>  [248,]   12    9   11    2    7   10    1    4    3     5     8     6
#>  [249,]    5    6    8   12    1    3   10    9   11     7     4     2
#>  [250,]    7    5    4    1   11    6   10    2    9     8    12     3
#>  [251,]   11   12    6    8    7    3    9    1    2     4     5    10
#>  [252,]    9    1    4    2    8    5    3   10    7     6    11    12
#>  [253,]   11    1    3    7   10    8    9    6   12     4     2     5
#>  [254,]    3    8    5   12    6    1    2   11    9    10     4     7
#>  [255,]    3    8   10   11    9    4    2    6   12     7     1     5
#>  [256,]    4    1   12   10   11    5    7    2    9     3     6     8
#>  [257,]    6    9    1    2    7    8   11    3    4     5    12    10
#>  [258,]   10    4    7    2   12    5    1    8   11     3     6     9
#>  [259,]   12    1    9   10    8    3    7    4    6     2    11     5
#>  [260,]   10    9    1    7    5    4    8    3    6     2    12    11
#>  [261,]    2    3   12    9    5    1    4    7   10     6     8    11
#>  [262,]    6    5    9   10   12    2    7    3    8    11     4     1
#>  [263,]    8   12   10   11    1    2    5    6    4     3     7     9
#>  [264,]    9    1    6    7    3    8    2    4    5    10    12    11
#>  [265,]    8    4   10   12    5   11    3    6    2     7     1     9
#>  [266,]    4   11    9   10    3    2    5    8    7     6    12     1
#>  [267,]    3   12    8   10    9    2    6    1    4     5    11     7
#>  [268,]    4   10   11    7    1    2   12    9    3     5     6     8
#>  [269,]   10    3    1    6   11    2   12    9    8     4     7     5
#>  [270,]    4   12    6    1   11    2    8    7   10     3     5     9
#>  [271,]    2    6    1   12    9   10   11    3    8     4     7     5
#>  [272,]    6    4   10    2   12    1    7    3    8    11     5     9
#>  [273,]    4   11    5    1   12    6    9    8    2     3     7    10
#>  [274,]    8   12    6    5    3    7   11    2    1    10     9     4
#>  [275,]    1    8    7   12    9    4    6    2   10     5     3    11
#>  [276,]   12    5    7   11    2   10    4    8    3     9     1     6
#>  [277,]    4    5    1    9    7    2   11   10   12     3     6     8
#>  [278,]    1   10    3    2    6   12    8    9    4    11     5     7
#>  [279,]   10    1    9    5    3   11    6   12    2     8     7     4
#>  [280,]    6   10   12    4    2    9   11    5    1     3     8     7
#>  [281,]    9    4    3   12    5    1    8    6    2     7    10    11
#>  [282,]    2   11   10   12    9    8    7    1    5     6     4     3
#>  [283,]    8    3    6    1    4    5    2    9   11    10     7    12
#>  [284,]    1   11    4    2   12    5   10    8    7     3     9     6
#>  [285,]   12    1   10    6    2    8    4    3    7    11     9     5
#>  [286,]   10    1    6    2    7    3   12    9    5    11     8     4
#>  [287,]    5    1   11    4   12    8    9    2    6     7    10     3
#>  [288,]    2    4    5   10    9    3    1   11    8     6     7    12
#>  [289,]    2   10    3    6    5    7    1    9   11     8    12     4
#>  [290,]   11    6    3    8   10    7    4    9    2    12     5     1
#>  [291,]   10   11    3    7    1    4    2   12    5     6     8     9
#>  [292,]    8    7    4    1    3   11    5    9   10     6    12     2
#>  [293,]    6   10    3   11    9    2    8    7    1     4    12     5
#>  [294,]   10   12    3    2    9    4   11    7    8     1     6     5
#>  [295,]    2    3    1   12    6   10    7   11    9     8     5     4
#>  [296,]    5    2    6    7   10    4   11    9    8    12     1     3
#>  [297,]    4    6    3   11    9    8    7    1   12     2     5    10
#>  [298,]    5   10    8   12    2    7    1    3    6    11     4     9
#>  [299,]    9    7    2    6   10   11    1    3    5     8     4    12
#>  [300,]    2    9   12    5    4    6    8    3    7    11    10     1
#>  [301,]    4    5    9   10    2   12    1    8    6    11     3     7
#>  [302,]    5    7   10   12   11    3    4    9    6     8     2     1
#>  [303,]    8   12   10    4    5    9    6    1    7     3     2    11
#>  [304,]    7    3    2   11    9   12    5    6   10     4     8     1
#>  [305,]    1   12    6    8    4    5   11    9    3     7     2    10
#>  [306,]   12    3    8    4    9    1    5    6   10     2     7    11
#>  [307,]   10    1    7    6    3    4    9    2   12     8     5    11
#>  [308,]    6    7    2   11   10    4    3    5    1    12     9     8
#>  [309,]    1    8   12    7    9    5    6    3    4    10     2    11
#>  [310,]    1    2   10    7    8    4   12    3    5    11     9     6
#>  [311,]    6   12   11    5    8    3    1    4    7     2     9    10
#>  [312,]    5   11    7    4   12    2    8    9    1     6    10     3
#>  [313,]    6   10   11   12    9    7    4    5    3     1     2     8
#>  [314,]    3    6    8   10   11    9    2    5    4    12     1     7
#>  [315,]    8    6   11   10    3    7    4    9    2     5    12     1
#>  [316,]    2    3    1    5   10    8    6    7    9    11     4    12
#>  [317,]   11    7    6   12    2    5    1    8    9     3    10     4
#>  [318,]    5   12    8   10    6    4   11    7    3     1     9     2
#>  [319,]    6    5    9    4    2   10    8    7    3     1    12    11
#>  [320,]    2   11    6    4    7    5   12    8    1    10     3     9
#>  [321,]    2    6   10    8    3    5    9   11    1     4    12     7
#>  [322,]    9    1   11    2    6    4    8    3   12    10     7     5
#>  [323,]    2    6    3    1    5    8   11   12    4    10     7     9
#>  [324,]    3    2    6    9    5   11   10    8    4     7    12     1
#>  [325,]    3    6    1   11    7    8    5    4   12     9    10     2
#>  [326,]    1    2    6   12    4    5    8   10   11     9     3     7
#>  [327,]   12   10    4    3    2    7    5    6   11     9     8     1
#>  [328,]    9    2    3    5    7    6    4   11   10     1     8    12
#>  [329,]    2    8    1    7    3    9    5   12   10     4     6    11
#>  [330,]    3    6    1   12    2    4    7    8    5    11     9    10
#>  [331,]    9    2   11   10    6    1    4   12    7     3     8     5
#>  [332,]   11    2    8   12    7    5    1   10    6     9     3     4
#>  [333,]   10    2    6    3   12    8    5    4    1    11     7     9
#>  [334,]    7   10    5   11    4   12    1    8    6     3     2     9
#>  [335,]    6   12    3    1    8    4   11   10    2     5     9     7
#>  [336,]   11    4    7    8   12    6   10    9    5     3     1     2
#>  [337,]   10   12    6    5   11    7    4    8    2     3     9     1
#>  [338,]    4    5   10    2    8   11   12    7    9     1     6     3
#>  [339,]    4    1    5    8    7    9    6   12    2    10     3    11
#>  [340,]    4   10    5   11    3    1    2    7    6     9     8    12
#>  [341,]    5   10   12    1   11    6    9    4    2     7     3     8
#>  [342,]    6    8    3    5    1    4   12   10    9     2    11     7
#>  [343,]    7    8    9    4    6    2   11    1   12     5    10     3
#>  [344,]    7    9   10    6    8    2   11    4    5    12     1     3
#>  [345,]    6    7    3    9    2   12   10   11    1     5     4     8
#>  [346,]    6    1    4    5    3    2    8   10    7    12     9    11
#>  [347,]    4    6    2    7    8    9   10   12    1    11     3     5
#>  [348,]    3    8    6    4    5    9    7   11    2    10    12     1
#>  [349,]    5    2    4   12    9    7    3    8    1    11    10     6
#>  [350,]    7    9   12   11    3    2    6    8    5     4     1    10
#>  [351,]    9    4   10    1    6    8    5    3    7     2    12    11
#>  [352,]    7   11    4    8    3    5    1    9    6    10    12     2
#>  [353,]   12    5    3   10    1   11    9    7    2     4     8     6
#>  [354,]    3    1    5    9    8    6    4   12   11     2    10     7
#>  [355,]    8    4    5    3   11   12    7    2    1    10     9     6
#>  [356,]    6    4   10   11    7    5    1    2    9     8     3    12
#>  [357,]    2   10   11    6    8    5    1    9   12     7     4     3
#>  [358,]    5   12    2   10    7    4    8    6    9     1     3    11
#>  [359,]    7    9    6    3   10    4    1    8    5    12     2    11
#>  [360,]   11    9   12   10    4    3    5    6    7     8     1     2
#>  [361,]   11    9    8   12    5    4    1    2   10     3     7     6
#>  [362,]    7    3   11    9    4   10    2    5    8     6     1    12
#>  [363,]   11    6   10    1    4    2    9    3    5     7     8    12
#>  [364,]    9    7    2    4    5   12   10    8   11     3     6     1
#>  [365,]    2    9   11    5    8    3    6    4    1    12     7    10
#>  [366,]    5    7    9    4   12    3    1   11   10     8     6     2
#>  [367,]    7    6    1   11    4    8    9    5   12     2    10     3
#>  [368,]    7    2    9    6    8   11    3   10    5    12     1     4
#>  [369,]    7    1    8    2   12    4    6   11    5    10     3     9
#>  [370,]    9    7    6    8   12   10    2    5    3    11     1     4
#>  [371,]    4    6    2    7   10   11   12    3    9     5     1     8
#>  [372,]   11    4    3    1    5    7    2   12   10     8     6     9
#>  [373,]    2    7    1   11   10   12    8    6    5     4     3     9
#>  [374,]    9    2    1   12    4    3    7   11    5     6    10     8
#>  [375,]    7    8    5    9    2   10   12    3   11     1     4     6
#>  [376,]    8    3    4    1    6    2   10    7    5    11     9    12
#>  [377,]   10    3    4    6   11    2    8    7    1     5    12     9
#>  [378,]   10    3    8    2    5    4    9   12    6     7    11     1
#>  [379,]    8    4   12    2    5    3    6    1    9     7    10    11
#>  [380,]    4    6   11    7    8    2   10    5    1     3     9    12
#>  [381,]    6   11    2    5    3   12    4    8   10     7     1     9
#>  [382,]   12   11    5    3    8    6    2    9    4     1    10     7
#>  [383,]    5    8    9   11    3   12    6    1    7     2    10     4
#>  [384,]    3    4   11    7    9   12    6    1    8    10     5     2
#>  [385,]   10   11    3    7    2    4    6    5    8     1     9    12
#>  [386,]    3    2   11    1    6    5   12    7    9     4    10     8
#>  [387,]   12   11    4    9    2    6    5    3    8    10     1     7
#>  [388,]    6    7    3    9    2    4    1    8   12     5    11    10
#>  [389,]    9    8    2    3    5   11    4    7   10     6     1    12
#>  [390,]   10    7    4    9    2    1    6    8    3    11    12     5
#>  [391,]    6    4   10    9   11    5    3    2    1     7    12     8
#>  [392,]    6   12    9    1   10    7    4    5    3    11     8     2
#>  [393,]    3    8    9    5    6   10    2    1   12     7     4    11
#>  [394,]    7    2   12    5   11   10    9    1    4     8     3     6
#>  [395,]   11    4    8   12    9   10    1    3    6     5     2     7
#>  [396,]    6    1    4   10   12    2    5   11    7     8     3     9
#>  [397,]    8    4    2    5    9    3   10    6    1    11    12     7
#>  [398,]   12    1   11   10    6    4    2    8    5     3     9     7
#>  [399,]    1    3    7    4    6   10    9    2   11     5     8    12
#>  [400,]    4   12    7    6    2   11    5    9    1    10     8     3
#>  [401,]   12    6    4    5    3    1    2    7   11    10     8     9
#>  [402,]    8    4    3   12    5    1   11    7   10     2     6     9
#>  [403,]    6    2    5   10   11    8    3    7    1     9     4    12
#>  [404,]    4    3    1   10    6    7    9    5    8    11     2    12
#>  [405,]   12    6    8    3    9    7   11    1   10     4     2     5
#>  [406,]    6   12    8    2   11    7    5    9    1    10     3     4
#>  [407,]   11    8    2    9   10    1    6   12    7     3     5     4
#>  [408,]    3   10   11    1    4    6    7    2    8     5    12     9
#>  [409,]    1   12    3    8   11   10    2    9    4     6     7     5
#>  [410,]   11    9    1    7    2    5    8    3    6    10     4    12
#>  [411,]   11    5    2    3    8   10    7    6   12     1     4     9
#>  [412,]   12    5    3    9    8   10    1    4    2    11     6     7
#>  [413,]    6   12    1    5   10    7    3    8    4    11     2     9
#>  [414,]    3   11    7    8    2    5    9    1   10     6    12     4
#>  [415,]    6   11    1    2    7   10    5    8    9     3    12     4
#>  [416,]   11    9    3    2   10    1    8    7    4    12     6     5
#>  [417,]    6    4   12    1    9    7    2    5    3    10    11     8
#>  [418,]    5    8    2    7   11   12    3   10    6     9     1     4
#>  [419,]    9    5    8    7   12   11    6    1    2     4     3    10
#>  [420,]   12    9   10    3    1    6    4    7   11     8     5     2
#>  [421,]   11    8   12    6    3    1    7    5    4     9     2    10
#>  [422,]    3   10    8   12    4    2    6    5    9     1     7    11
#>  [423,]    2    6   12    3    5   11   10    9    7     8     1     4
#>  [424,]    7    1   10   11    3    5    2    9    6     4    12     8
#>  [425,]    5   12    9   11    4    3   10    2    8     7     6     1
#>  [426,]    4    7   12    3    8    1    5    6   10     9    11     2
#>  [427,]    1   12   11    5    4    7   10    8    6     2     9     3
#>  [428,]    3    2    5   11    4   12    8    9    7     1    10     6
#>  [429,]    2   11    5    4    1    3    7    8   12    10     6     9
#>  [430,]    8    6    5    4    2    7   11   10   12     3     1     9
#>  [431,]   10    6    5    4   11    9    3    7    2     1    12     8
#>  [432,]    9    4    7   10    1    5    2    8    3    11     6    12
#>  [433,]    3    1    2    6   12    8   10    7   11     5     4     9
#>  [434,]    3   12   10    7    8    6    2    9   11     5     1     4
#>  [435,]    1    8    2   12    7   11    5   10    6     9     4     3
#>  [436,]   12    1    4    3    5   11    2    7    8    10     6     9
#>  [437,]    4    1    7    9   12    2    3    8    5    11     6    10
#>  [438,]    3    6   12    2    9    1    5    7    4     8    11    10
#>  [439,]   12   10    9    6    3    1    4    7    5     2     8    11
#>  [440,]    9    1    5    2    4    7    8   10   12     3    11     6
#>  [441,]    9    8    3   10    4    5    2    7    6    11    12     1
#>  [442,]    9    4    3    1    8   12   10    6    2     5     7    11
#>  [443,]    2   11    3    4    8    9    6   10    5     1     7    12
#>  [444,]    9    7   11   12    5    3    1    8    6     4    10     2
#>  [445,]    3    4   12   11    9    6    7   10    5     1     8     2
#>  [446,]    5   10    8    7    1    3   12    4    9     2    11     6
#>  [447,]    8    6    9   10    3    7    1    2    4    12     5    11
#>  [448,]    7    1    3   10    8    9    4    6    5    12    11     2
#>  [449,]    7    5    6    2   10   12   11    4    9     1     8     3
#>  [450,]    5    8    6    7    4    9    2   12   11    10     3     1
#>  [451,]    8    5    3    1    6    9    7   10    4    12    11     2
#>  [452,]    6    4   10   12    1    9    8   11    5     7     3     2
#>  [453,]    2   10    9    3   11    4    7    5   12     1     8     6
#>  [454,]    5    6    7    4   10    8    2    3   12    11     1     9
#>  [455,]   12    9   10    8    5    6    7    4    1     3    11     2
#>  [456,]    1    8   10    9   12    5    4   11    2     3     7     6
#>  [457,]    5   10   11    1   12    9    2    4    8     7     3     6
#>  [458,]    6    1    8   12    9    2   11    7    3     4     5    10
#>  [459,]    7   12    9    4    3    5    1    6   11    10     2     8
#>  [460,]    1    9    7   12    5    3    8    2    6    10    11     4
#>  [461,]    2    3   10    9    5    8   12    1    6     7    11     4
#>  [462,]   11    7    9   10    5    4    1   12    3     8     2     6
#>  [463,]    5    1   11   10   12    9    6    2    4     8     7     3
#>  [464,]   10    1    5    3    6    9   12    2   11     7     4     8
#>  [465,]   10    2    8    4    9    3    1    5    6    11    12     7
#>  [466,]    7    5    3    6    2   10   12    1    8    11     4     9
#>  [467,]    9   11    4    1    5   12    2    7    3    10     8     6
#>  [468,]    6    9    4    3   12    7    1    8    5     2    10    11
#>  [469,]    9    1   10    6    8    5   12    7    3     2    11     4
#>  [470,]    3   11   12    6    1    4   10    9    5     2     7     8
#>  [471,]    4   10    5    6   11    1    8    3    9     7     2    12
#>  [472,]   12    5    4    6   10    8    9    7    1     2     3    11
#>  [473,]    1    2    5    7    9   10    6    3   12     8    11     4
#>  [474,]    2   12    3   11    1    5    8    4    6     7     9    10
#>  [475,]    4    5    6    1    9    2   10    7   12    11     8     3
#>  [476,]    7   12    2    8    9    6    1    5    4    10    11     3
#>  [477,]    9    5    1   12    7   10   11    4    2     8     6     3
#>  [478,]    8    7    6   12    4    5    2    3   10    11     9     1
#>  [479,]    4    6    3    1   12   10   11    8    7     2     9     5
#>  [480,]    8    1   11   12    7    3    6    5   10     4     2     9
#>  [481,]   10   12    8    4    5    3   11    2    7     6     1     9
#>  [482,]    6   10    9    2   12    7   11    4    3     1     8     5
#>  [483,]    7    5   11    3    8    1   12    4    6    10     9     2
#>  [484,]    4    5    3    2   10    7    9    6   11    12     1     8
#>  [485,]   11   10    1    6    2    5    7    3    8     4    12     9
#>  [486,]    5    6    3   11   10    8    4    7    2     9    12     1
#>  [487,]   10   12    3    6    1    5    7    9   11     4     8     2
#>  [488,]    1    4   10    2    6    5    8   12    7     3     9    11
#>  [489,]    2    6   11    3   12    7    5    9    8     4    10     1
#>  [490,]    8    6    1    7   12    9    2    3   11     4    10     5
#>  [491,]    5    9    6    8   11    7    1   10    3    12     2     4
#>  [492,]    2    8    5    4    3    6    9    1   11    10     7    12
#>  [493,]    5   11    2    7    9    3    6   12   10     4     8     1
#>  [494,]    4    1    7   11    8    3    5    9   10     6    12     2
#>  [495,]    1   10    9    3    7    8    5   11    6     4    12     2
#>  [496,]    8    7    2    6    3    9    4   10    1    11     5    12
#>  [497,]    7   10    9    6    1    8   11    3   12     4     2     5
#>  [498,]    6   11    4    1    7    3   10    9    2    12     8     5
#>  [499,]   12    9   10    5    4    2    3   11    7     6     1     8
#>  [500,]    4    1    8   11    5    3    9    7   12     2    10     6
#>  [501,]    2   11    8   12    3    5    6    4    1    10     7     9
#>  [502,]    9   12    5    2    6    1   10    4    8     7     3    11
#>  [503,]    6    1    3    4   10    8   11    9   12     2     5     7
#>  [504,]    1   11    2   12    3    5    4    8   10     9     7     6
#>  [505,]    3    5    6    2    1   12   10    7   11     9     8     4
#>  [506,]    8    6    9    4   10   11    2    5    3    12     7     1
#>  [507,]    8    9    7    4   11    5   12    2    1    10     6     3
#>  [508,]   10    1    3    2    9   12    5    7   11     4     8     6
#>  [509,]    2    6    3   11    8    5   10    7   12     4     9     1
#>  [510,]   11    4   12    2    5    8    9    3    6     1    10     7
#>  [511,]    2    8    6    7   12    4   10    5    3     9     1    11
#>  [512,]    7   11    4    6    5   10    3    2    8     9     1    12
#>  [513,]    8    9    1    5    4    6    7   10    3     2    12    11
#>  [514,]    7   11    5    1    8    6   10   12    2     4     9     3
#>  [515,]    3    2    7    1    4    8   11    5    9    12     6    10
#>  [516,]   10    9    2    7    5    1   12    3    4    11     6     8
#>  [517,]    4    8    9   10   12    1    6   11    3     2     5     7
#>  [518,]    2    9   12    3   10    4    6   11    5     8     1     7
#>  [519,]    9    5   10   11    1    7    3    2   12     8     6     4
#>  [520,]    8    6    7   10   11    2    5    4    3     1    12     9
#>  [521,]    1    9   11   12    7    4    5    2   10     3     6     8
#>  [522,]    9    7    4   11    8    3    5    2   12     6    10     1
#>  [523,]   11    4    2    5   10    9    6   12    1     7     3     8
#>  [524,]    4   11   10    6    7    2    8    9    5     3     1    12
#>  [525,]   10    9   12    4    7    6    1    2    8     3    11     5
#>  [526,]    2   11    7    8    3    5    4    1   10     9    12     6
#>  [527,]   11    9    5    2    7    8   12    4    3     1     6    10
#>  [528,]    9    4    3    7    1   10   12    8    2    11     6     5
#>  [529,]    7    3    1   11    2    8   12    4    9     5    10     6
#>  [530,]    8    5    9   12    3   11    4    7    1     2     6    10
#>  [531,]    7    2   11    6    3    4    5   10    9     1    12     8
#>  [532,]    1   10    7    2    6   12   11    8    3     5     9     4
#>  [533,]    9   12    3   11    7    5    8    2    6     1     4    10
#>  [534,]    3    5    7    1    9    4    8   12   10     6     2    11
#>  [535,]   10    5    8    1    2    7    6    4    9     3    12    11
#>  [536,]    8    5   10    7   11    4    1   12    9     6     2     3
#>  [537,]    5    9    8    7   12    4   11   10    3     2     1     6
#>  [538,]    6    7   12    5   11    3    9    4    2    10     1     8
#>  [539,]    7    6   10    5   11    1    4    3    2     9    12     8
#>  [540,]   11   12    8    2    1    5    9    3    4     6    10     7
#>  [541,]   12    6    7   11    4    5    1    2   10     9     8     3
#>  [542,]    2   10   11    9    6    1    5    4    8    12     7     3
#>  [543,]    2    3    6   11    4   12    5    8    7     1    10     9
#>  [544,]   10    4    6    8    2    9    3    5    1    12     7    11
#>  [545,]    2    5    9    7    4    3   12    1    8    11     6    10
#>  [546,]    6    8    1    7   10    2   12   11    5     3     4     9
#>  [547,]    5   11    1    8    3    6    9    7   10     4    12     2
#>  [548,]   11    4    2    8   10    9    7   12    5     1     6     3
#>  [549,]   12    4    2    8   10    6    1    3    9     5     7    11
#>  [550,]    3   10    6    9    1   12   11    4    5     2     7     8
#>  [551,]    5    4   10    2    1   12    3    7    8     6     9    11
#>  [552,]    9    1    6   10    7    3    8    4    5    12     2    11
#>  [553,]    5    4   10    1    6    7   11    8    9    12     3     2
#>  [554,]    2    6    5    8    7   11    4    1    9     3    12    10
#>  [555,]    3   12    5    7    1    8    9   10    2    11     6     4
#>  [556,]    6    3   11    9    5    1   12   10    4     2     7     8
#>  [557,]    3   11   10    2    6   12    4    5    7     8     9     1
#>  [558,]    8   11    1    2    3   10    6    7    9    12     4     5
#>  [559,]   11    5    7    9    6   12    1   10    8     4     3     2
#>  [560,]    2   11    3   12    5    4    7    6    8     9     1    10
#>  [561,]    2    4   11   10    1    5    6   12    8     3     7     9
#>  [562,]   12    6   10    9    8   11    2    1    5     7     4     3
#>  [563,]   10    4    3   11    5    1    7    2    6     8    12     9
#>  [564,]    7   12    6    2    9    1    4    8    3    11    10     5
#>  [565,]    6   12    9    3    8    5    4   10    1    11     7     2
#>  [566,]    3    9    8    1   10    4    5    7   11     6    12     2
#>  [567,]    8    1    2    9    7    3    6    5   11    12     4    10
#>  [568,]    6   11    4    9    7    2   12    8    1     3    10     5
#>  [569,]    2    9   12    4   11    5    7    3    1     8     6    10
#>  [570,]    4    6   10    8    9    5   12    7    3     2     1    11
#>  [571,]    9   11   10    5    3    2    1    8   12     4     6     7
#>  [572,]    7    9   10   11   12    6    5    1    4     3     2     8
#>  [573,]   10    8   12    9    6    5   11    2    4     3     1     7
#>  [574,]   11    8    1    3    6    9    5    2   12    10     4     7
#>  [575,]   12    9    3    6   11    4    8    1    5     7     2    10
#>  [576,]    1    6   12    7    8    9    2    4    5    10    11     3
#>  [577,]    8    9    5   11    6    2    7    4    3    12     1    10
#>  [578,]    1    3    5    9    6   10    8    2   11    12     7     4
#>  [579,]    7    4    8    5   10    2   12   11    1     3     9     6
#>  [580,]   10    1    8   11    4    3    7   12    2     6     5     9
#>  [581,]    8    5    1    4    6   12    3    2   10    11     7     9
#>  [582,]    9    1    4   12    7    3    2    5   11    10     6     8
#>  [583,]    8    2    4   12    7    6    9    3   10    11     5     1
#>  [584,]    5    8    7    1    9    6    4    3   12    11     2    10
#>  [585,]    6   10   12   11    4    9    7    8    2     5     3     1
#>  [586,]    1    3   12    2    5   11    6    8    7     9    10     4
#>  [587,]    1    2    3    5    6    8   10   12    4     7    11     9
#>  [588,]    3    4    7    9   12    8    2   11    1    10     6     5
#>  [589,]    7    8   12   10   11    5    3    6    2     1     9     4
#>  [590,]   11   10    1    2    8    3    5    4   12     7     6     9
#>  [591,]    4    6   12    3    9    8    1    7   10    11     5     2
#>  [592,]    7    8    6    9    4    2   11   12    1     5     3    10
#>  [593,]    3   11    7    4    1    6    2   10    9    12     8     5
#>  [594,]    4   11    3    8    2    5    6    7    9    10     1    12
#>  [595,]   11    5    8   10    6    7   12    4    2     9     3     1
#>  [596,]   11    9    6   10    3   12    4    7    5     8     1     2
#>  [597,]    5   12    1    8   10   11    2    7    6     9     4     3
#>  [598,]   12    7    9   10    4    5    2    3   11     6     1     8
#>  [599,]    9    1    2    6   10    8   11   12    3     4     7     5
#>  [600,]   10    4    3    6    5    9    1    7    2    12    11     8
#>  [601,]    5    8   11    4    1   10    9    2   12     7     3     6
#>  [602,]    5   12   11    3    2    7    1    9    8    10     4     6
#>  [603,]   12    3    1    6   11    2    5    4    8    10     7     9
#>  [604,]    2   12   11    8    7    6    1    5    4     3    10     9
#>  [605,]    2    3    8    5   11    6    9   10    7    12     1     4
#>  [606,]    9    2   10   12    3    6    5    8    4     7    11     1
#>  [607,]   12    2    3    8    5    7   11   10    4     1     9     6
#>  [608,]   10   11    8    3    7    9    2    4    1    12     5     6
#>  [609,]   12   10    4    8    5    9   11    6    7     1     2     3
#>  [610,]   12    5    6    7    4   11    1   10    3     9     8     2
#>  [611,]    1   11    4    7    3   10    2   12    5     9     8     6
#>  [612,]    1    5   11    4   12   10    8    3    2     9     6     7
#>  [613,]    4    6   11    5    9    8    7   12    1    10     2     3
#>  [614,]    8    5   12   11   10    2    9    6    3     1     7     4
#>  [615,]    4    3   11    7    2   12    6    5    1    10     9     8
#>  [616,]   10    9    4   11    6    3    8    2    1    12     5     7
#>  [617,]    9   11   12    8    2    1    5   10    6     3     7     4
#>  [618,]    6    5    9    3    4    1    7   11   10    12     8     2
#>  [619,]    2   10   11    5    4    6    7    8   12     9     1     3
#>  [620,]    5    6   11    9    2   12   10    7    1     3     8     4
#>  [621,]    1    6    3   10   12    9   11    8    7     2     4     5
#>  [622,]    5    9    4   11   12    2    8    6    7     3     1    10
#>  [623,]    4    6    5    8    2   12    9    1    3     7    11    10
#>  [624,]    1    5    3    6   10   11    7    8   12     2     4     9
#>  [625,]    3    4   10    1    2    8    7    6    5    12     9    11
#>  [626,]    5   12    3    9   10    2    1    7    8    11     6     4
#>  [627,]   10    6   12   11    7    3    2    5    8     9     4     1
#>  [628,]    3    1    6   10    8    9   11   12    7     2     5     4
#>  [629,]    1    5    2    7    3   10    6   11    9     8    12     4
#>  [630,]   10    3    4   12    5    2    8    1    7    11     6     9
#>  [631,]    2   11   12    3    8    9   10    6    4     5     1     7
#>  [632,]   12   11    1   10    5    7    3    2    6     9     4     8
#>  [633,]   10    2   12    8   11    7    6    9    1     4     5     3
#>  [634,]    8    1    3    2    9   12    4   10   11     7     6     5
#>  [635,]    8    3    2   10    5    9   12    4    1     6    11     7
#>  [636,]    2   10    1    4   11    9   12    3    5     8     6     7
#>  [637,]    7   10    2    3    5   12    8    9    4     6     1    11
#>  [638,]    7   12    5   10    8    2    4    1    6     3     9    11
#>  [639,]    5    6   10    1    8   12   11    4    9     7     2     3
#>  [640,]    9    2    6   12    1    8   10   11    7     4     5     3
#>  [641,]    6    4    7    5   10   11    2    1   12     8     3     9
#>  [642,]    4    2    8    6    1   10    5   12   11     3     7     9
#>  [643,]    6    3   12   11    8    4    7    5   10     9     2     1
#>  [644,]    8   10   11    3    6    5    7    2    9    12     1     4
#>  [645,]    2   12    7   11    1    8   10    5    3     9     6     4
#>  [646,]    5    8   12    4   10    9    7    1    6    11     2     3
#>  [647,]    7   11   12    9    2    1    8    4    6     3    10     5
#>  [648,]   11    6    4    9    5    7    3    1   12     8    10     2
#>  [649,]    1    8    9   11    7   10    3    5   12     6     2     4
#>  [650,]    9    1   11    6   12    3   10    8    7     2     5     4
#>  [651,]    1   12    5    4   11    2    8    3   10     6     7     9
#>  [652,]    4    6    9    2    8    7    5   11   10     1    12     3
#>  [653,]   10    1    6    5    9    7    3    4   11     2    12     8
#>  [654,]    1    7   11   10    6    8    5    9    3     2    12     4
#>  [655,]   11    1    8    4    2   12   10    9    6     7     3     5
#>  [656,]   10    1    5    7    3    8    4    9   12     6    11     2
#>  [657,]    5    3   12   11    7    8    6    1   10     4     2     9
#>  [658,]   10    7    9    6   12    5    1    4    8     2    11     3
#>  [659,]    2   11    9    7    5    1    8   12    6    10     4     3
#>  [660,]    9    4   11    5   12   10    2    6    1     8     7     3
#>  [661,]   11    8    5    7    2   12    4    1   10     9     3     6
#>  [662,]   11    8    4    7    9   10    3    5    6     1    12     2
#>  [663,]    5    2   10    3    1   12   11    6    4     7     8     9
#>  [664,]    1    7    2   10    9   12    3   11    5     8     4     6
#>  [665,]    4    1    2    8   10    9    5    6   11    12     7     3
#>  [666,]    1    3    4    2   10    5   12    8    6     9    11     7
#>  [667,]    8    4    9    1   11    7   10    2   12     5     3     6
#>  [668,]    3    8   11    7   12    2    6    9    5     4    10     1
#>  [669,]   10   11    2    5    3   12    8    4    7     9     6     1
#>  [670,]   11   12    3    7    8    6    2   10    5     9     4     1
#>  [671,]    4   10   12    2    7   11    6    8    1     3     5     9
#>  [672,]    5    4    2    6   10   12    8    3   11     1     7     9
#>  [673,]   10   12    3   11    1    7    9    6    5     4     8     2
#>  [674,]   12    3    4    6    7   11    2    1    8    10     5     9
#>  [675,]    7    9   11   12    1    3    2    6    8    10     4     5
#>  [676,]   11   12    4   10    2    7    3    1    6     5     9     8
#>  [677,]   12    5    8    7    9    4    2    3   10     1    11     6
#>  [678,]    5   12    6   11   10    2    1    9    3     7     8     4
#>  [679,]    6   10    9   11    8    2    1    7    4     3     5    12
#>  [680,]   12    9    7    4    8    5    2   10   11     1     6     3
#>  [681,]    8    7   12    2    1   10    9   11    4     5     3     6
#>  [682,]    2    8   10   12    5    9   11    1    4     6     7     3
#>  [683,]   10    8    7    3    9    6    5    1    4    12     2    11
#>  [684,]    3    9   12    5    6   10    7    4    8     1    11     2
#>  [685,]    4   12    5    3    9    7    8   10    1     6    11     2
#>  [686,]    1    3   11    4    5    7   12    2    8    10     9     6
#>  [687,]   11    8   10    3   12    5    9    4    7     1     6     2
#>  [688,]    3    5    1    7   12    9    8    4    2    10    11     6
#>  [689,]    9    7   12    4    1    8   11    2   10     3     6     5
#>  [690,]    7    1    6   11    8   12    3   10    2     9     4     5
#>  [691,]    9   10    6   12    2    1   11    8    7     3     4     5
#>  [692,]    6    8    3    5    7   11    1    2    4    12     9    10
#>  [693,]   10    3    8    6    9   11    5    7    1    12     2     4
#>  [694,]   10    7    2   12    5    9    6    8    1    11     3     4
#>  [695,]    9    4    7    5    2    6    1   12    8     3    10    11
#>  [696,]    3    7   10    8    5    6    1   11    2     4    12     9
#>  [697,]   12    9    1   11    4    5    6    3    8     7    10     2
#>  [698,]   10    2    1    9    3    5    8    4    7    12     6    11
#>  [699,]   11   10   12    7    2    9    5    6    3     8     1     4
#>  [700,]    3    7    8   12   11    9    2    4    5    10     6     1
#>  [701,]   11    7    6    9    1    2    3    8   12    10     5     4
#>  [702,]    9   11    6    4    1   12    7    8    5    10     3     2
#>  [703,]    7    8   10    3    6   12   11    5    2     4     9     1
#>  [704,]    9    6    5   10   12    1    7    8    3     2     4    11
#>  [705,]    7   10    1    5    9   12   11    4    6     3     2     8
#>  [706,]    3   11    7    1    8    2   12    5    9     4    10     6
#>  [707,]    4    2    3   12    5    6    9   10   11     1     8     7
#>  [708,]    4    5    1    7    3   10    9    8    2    11    12     6
#>  [709,]    1   11    2    5    7    8   12    6    9    10     3     4
#>  [710,]    4    2    9    5   10    7    3   11    8    12     6     1
#>  [711,]    6    1   10    3    2    8    7    9   12     5     4    11
#>  [712,]    2    3   11    6    9   10    4    5    7     1    12     8
#>  [713,]    2   10    9    5   12    7    3    6    4     1    11     8
#>  [714,]    2   10    7   11    8    4    9    3   12     5     1     6
#>  [715,]    8    5    1   12    9    4    3    6   11     7    10     2
#>  [716,]    8    5   10    4    9    2    1    6   11     3    12     7
#>  [717,]   11    5    7   12    3    2    1    8    9     6    10     4
#>  [718,]   11    4    6    2   10   12    9    8    1     5     7     3
#>  [719,]    4    9   11    7    3   12    1    2   10     8     6     5
#>  [720,]    5   10    6    9    1    2    3   12    7    11     8     4
#>  [721,]   10    1    9    2    3    5    8    4    7    12    11     6
#>  [722,]    5   11   12    6    7    8    2    3   10     9     1     4
#>  [723,]    1   12    5    6    7    2   11    4    8     3    10     9
#>  [724,]    9    6   12    5   11    8    3    7    2    10     1     4
#>  [725,]    8    9    4   11   10    6    7   12    1     5     3     2
#>  [726,]    6    9    8    7    1    2   12   11    3     4    10     5
#>  [727,]    3    2   12   10    7    9    4   11    6     5     8     1
#>  [728,]   11    1   10    6   12    9    7    8    2     3     4     5
#>  [729,]   12    6    4    2    5    8    7    3    9    11    10     1
#>  [730,]    8    9    3    4    2    7    6   12   10     1     5    11
#>  [731,]   12   10    6    7    9    1    2    5    3    11     4     8
#>  [732,]    8    1    6   12   11    4    7    9    2    10     3     5
#>  [733,]    8    6    7   10    5    1    4    3    2    11    12     9
#>  [734,]   10    6    2    8    5    1    4    7   12     9    11     3
#>  [735,]    3    7    1    6   10   11    8   12    2     9     4     5
#>  [736,]    3    9   11    7    1    5    2   10    8     4    12     6
#>  [737,]   10    9    5    6   11    4    3    1    2    12     7     8
#>  [738,]   10    4    9    3    1   11    5    2    8     7    12     6
#>  [739,]    8    2   10    3    9   12    5    7    1     6     4    11
#>  [740,]    2    6    8   11    1    3    9    5   10     4     7    12
#>  [741,]    9    5    4    2   12   11    6    3    7    10     1     8
#>  [742,]    9    8   11   12    5    1   10    3    7     6     4     2
#>  [743,]    5   12    6    2    3    8    7    9   11     1     4    10
#>  [744,]    8    5    6    2    9    4    1   10   11     3     7    12
#>  [745,]    3   11    2    1   10    7   12    8    6     4     9     5
#>  [746,]    2   12    8    9    7   10    5    4    3    11     1     6
#>  [747,]    3    9    7   10    6   12    8    1   11     4     2     5
#>  [748,]   11    8    3   12    2    9    6    5    7     1     4    10
#>  [749,]   11    8    6    4    5    7   10    2    3     9     1    12
#>  [750,]    6    1    3   11    7    4    9    2    5     8    10    12
#>  [751,]    7    5   11    4   12    9    6    8    1     3    10     2
#>  [752,]    2   12    7    8    9    6    1    3   11    10     4     5
#>  [753,]    3    4   10    2    8   12    7   11    6     1     5     9
#>  [754,]    9    1    8    6    5   12    7    4   10     3     2    11
#>  [755,]    6    3    7    2   11    5    4    8    1    12    10     9
#>  [756,]    9    4    8    5   12   11    6    1    2     3    10     7
#>  [757,]    7    9    8    3   12    5    2   11    4     6    10     1
#>  [758,]   11    5    8    3    7    1   12    6    2    10     9     4
#>  [759,]    7    2    3    9    1    5   10   12    6    11     8     4
#>  [760,]    2    9   11    5    7    1    3    6   10     4     8    12
#>  [761,]    6    9   12   10    7    5    1    3    2     4     8    11
#>  [762,]    2    7    6    4    8    5   11   10    3    12     1     9
#>  [763,]    6    3    1    9    2   12    4    7    8    11    10     5
#>  [764,]   12    4    5   11    3    7    8   10    2     9     1     6
#>  [765,]    5   10    7    8   11    6    9    2    4     3     1    12
#>  [766,]    4    7   11    5    2    6   12    9    3    10     1     8
#>  [767,]    1   10    2    8    4   12   11    3    9     6     7     5
#>  [768,]    4    3    5    1    2   11    8   12   10     7     6     9
#>  [769,]   10    6   11    1    9   12    2    4    8     3     5     7
#>  [770,]    7    5    4    1    9   12    3    2   11     8    10     6
#>  [771,]    4   12    7    5    9    8   11    2   10     1     6     3
#>  [772,]    9    8    1    4    6    3    5   11    2    10    12     7
#>  [773,]    1   11   12    6   10    5    4    9    2     8     7     3
#>  [774,]   12   11    5    8    2    7    9    6   10     4     1     3
#>  [775,]    8    9   11    2    6    3   12    1   10     4     5     7
#>  [776,]   10    7    5   12    1   11    9    8    2     6     4     3
#>  [777,]    4   10    7    1    8    9    5   12   11     2     3     6
#>  [778,]    2    6    7    5    1   11   12   10    9     4     8     3
#>  [779,]   11    4    2    7    5   12   10    6    8     3     1     9
#>  [780,]    6    5    9    2    3    4   11    8    1    12     7    10
#>  [781,]    4    7    3   10    5   11    2    9   12     8     1     6
#>  [782,]    8    9    7    3    6    5    4   11    2    10     1    12
#>  [783,]    4    8    1    2    3    5   12   10    7     9    11     6
#>  [784,]    8    7    5   10    3    1    4   11    6    12     2     9
#>  [785,]    1   11    5   12    4    7    9    6   10     2     8     3
#>  [786,]    8   10    4    3    1    9    2    7    5    12    11     6
#>  [787,]   11    3    5    9    8    4   12    6   10     2     1     7
#>  [788,]    2    7   12    1    8    5    4   10    6    11     3     9
#>  [789,]    3    6    7    4    9   11    2    8    1    10     5    12
#>  [790,]    8    7   10    9    4    5   11    6    3     2     1    12
#>  [791,]    7    8    1    4   10    2    3   12    5     6     9    11
#>  [792,]   12    4    8    2    5    6    9   10    3    11     7     1
#>  [793,]   10    6    8    3   12   11    7    4    2     5     9     1
#>  [794,]    5   12    8   11    4    9    3   10    7     2     6     1
#>  [795,]    2   11    6    4   12    3    8    9    7     1    10     5
#>  [796,]   11    9   10    6    1    3    4    7    5     8     2    12
#>  [797,]    3    9    2   11    7    5   12    4   10     1     6     8
#>  [798,]   11    8    4   10   12    3    2    1    6     9     7     5
#>  [799,]    8   10    2   11    1    3   12    9    5     4     6     7
#>  [800,]    5    8   11   10    9    1    4    3   12     6     7     2
#>  [801,]   12   10    7    2    9    5    3    8    1    11     6     4
#>  [802,]    4    2   12    7    5    3   10    8   11     6     9     1
#>  [803,]   10    8   12    4    7    9   11    3    2     1     5     6
#>  [804,]    5    2    7    8   12    1    3    9   10     6    11     4
#>  [805,]   11    2    1    6   10    9    4    5    7     3    12     8
#>  [806,]    9   11   10    3    8    4    7    5    6     1     2    12
#>  [807,]    4   12    5    3    8    6   10    9    1     7    11     2
#>  [808,]    3    7    8   12    5   10   11    2    1     4     6     9
#>  [809,]    7   10   12    9    1    3   11    2    4     8     6     5
#>  [810,]   11    9    4    5   12    6    7    3    1     2     8    10
#>  [811,]    3    1   12    2    8   11   10    9    5     6     4     7
#>  [812,]    5    9    8   10   12    1    4    3   11     2     6     7
#>  [813,]   10   11    2    9    1    7    8   12    6     4     5     3
#>  [814,]    2    4    7    1    3    8    9   10    5    12    11     6
#>  [815,]    2   12    9    1    4    5    6    3   10    11     8     7
#>  [816,]    9    5   11    1    3    4    6    2    7    10    12     8
#>  [817,]    8   10   12    7    5    3    9    2    4     6     1    11
#>  [818,]    9    6   11    1    2    4    8    5   10     7    12     3
#>  [819,]    3    9    2    4    5   11   10   12    6     7     8     1
#>  [820,]    3    8   12   10   11    6    4    1    7     5     2     9
#>  [821,]    3   10   11    4    1    8    6    5    7     2    12     9
#>  [822,]    8    6    2    7   11    4   10    3    1    12     5     9
#>  [823,]    3    5    6   10    8    2   11   12    7     4     9     1
#>  [824,]   11    9    7    6    2    5   12    3   10     4     8     1
#>  [825,]    8    6    5    9   11    1   12    7   10     4     2     3
#>  [826,]   11    2   12    3    1   10    5    8    9     7     6     4
#>  [827,]    8    3   11   12    6   10    2    4    5     7     9     1
#>  [828,]    1    9   12   11    5    8    6   10    4     2     3     7
#>  [829,]    6    1    7   12    8   10    3    2    4     5    11     9
#>  [830,]    8    2    9    3    6    5    4   11    7    10     1    12
#>  [831,]   11   12    9    4    6    8    3    7    1    10     2     5
#>  [832,]    7   12    6    5    1    4    3    8    9    11     2    10
#>  [833,]    7    1    5    9    3    2   10   12   11     6     8     4
#>  [834,]   12    5    4    3    8    6   10    1    2     9    11     7
#>  [835,]    4    1    3    9    6    5    7   12   10    11     8     2
#>  [836,]    9   11    1    5    7    8   10    6    3    12     2     4
#>  [837,]    8   11    7   12   10    2    5    6    4     1     3     9
#>  [838,]    8   11    2    7    5    3    9    1   12    10     6     4
#>  [839,]    9    2   12    3    6    8    4   10    5    11     7     1
#>  [840,]    2    1    8   11    5    3    7    6    9    12     4    10
#>  [841,]    6   10    4    3   11    1    2    5    8     7    12     9
#>  [842,]   12   11    2    6    9    1    3   10    4     5     8     7
#>  [843,]    8   11    4    3    1    2    5   12    9    10     7     6
#>  [844,]    7    6    5    3    8    4    1   11    2     9    10    12
#>  [845,]   10    9    8   12    2    1    4    7    5     6    11     3
#>  [846,]   10    7   12    3    5    9    8    6    1     4     2    11
#>  [847,]    3   10    4   11    5    8    1   12    6     2     7     9
#>  [848,]    1    2    7    5    8   12    3    6    4    10     9    11
#>  [849,]   10    6    9    7    5    2    1    3    4     8    12    11
#>  [850,]    1    6   12    5    3   10   11    2    9     8     4     7
#>  [851,]    9    8    5    7    3   12    2    6    4     1    11    10
#>  [852,]    8   12    9   10    3    5    4    6   11     7     1     2
#>  [853,]    2    7    8    3    9   12   10    6    5     4     1    11
#>  [854,]    9   10    5    6    2   11    3    8    7     1     4    12
#>  [855,]    2   11   12    3    5   10    4    7    8     6     1     9
#>  [856,]   11    5    6    9    3    2    4   12    8     1     7    10
#>  [857,]    3    1   11   10    7   12    5    9    8     6     4     2
#>  [858,]   10    8    5    2    1   12    4    7    9     3    11     6
#>  [859,]    3    8    4    9   10   11    7    2    1    12     5     6
#>  [860,]    3    9   10    7    5   11    8    4    6    12     2     1
#>  [861,]    5    1   10    7    3    6    8    2   11     4    12     9
#>  [862,]    3    9    6   12   11    8    2   10    7     4     1     5
#>  [863,]   11    4    7    3    1   12    9    8   10     2     5     6
#>  [864,]    9   12    3    4    1    5   10    8    7     2     6    11
#>  [865,]    2    1    4    5   10    7    3   11    8    12     6     9
#>  [866,]   11    3    8    2   10    9    7    4    1     6     5    12
#>  [867,]    2    7    4    8    9   10    3    6    5     1    12    11
#>  [868,]   11    6    2    9   12    1    7   10    5     3     4     8
#>  [869,]    1    9    6   10    8   11    7    2    5     4    12     3
#>  [870,]   11   10    5    2    3    9    6    4    7     8    12     1
#>  [871,]    8   12    9    5    1    7   11   10    4     2     6     3
#>  [872,]    6   10    1    9   12    3    8    5    2     7     4    11
#>  [873,]   11    2    6    3   12    4    5    7    8     1    10     9
#>  [874,]    4    9   11    8    5    6    2    1    7     3    12    10
#>  [875,]    2   11    7    4    8   12    1    9    5    10     3     6
#>  [876,]    1   12   10    8    5    4    6    9   11     7     3     2
#>  [877,]    8   12    5    4    7    9    1    2    6    10     3    11
#>  [878,]    9    7    6   11   12    5   10    4    2     1     8     3
#>  [879,]   10    4   12    7    6    8    5    3   11     2     1     9
#>  [880,]    5    4    9    8   11    6    1    2   12     3     7    10
#>  [881,]    9    4    2    6   12    7   11    5    8    10     3     1
#>  [882,]    1    9   10    3    7    5    6    4   11     2    12     8
#>  [883,]   10    2    8    3    6    7    1   11   12     9     4     5
#>  [884,]    7    9   10    3    6    4    2    8   11    12     1     5
#>  [885,]    3    8   10    9    4    5   12    1    6     7     2    11
#>  [886,]    7   11    1    2    9   12    6    5    4    10     8     3
#>  [887,]    2   12   11    3    4    9   10    7    1     5     8     6
#>  [888,]    3    7    1    4    6   10   11    5   12     8     9     2
#>  [889,]    3   11   10    5    7    9   12    8    2     4     6     1
#>  [890,]    2   11    9   12    4   10    1    3    6     5     7     8
#>  [891,]    9   12    5   11    7    2    6   10    8     4     1     3
#>  [892,]    9    5   11   12    6   10    2    7    4     8     1     3
#>  [893,]   12    2    8    7    9   10    6    4   11     5     1     3
#>  [894,]   10    6   11    2    4    8    3    1    7     9     5    12
#>  [895,]    3    4    8    9   11    1   10   12    2     6     7     5
#>  [896,]    9    6   11   10   12    8    4    3    1     2     7     5
#>  [897,]    7    4    9    1    2   11   10   12    5     8     6     3
#>  [898,]    4    6    2    7    9    1   12    8    3     5    11    10
#>  [899,]    9    2    7   12    6    5   11   10    3     4     8     1
#>  [900,]    7    9    8    3   12    1    4   10    5     2    11     6
#>  [901,]    2    9    3   11    5    1    7    6   12     4     8    10
#>  [902,]    5    1    6    2    8    4    3    7   11    10     9    12
#>  [903,]    5    8   12    2   11    6    1    3    4     7     9    10
#>  [904,]   12    8    7   10    3    4    5   11    9     2     6     1
#>  [905,]   10    9    8    6    7    2   12   11    3     1     5     4
#>  [906,]    5    2   11    1    4    3   12    6   10     7     9     8
#>  [907,]    8    5   10    2   12    3   11    9    1     6     4     7
#>  [908,]    4    9   12    3    6   11    8    7    5     2     1    10
#>  [909,]    9    1    4    2   11    7   12    8    5     3     6    10
#>  [910,]    8    1   12    9    5    3   11    6   10     2     4     7
#>  [911,]    1   11    7    9    3   10   12    2    8     6     5     4
#>  [912,]    6    2    4    9    8   10    7    1    3     5    11    12
#>  [913,]    7    5    8    4   10    9    2   12    3     6     1    11
#>  [914,]    6   12    1    8    9    3    4   11    2     5     7    10
#>  [915,]    1   11    8   12    9    7    2    6    5     4     3    10
#>  [916,]    6    8    3    9    2    1    4   10   12    11     7     5
#>  [917,]    7    8    1   10   11    6   12    9    4     5     2     3
#>  [918,]    4    3    5   11   12    7    9    6    2    10     1     8
#>  [919,]    7    1   12    3   10    8   11    9    6     2     4     5
#>  [920,]    4   10    6    7    3    8    1    2   11     5    12     9
#>  [921,]    5    7   10    8    1    6    2    3   11    12     9     4
#>  [922,]    9    6   10    7    4    3    8    2    5    11     1    12
#>  [923,]    8    9   10    1    3    6   12    2    4     5    11     7
#>  [924,]    9    8    7    6    4   12    5    2    3    11    10     1
#>  [925,]    9   11    8    7    6    3   12    1    2     5    10     4
#>  [926,]    1    6    3   10    9    8   11    5   12     4     2     7
#>  [927,]    4    6   10    7    1   11    5    8    2     9    12     3
#>  [928,]    5    7    3    8   11    9   10    6   12     4     2     1
#>  [929,]    1    3   10    6    4   12    2    8    9     5    11     7
#>  [930,]   10    8    1    9    4   11    7    6    2    12     5     3
#>  [931,]    5   10    1    3    7    2    4   11    6     9    12     8
#>  [932,]    9    7    4    6    3   10    2    8    5    11     1    12
#>  [933,]   10    2    4    3    8    9    6    7    1    11    12     5
#>  [934,]   11    1    4   10    5    2    8    6    9    12     7     3
#>  [935,]    1    6    9    5    2   12    4    8   11     7     3    10
#>  [936,]    7    3    8    5    9   12    1    2    4    10    11     6
#>  [937,]   12    3    2    5    6   10    4    7   11     8     1     9
#>  [938,]    6    2    9    1    4    3   11    5    8    10     7    12
#>  [939,]    5   10    4    9    8    7    3   11    6     1    12     2
#>  [940,]    8    5    6    3   11    2    4    9    1     7    12    10
#>  [941,]    8   12    1    7    3    2    6    4   10     5    11     9
#>  [942,]    5    4   10    1    2    3    8   11    6     7    12     9
#>  [943,]    8   10    4    2    1    7   11    6    5     3    12     9
#>  [944,]    4   11    6    5   10    9    7    3    2    12     8     1
#>  [945,]   11    9    1    4   10   12    6    5    8     3     2     7
#>  [946,]    3    5   11    8    9    2    4    1    6    12     7    10
#>  [947,]    4    3    5    9    1    7    2    6    8    10    12    11
#>  [948,]    8    9    2    3   11    6    1    5    4     7    12    10
#>  [949,]    3   12    5   10    6    8    2   11    9     4     7     1
#>  [950,]    2    5    8    6   12    4    7    9   11     1    10     3
#>  [951,]    2   11    3   12    6    1   10    4    5     7     9     8
#>  [952,]    7    4    1    9    6    3    8   12   10     5    11     2
#>  [953,]    7   10   12    9    6    2   11    5    1     3     8     4
#>  [954,]    6    7   10    8    4    3   11    1   12     9     5     2
#>  [955,]   12    1   11    9    4    6    5    7    8     2    10     3
#>  [956,]    4    8    6    3    5    1   10    9   11     2    12     7
#>  [957,]    5    2    7    6    8    3    1   10   12     4     9    11
#>  [958,]    4    3    5    7   12   11    1    8    6    10     9     2
#>  [959,]   10    6    9   11    4    5    7    1    3     2     8    12
#>  [960,]   12    4    6   11    3    8    7    5   10     2     1     9
#>  [961,]   10    1    5    2    4    9    8   11    3     6    12     7
#>  [962,]    7    5    6    1    4    2   12   11    9     3    10     8
#>  [963,]   12    7    4    9   10    1    2   11    8     3     5     6
#>  [964,]    7    9   11    4    2   12    1    5   10     6     3     8
#>  [965,]   10    7    2   12    5    3    6   11    9     4     1     8
#>  [966,]    3   11   10    9    5    8    4    6    2     1     7    12
#>  [967,]    3    4    8    5   10    9    6   11   12     1     2     7
#>  [968,]    8    9    3   12    6    1   10    5    4     2     7    11
#>  [969,]    6    9    8    1    5    2   12   10    4    11     3     7
#>  [970,]    1    9    2    4    7    5   12    3    6    11     8    10
#>  [971,]    8    6    4    9   12    7    2    1    5     3    11    10
#>  [972,]    5    7    9    6   11   12    3   10    4     2     8     1
#>  [973,]    3    7    6   12    9    5    8   10    2     1     4    11
#>  [974,]    2   10    4    8   11    5   12    1    9     3     7     6
#>  [975,]   12    3   11   10    9    7    4    5    8     2     1     6
#>  [976,]    8    9    4   12    2    3    6   10    1     5     7    11
#>  [977,]    9    1    5    2    8   12    3   10   11     6     4     7
#>  [978,]    4    2    5    1    9   10   12   11    6     7     3     8
#>  [979,]    9   11    2    1    6    4   12    5    7    10     8     3
#>  [980,]    8    3   10    7   11    9    6    2   12     1     5     4
#>  [981,]    6   10    1    3    4   12   11    7    2     8     9     5
#>  [982,]    4   10    7    5    9    1    6    8   12    11     2     3
#>  [983,]   10    1    2    5    9    6    7    3    8    12    11     4
#>  [984,]   10   12    8    5    6   11    1    9    7     4     3     2
#>  [985,]   10    4   12    6   11    5    9    8    3     1     7     2
#>  [986,]    7    3   12    4   11    9   10    8    5     1     6     2
#>  [987,]    8    1    6    7    9   12    5   10    3    11     4     2
#>  [988,]    2    5   11    8    3   12    4    7    1     6    10     9
#>  [989,]   12   11    5    3    1    2    8    9    7     6    10     4
#>  [990,]    4    1    3    2   11   10    5    7   12     8     6     9
#>  [991,]   12    4    5    6    1    2   11   10    9     3     8     7
#>  [992,]    3    1    6    9    7   12   10    5    4    11     8     2
#>  [993,]    6    4   11    2    3    5    7    8    9    12     1    10
#>  [994,]    8    2    9    4    6   12    5    3   10     1    11     7
#>  [995,]    4    8    7    3   12   11   10    9    1     6     2     5
#>  [996,]    5    8    2    3    7   10   12   11    4     6     9     1
#>  [997,]    3   12    9   10    8    4    7    2    1    11     5     6
#>  [998,]   11   12    1   10    9    5    7    4    8     3     2     6
#>  [999,]   10    8    4    3    2   12    9    5    7     1    11     6
#> [1000,]    4   10    9    8    2   11    6    5    3     7     1    12


# Graph of bootstrap distributions
boxplot(as.vector(Cornell.bootYX$t[,-1])~factor(rep(1:7,rep(1000,7))),
main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.bootYX$t0[-1],col="red",pch=19)

# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.bootYX,indices=2:8)



library(boot)
plot(Cornell.bootYX,index=2)


qqnorm(Cornell.bootYX$t[,2],ylim=c(-1,1))
abline(h=Cornell.bootYX$t0[2],lty=2)

(sum(abs(Cornell.bootYX$t[,2])>=abs(Cornell.bootYX$t0[2]))+1)/(length(Cornell.bootYX$t[,2])+1)
#> [1] 0.5584416

rm(Cornell.bootYX)
# }
```
