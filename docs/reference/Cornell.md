# Cornell dataset

The famous Cornell dataset. A mixture experiment on `X1`, `X2`, `X3`,
`X4`, `X5`, `X6` and `X7` to analyse octane degree (`Y`) in gazoline.

## Format

A data frame with 12 observations on the following 8 variables.

- X1:

  a numeric vector

- X2:

  a numeric vector

- X3:

  a numeric vector

- X4:

  a numeric vector

- X5:

  a numeric vector

- X6:

  a numeric vector

- X7:

  a numeric vector

- Y:

  response value: a numeric vector

## Source

M. Tenenhaus. (1998). *La regression PLS, Theorie et pratique*. Editions
Technip, Paris.

## References

N. Kettaneh-Wold. Analysis of mixture data with partial least squares.
(1992). *Chemometrics and Intelligent Laboratory Systems*, 14(1):57-69.

## Examples

``` r
data(Cornell)
str(Cornell)
#> 'data.frame':    12 obs. of  8 variables:
#>  $ X1: num  0 0 0 0 0 0 0.17 0.17 0.17 0.17 ...
#>  $ X2: num  0.23 0.1 0 0.49 0 0.62 0.27 0.19 0.21 0.15 ...
#>  $ X3: num  0 0 0 0 0 0 0.1 0.1 0.1 0.1 ...
#>  $ X4: num  0 0 0.1 0 0.62 0 0.38 0.38 0.38 0.38 ...
#>  $ X5: num  0 0.12 0.12 0.12 0.12 0 0 0.02 0 0.02 ...
#>  $ X6: num  0.74 0.74 0.74 0.37 0.18 0.37 0 0.06 0.06 0.1 ...
#>  $ X7: num  0.03 0.04 0.04 0.02 0.08 0.01 0.08 0.08 0.08 0.08 ...
#>  $ Y : num  98.7 97.8 96.6 92 86.6 91.2 81.9 83.1 82.4 83.2 ...
```
