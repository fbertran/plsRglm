# Correlation matrix for simulating plsR datasets

A correlation matrix to simulate datasets

## Format

A data frame with 17 observations on the following 17 variables.

- y:

  a numeric vector

- x11:

  a numeric vector

- x12:

  a numeric vector

- x13:

  a numeric vector

- x21:

  a numeric vector

- x22:

  a numeric vector

- x31:

  a numeric vector

- x32:

  a numeric vector

- x33:

  a numeric vector

- x34:

  a numeric vector

- x41:

  a numeric vector

- x42:

  a numeric vector

- x51:

  a numeric vector

- x61:

  a numeric vector

- x62:

  a numeric vector

- x63:

  a numeric vector

- x64:

  a numeric vector

## Source

Handmade.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## Examples

``` r
data(CorMat)
str(CorMat)
#> 'data.frame':    17 obs. of  17 variables:
#>  $ y  : num  1 0.9 0.88 0.92 0.77 0.8 0.65 0.7 0.66 0.6 ...
#>  $ x11: num  0.9 1 0.75 0.8 0.1 0.1 0.05 0.1 0.1 0.05 ...
#>  $ x12: num  0.88 0.75 1 0.65 0.1 0.05 0.1 0.1 0.05 0.1 ...
#>  $ x13: num  0.92 0.8 0.65 1 0.05 0.1 0.15 0.05 0.1 0.15 ...
#>  $ x21: num  0.77 0.1 0.1 0.05 1 0.95 0.1 0.05 0.1 0.1 ...
#>  $ x22: num  0.8 0.1 0.05 0.1 0.95 1 0.05 0.1 0.15 0.05 ...
#>  $ x31: num  0.65 0.05 0.1 0.15 0.1 0.05 1 0.75 0.8 0.92 ...
#>  $ x32: num  0.7 0.1 0.1 0.05 0.05 0.1 0.75 1 0.65 0.55 ...
#>  $ x33: num  0.66 0.1 0.05 0.1 0.1 0.15 0.8 0.65 1 0.7 ...
#>  $ x34: num  0.6 0.05 0.1 0.15 0.1 0.05 0.92 0.55 0.7 1 ...
#>  $ x41: num  0.2 0.1 0.1 0.05 0.1 0.1 0.1 0.1 0.05 0.1 ...
#>  $ x42: num  0.15 0.1 0.05 0.1 0.1 0.05 0.1 0.05 0.1 0.1 ...
#>  $ x51: num  0.1 0.05 0.1 0.15 0.05 0.1 0.05 0.1 0.1 0.1 ...
#>  $ x61: num  0.05 0.1 0.1 0.05 0.1 0.1 0.1 0.1 0.1 0.1 ...
#>  $ x62: num  0.07 0.1 0.05 0.1 0.1 0.05 0.1 0.05 0.05 0.1 ...
#>  $ x63: num  0.08 0.05 0.1 0.15 0.05 0.1 0.05 0.1 0.1 0.05 ...
#>  $ x64: num  0.02 0.15 0.1 0.05 0.1 0.1 0.1 0.1 0.05 0.1 ...
```
