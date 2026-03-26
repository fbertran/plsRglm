# Quality of wine dataset

Quality of Bordeaux wines (`Quality`) and four potentially predictive
variables (`Temperature`, `Sunshine`, `Heat` and `Rain`).

## Format

A data frame with 34 observations on the following 5 variables.

- Temperature:

  a numeric vector

- Sunshine:

  a numeric vector

- Heat:

  a numeric vector

- Rain:

  a numeric vector

- Quality:

  an ordered factor with levels `1` \< `2` \< `3`

## Source

P. Bastien, V. Esposito-Vinzi, and M. Tenenhaus. (2005). PLS generalised
linear regression. *Computational Statistics & Data Analysis*,
48(1):17-46.

## References

M. Tenenhaus. (2005). La regression logistique PLS. In J.-J. Droesbeke,
M. Lejeune, and G. Saporta, editors, Modeles statistiques pour donnees
qualitatives. Editions Technip, Paris.

## Examples

``` r
data(bordeaux)
str(bordeaux)
#> 'data.frame':    34 obs. of  5 variables:
#>  $ Temperature: int  3064 3000 3155 3085 3245 3267 3080 2974 3038 3318 ...
#>  $ Sunshine   : int  1201 1053 1133 970 1258 1386 966 1189 1103 1310 ...
#>  $ Heat       : int  10 11 19 4 36 35 13 12 14 29 ...
#>  $ Rain       : int  361 338 393 467 294 225 417 488 677 427 ...
#>  $ Quality    : Ord.factor w/ 3 levels "1"<"2"<"3": 2 3 2 3 1 1 3 3 3 2 ...
```
