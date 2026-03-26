# Complete Pine dataset

This is a supplementary dataset (used as a test set for the `pine`
dataset) that was extracted from a 1973 study on pine_sup processionary
caterpillars. It assesses the influence of some forest settlement
characteristics on the development of caterpillar colonies. The response
variable is the logarithmic transform of the average number of nests of
caterpillars per tree in an area of 500 square meters (`x11`). There are
k=10 potentially explanatory variables defined on n=22 areas.

## Format

A data frame with 22 observations on the following 11 variables.

- x1:

  altitude (in meters)

- x2:

  slope (en degrees)

- x3:

  number of pine_sups in the area

- x4:

  height (in meters) of the tree sampled at the center of the area

- x5:

  diameter (in meters) of the tree sampled at the center of the area

- x6:

  index of the settlement density

- x7:

  orientation of the area (from 1 if southbound to 2 otherwise)

- x8:

  height (in meters) of the dominant tree

- x9:

  number of vegetation strata

- x10:

  mix settlement index (from 1 if not mixed to 2 if mixed)

- x11:

  logarithmic transform of the average number of nests of caterpillars
  per tree

## Source

Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992), “La
régression, nouveaux regards sur une ancienne méthode statistique”,
INRA, *Actualités Scientifiques et Agronomiques*, Masson, Paris.

## Details

These caterpillars got their names from their habit of moving over the
ground in incredibly long head-to-tail processions when leaving their
nest to create a new colony.  

The `pine_sup` dataset can be used as a test set to assess model
prediction error of a model trained on the `pine` dataset.

## References

J.-M. Marin, C. Robert. (2007). *Bayesian Core: A Practical Approach to
Computational Bayesian Statistics*. Springer, New-York, pages 48-49.

## Examples

``` r
data(pine_sup)
str(pine_sup)
#> 'data.frame':    25 obs. of  11 variables:
#>  $ x1 : int  1107 1116 1174 1131 1150 1132 1258 1114 1177 1146 ...
#>  $ x2 : int  31 34 32 30 34 22 14 26 36 26 ...
#>  $ x3 : int  23 6 22 6 12 18 0 9 3 18 ...
#>  $ x4 : num  6 2.5 3.9 4.7 3.1 7 3.8 3.4 4 4.3 ...
#>  $ x5 : num  22.2 6.6 11.9 15.3 9.4 37 8.5 9.6 14.2 17.9 ...
#>  $ x6 : num  2.6 1.3 2.3 1.5 1.7 2.5 1 1.6 1.2 2.3 ...
#>  $ x7 : num  1 1.8 1.7 1.5 1.8 1.5 1.2 1.6 1.3 1.6 ...
#>  $ x8 : num  9 3.9 6.1 6.5 4.8 9 5.6 5.1 5.9 7.7 ...
#>  $ x9 : num  3 1.2 1.8 1.4 1.6 2 1 1.5 1.3 2 ...
#>  $ x10: num  1.4 1.5 1.5 1.3 1.3 1.5 1 1.3 1.6 1.4 ...
#>  $ x11: num  1.17 0.67 0.9 2.32 3.89 6 3.18 0.9 2.5 2 ...
```
