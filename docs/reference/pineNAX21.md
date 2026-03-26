# Incomplete dataset from the pine caterpillars example

The caterpillar dataset was extracted from a 1973 study on pine
processionary caterpillars. It assesses the influence of some forest
settlement characteristics on the development of caterpillar colonies.
There are k=10 potentially explanatory variables defined on n=33
areas.  
The value of x2 for the first observation was removed from the matrix of
predictors on purpose.

## Format

A data frame with 33 observations on the following 11 variables and one
missing value.

- x1:

  altitude (in meters)

- x2:

  slope (en degrees)

- x3:

  number of pines in the area

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

Tomassone R., Audrain S., Lesquoy-de Turckeim E., Millier C. (1992). “La
régression, nouveaux regards sur une ancienne méthode statistique”,
INRA, *Actualités Scientifiques et Agronomiques*, Masson, Paris.

## Details

These caterpillars got their names from their habit of moving over the
ground in incredibly long head-to-tail processions when leaving their
nest to create a new colony.  
The `pineNAX21` is a dataset with a missing value for testing purpose.

## Examples

``` r
data(pineNAX21)
str(pineNAX21)
#> 'data.frame':    33 obs. of  11 variables:
#>  $ x1 : int  1200 1342 1231 1254 1357 1250 1422 1309 1127 1075 ...
#>  $ x2 : int  NA 28 28 28 32 27 37 46 24 34 ...
#>  $ x3 : int  1 8 5 18 7 1 22 7 2 9 ...
#>  $ x4 : num  4 4.4 2.4 3 3.7 4.4 3 5.7 3.5 4.3 ...
#>  $ x5 : num  14.8 18 7.8 9.2 10.7 14.8 8.1 19.6 12.6 12 ...
#>  $ x6 : num  1 1.5 1.3 2.3 1.4 1 2.7 1.5 1 1.6 ...
#>  $ x7 : num  1.1 1.5 1.6 1.7 1.7 1.7 1.9 1.3 1.7 1.8 ...
#>  $ x8 : num  5.9 6.4 4.3 6.9 6.6 5.8 8.3 7.8 4.9 6.8 ...
#>  $ x9 : num  1.4 1.7 1.5 2.3 1.8 1.3 2.5 1.8 1.5 2 ...
#>  $ x10: num  1.4 1.7 1.4 1.6 1.3 1.4 2 1.6 2 2 ...
#>  $ x11: num  2.37 1.47 1.13 0.85 0.24 1.49 0.3 0.07 3 1.21 ...
```
