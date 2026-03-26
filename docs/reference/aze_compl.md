# As aze without missing values

This is a single imputation of the
[`aze`](https://fbertran.github.io/plsRglm/reference/aze.md) dataset
which was collected on patients carrying a colon adenocarcinoma. It has
104 observations on 33 binary qualitative explanatory variables and one
response variable `y` representing the cancer stage according to the to
Astler-Coller classification (Astler and Coller, 1954). A microsattelite
is a non-coding DNA sequence.

## Format

A data frame with 104 observations on the following 34 variables.

- y:

  the response: a binary vector (Astler-Coller score).

- D2S138:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D18S61:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D16S422:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D17S794:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D6S264:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D14S65:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D18S53:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D17S790:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D1S225:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D3S1282:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D9S179:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D5S430:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D8S283:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D11S916:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D2S159:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D16S408:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D5S346:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D10S191:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D13S173:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D6S275:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D15S127:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D1S305:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D4S394:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D20S107:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D1S197:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D1S207:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D10S192:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D3S1283:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D4S414:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D8S264:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D22S928:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- TP53:

  a binary vector that indicates whether this microsatellite is altered
  or not.

- D9S171:

  a binary vector that indicates whether this microsatellite is altered
  or not.

## Source

Weber *et al.* (2007). Allelotyping analyzes of synchronous primary and
metastasis CIN colon cancers identified different subtypes. *Int J
Cancer*, 120(3), pages 524-32.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Société
Française de Statistique*, 151(2), pages 1-18.

## Examples

``` r
data(aze_compl)
str(aze_compl)
#> 'data.frame':    104 obs. of  34 variables:
#>  $ y      : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ D2S138 : int  0 0 1 1 0 0 0 0 1 0 ...
#>  $ D18S61 : int  1 0 0 0 0 0 1 1 1 1 ...
#>  $ D16S422: int  1 1 0 0 1 0 0 1 1 0 ...
#>  $ D17S794: int  1 1 1 1 1 0 1 0 0 1 ...
#>  $ D6S264 : int  0 1 0 0 1 0 1 1 0 0 ...
#>  $ D14S65 : int  1 1 1 1 0 1 1 1 0 0 ...
#>  $ D18S53 : int  1 0 1 1 0 0 1 1 1 0 ...
#>  $ D17S790: int  0 1 0 0 1 0 1 0 0 1 ...
#>  $ D1S225 : int  0 1 0 0 0 0 1 0 1 0 ...
#>  $ D3S1282: int  0 0 1 0 0 0 1 1 0 0 ...
#>  $ D9S179 : int  1 1 1 0 0 0 0 0 1 1 ...
#>  $ D5S430 : int  1 1 0 1 1 0 1 0 1 0 ...
#>  $ D8S283 : int  1 1 0 1 0 0 0 1 1 0 ...
#>  $ D11S916: int  0 0 0 0 0 0 0 0 1 1 ...
#>  $ D2S159 : int  1 0 0 0 0 0 0 0 1 0 ...
#>  $ D16S408: int  1 1 0 1 0 0 1 1 1 0 ...
#>  $ D5S346 : int  0 1 0 1 0 0 0 1 1 0 ...
#>  $ D10S191: int  0 0 0 0 0 0 0 1 1 0 ...
#>  $ D13S173: int  1 1 1 1 1 0 0 1 0 1 ...
#>  $ D6S275 : int  0 1 1 0 0 1 1 0 1 1 ...
#>  $ D15S127: int  1 0 0 0 0 1 1 0 1 0 ...
#>  $ D1S305 : int  0 0 0 0 0 0 1 0 1 0 ...
#>  $ D4S394 : int  0 0 1 0 1 0 1 0 1 1 ...
#>  $ D20S107: int  1 1 1 1 0 0 1 1 1 1 ...
#>  $ D1S197 : int  0 0 1 0 0 1 1 1 1 0 ...
#>  $ D1S207 : int  0 0 0 0 0 0 1 0 1 0 ...
#>  $ D10S192: int  1 1 1 1 0 0 1 0 1 1 ...
#>  $ D3S1283: int  0 0 1 1 0 0 1 1 1 0 ...
#>  $ D4S414 : int  0 1 0 0 0 0 1 0 1 1 ...
#>  $ D8S264 : int  0 1 0 0 0 0 0 1 1 0 ...
#>  $ D22S928: int  0 1 0 0 0 0 0 1 1 0 ...
#>  $ TP53   : int  1 1 0 0 1 0 1 1 1 1 ...
#>  $ D9S171 : int  0 1 0 1 1 0 0 1 1 1 ...
```
