# Graphical assessment of the stability of selected variables

This fonctions plots, for each of the model, the

## Usage

``` r
signpred(
  matbin,
  pred.lablength = max(sapply(rownames(matbin), nchar)),
  labsize = 1,
  plotsize = 12
)
```

## Arguments

- matbin:

  Matrix with 0 or 1 entries. Each row per predictor and a column for
  every model. 0 means the predictor is not significant in the model and
  1 that, on the contrary, it is significant.

- pred.lablength:

  Maximum length of the predictors labels. Defaults to full label
  length.

- labsize:

  Size of the predictors labels.

- plotsize:

  Global size of the graph.

## Value

A plot window.

## Details

This function is based on the
[`visweb`](https://rdrr.io/pkg/bipartite/man/visweb.html) function from
the bipartite package.

## References

Vazquez, P.D., Chacoff, N.,P. and Cagnolo, L. (2009) Evaluating multiple
determinants of the structure of plant-animal mutualistic networks.
*Ecology*, 90:2039-2046.

## See also

See Also [`visweb`](https://rdrr.io/pkg/bipartite/man/visweb.html)

## Author

Bernd Gruber with minor modifications from Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
signpred(matrix(rbinom(160,1,.2),ncol=8,dimnames=list(as.character(1:20),as.character(1:8))))

```
