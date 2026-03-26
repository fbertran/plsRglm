# plsRglm: Partial Least Squares Regression for Generalized Linear Models

Provides (weighted) Partial least squares Regression for generalized
linear models and repeated k-fold cross-validation of such models using
various criteria
[doi:10.48550/arXiv.1810.01005](https://doi.org/10.48550/arXiv.1810.01005)
. It allows for missing data in the explanatory variables. Bootstrap
confidence intervals constructions are also available.

## References

A short paper that sums up some of features of the package is available
on <https://arxiv.org/>, Frédéric Bertrand and Myriam Maumy-Bertrand
(2018), "plsRglm: Partial least squares linear and generalized linear
regression for processing incomplete datasets by cross-validation and
bootstrap techniques with R", \*arxiv\*,
<https://arxiv.org/abs/1810.01005>,
<https://github.com/fbertran/plsRglm/> et
<https://fbertran.github.io/plsRglm/>

## See also

Useful links:

- <https://fbertran.github.io/plsRglm/>

- <https://github.com/fbertran/plsRglm>

- Report bugs at <https://github.com/fbertran/plsRglm/issues>

## Author

**Maintainer**: Frederic Bertrand <frederic.bertrand@lecnam.net>
([ORCID](https://orcid.org/0000-0002-0837-8281))

Authors:

- Myriam Maumy-Bertrand <myriam.maumy@ehesp.fr>
  ([ORCID](https://orcid.org/0000-0002-4615-1512))

## Examples

``` r
set.seed(314)
library(plsRglm)
data(Cornell)
cv.modpls<-cv.plsR(Y~.,data=Cornell,nt=6,K=6)
#> NK: 1 
#> Number of groups : 6 
#> 1 
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ****________________________________________________****
#> 
#> 2 
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ****________________________________________________****
#> 
#> 3 
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ****________________________________________________****
#> 
#> 4 
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> Warning : 1 2 3 4 5 6 7 < 10^{-12}
#> Warning only 5 components could thus be extracted
#> ****________________________________________________****
#> 
#> 5 
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ****________________________________________________****
#> 
#> 6 
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ****________________________________________________****
#> 
res.cv.modpls<-cvtable(summary(cv.modpls))
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> 
#> CV Q2 criterion:
#> 0 1 
#> 0 1 
#> 
#> CV Press criterion:
#> 1 2 3 
#> 0 0 1 
```
