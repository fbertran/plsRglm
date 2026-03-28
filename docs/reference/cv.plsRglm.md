# Partial least squares regression glm models with k-fold cross validation

This function implements k-fold cross-validation on complete or
incomplete datasets for partial least squares regression generalized
linear models

## Usage

``` r
cv.plsRglm(object, ...)
# Default S3 method
cv.plsRglmmodel(object,dataX,nt=2,limQ2set=.0975,
modele="pls", family=NULL, K=5, NK=1, grouplist=NULL, random=TRUE, 
scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, 
keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights, method,
fit_backend="stats",verbose=TRUE,...)
# S3 method for class 'formula'
cv.plsRglmmodel(object,data=NULL,nt=2,limQ2set=.0975,
modele="pls", family=NULL, K=5, NK=1, grouplist=NULL, random=TRUE, 
scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, 
keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12),weights,subset,
start=NULL,etastart,mustart,offset,method,control= list(),contrasts=NULL,
fit_backend="stats",verbose=TRUE,...)
PLS_glm_kfoldcv(dataY, dataX, nt = 2, limQ2set = 0.0975, modele = "pls", 
family = NULL, K = 5, NK = 1, grouplist = NULL, random = TRUE, 
scaleX = TRUE, scaleY = NULL, keepcoeffs = FALSE, keepfolds = FALSE, 
keepdataY = TRUE, keepMclassed=FALSE, tol_Xi = 10^(-12), weights, method,
fit_backend="stats",verbose=TRUE)
PLS_glm_kfoldcv_formula(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls",
family=NULL, K=5, NK=1, grouplist=NULL, random=TRUE, 
scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12),weights,subset,start=NULL,etastart,
mustart,offset,method,control= list(),contrasts=NULL, fit_backend="stats",
verbose=TRUE)
```

## Arguments

- object:

  response (training) dataset or an object of class
  "[`formula`](https://rdrr.io/r/stats/formula.html)" (or one that can
  be coerced to that class): a symbolic description of the model to be
  fitted. The details of model specification are given under 'Details'.

- dataY:

  response (training) dataset

- dataX:

  predictor(s) (training) dataset

- formula:

  an object of class "[`formula`](https://rdrr.io/r/stats/formula.html)"
  (or one that can be coerced to that class): a symbolic description of
  the model to be fitted. The details of model specification are given
  under 'Details'.

- data:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in `data`,
  the variables are taken from `environment(formula)`, typically the
  environment from which `plsRglm` is called.

- nt:

  number of components to be extracted

- limQ2set:

  limit value for the Q2

- modele:

  name of the PLS glm model to be fitted (`"pls"`, `"pls-glm-Gamma"`,
  `"pls-glm-gaussian"`, `"pls-glm-inverse.gaussian"`,
  `"pls-glm-logistic"`, `"pls-glm-poisson"`, `"pls-glm-polr"`). Use
  `"modele=pls-glm-family"` to enable the `family` option.

- family:

  a description of the error distribution and link function to be used
  in the model. This can be a character string naming a family function,
  a family function or the result of a call to a family function. (See
  [`family`](https://rdrr.io/r/stats/family.html) for details of family
  functions.) To use the family option, please set
  `modele="pls-glm-family"`. User defined families can also be defined.
  See details.

- K:

  number of groups. Defaults to 5.

- NK:

  number of times the group division is made

- grouplist:

  to specify the members of the `K` groups

- random:

  should the `K` groups be made randomly. Defaults to `TRUE`

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since non always possible for glm
  responses.

- keepcoeffs:

  shall the coefficients for each model be returned

- keepfolds:

  shall the groups' composition be returned

- keepdataY:

  shall the observed value of the response for each one of the predicted
  value be returned

- keepMclassed:

  shall the number of miss classed be returned (unavailable)

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- fit_backend:

  backend used for repeated non-ordinal score-space GLM fits during
  cross-validation. Use `"stats"` for the compatibility path or
  `"fastglm"` to opt into the accelerated complete-data backend.
  Unsupported cases fall back to `"stats"` with a warning.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- start:

  starting values for the parameters in the linear predictor.

- etastart:

  starting values for the linear predictor.

- mustart:

  starting values for the vector of means.

- offset:

  this can be used to specify an *a priori* known component to be
  included in the linear predictor during fitting. This should be `NULL`
  or a numeric vector of length equal to the number of cases. One or
  more [`offset`](https://rdrr.io/r/stats/offset.html) terms can be
  included in the formula instead or as well, and if more than one is
  specified their sum is used. See
  [`model.offset`](https://rdrr.io/r/stats/model.extract.html).

- method:

  For non-ordinal GLM modes this argument is kept for backward
  compatibility; use `fit_backend` to choose the score-space fitting
  backend. For `pls-glm-polr`, use `logistic`, `probit`, complementary
  log-log or `cauchit`.

- control:

  a list of parameters for controlling the fitting process. For
  `glm.fit` this is passed to
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html).

- contrasts:

  an optional list. See the `contrasts.arg` of `model.matrix.default`.

- verbose:

  should info messages be displayed ?

- ...:

  arguments to pass to `cv.plsRglmmodel.default` or to
  `cv.plsRglmmodel.formula`

## Details

Predicts 1 group with the `K-1` other groups. Leave one out cross
validation is thus obtained for `K==nrow(dataX)`.

There are seven different predefined models with predefined link
functions available :

- `"pls"`:

  ordinary pls models

- `"pls-glm-Gamma"`:

  glm gaussian with inverse link pls models

- `"pls-glm-gaussian"`:

  glm gaussian with identity link pls models

- `"pls-glm-inverse-gamma"`:

  glm binomial with square inverse link pls models

- `"pls-glm-logistic"`:

  glm binomial with logit link pls models

- `"pls-glm-poisson"`:

  glm poisson with log link pls models

- `"pls-glm-polr"`:

  glm polr with logit link pls models

Using the `"family="` option and setting `"modele=pls-glm-family"`
allows changing the family and link function the same way as for the
[`glm`](https://rdrr.io/r/stats/glm.html) function. As a consequence
user-specified families can also be used.

- The `gaussian` family:

  accepts the links (as names) `identity`, `log` and `inverse`.

- The `binomial` family:

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- The `Gamma` family:

  accepts the links `inverse`, `identity` and `log`.

- The `poisson` family:

  accepts the links `log`, `identity`, and `sqrt`.

- The `inverse.gaussian` family:

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- The `quasi` family:

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- The function `power`:

  can be used to create a power link function.

- ...:

  arguments to pass to `cv.plsRglmmodel.default` or to
  `cv.plsRglmmodel.formula`

A typical predictor has the form response ~ terms where response is the
(numeric) response vector and terms is a series of terms which specifies
a linear predictor for response. A terms specification of the form
first + second indicates all the terms in first together with all the
terms in second with any duplicates removed.

A specification of the form first:second indicates the the set of terms
obtained by taking the interactions of all terms in first with all terms
in second. The specification first\*second indicates the cross of first
and second. This is the same as first + second + first:second.

The terms in the formula will be re-ordered so that main effects come
first, followed by the interactions, all second-order, all third-order
and so on: to avoid this pass a terms object as the formula.

Non-NULL weights can be used to indicate that different observations
have different dispersions (with the values in weights being inversely
proportional to the dispersions); or equivalently, when the elements of
weights are positive integers w_i, that each response y_i is the mean of
w_i unit-weight observations.

## Value

An object of class `"cv.plsRglmmodel"`.  

- results_kfolds:

  list of `NK`. Each element of the list sums up the results for a group
  division:

  list

  :   of `K` matrices of size about `nrow(dataX)/K * nt` with the
      predicted values for a growing number of components

  ...

  :   ...

  list

  :   of `K` matrices of size about `nrow(dataX)/K * nt` with the
      predicted values for a growing number of components

- folds:

  list of `NK`. Each element of the list sums up the informations for a
  group division:

  list

  :   of `K` vectors of length about `nrow(dataX)` with the numbers of
      the rows of `dataX` that were used as a training set

  ...

  :   ...

  list

  :   of `K` vectors of length about `nrow(dataX)` with the numbers of
      the rows of `dataX` that were used as a training set

- dataY_kfolds:

  list of `NK`. Each element of the list sums up the results for a group
  division:

  list

  :   of `K` matrices of size about `nrow(dataX)/K * 1` with the
      observed values of the response

  ...

  :   ...

  list

  :   of `K` matrices of size about `nrow(dataX)/K * 1` with the
      observed values of the response

- fit_backend:

  backend used for repeated non-ordinal score-space GLM fits during
  cross-validation

- call:

  the call of the function

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frederic Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.

## Author

Frederic Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Note

Work for complete and incomplete datasets.

## See also

Summary method `summary.cv.plsRglmmodel`.
[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.md),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.md),
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.md)
and [`summary`](https://rdrr.io/r/base/summary.html) to extract and
transform results from k-fold cross validation.

## Examples

``` r
data(Cornell)
bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10)
#> 
#> Model: pls 
#> 
#> NK: 1 
#> Number of groups : 5 
#> 1 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
#> 2 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
#> 3 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> Warning :  < 10^{-12}
#> Warning only 5 components could thus be extracted
#> ****________________________________________________****
#> 
#> 4 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
#> 5 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
(sum1<-summary(bbb))
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X or Y____
#> Loading required namespace: plsdof
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y       Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA         NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.8873612  0.0975  0.8873612 52.69207  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8551631  0.0975 -0.2858523 45.95956  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.6415328  0.0975 -1.4749712 27.38953   4.418081 0.9905556
#> Nb_Comp_4 33.76477 -0.7878496  0.0975 -3.9874852 22.03512   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -8.5815352  0.0975 -4.3592512 23.09427   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975         NA       NA   3.496074 0.9925265
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000002    0.7633343  0.9711321  1.1359501 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"
cvtable(sum1)
#> 
#> CV Q2 criterion:
#> 0 1 
#> 0 1 
#> 
#> CV Press criterion:
#> 1 2 3 4 
#> 0 0 0 1 

bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=3,
modele="pls-glm-family",family=gaussian(),K=12,verbose=FALSE)
(sum2<-summary(bbb2))
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
cvtable(sum2)
#> 
#> CV Q2Chi2 criterion:
#> 0 
#> 0 
#> 
#> CV PreChi2 criterion:
#> 1 
#> 0 

# \donttest{
#random=TRUE is the default to randomly create folds for repeated CV
bbb3 <- cv.plsRglm(Y~.,data=Cornell,nt=3,
modele="pls-glm-family",family=gaussian(),K=6,NK=10, verbose=FALSE)
(sum3<-summary(bbb3))
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1,  2,  3,  4,  5,  6,  7,  8,  9,  10
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[2]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[3]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[4]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[5]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[6]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[7]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[8]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[9]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> [[10]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 53.15173 54.60645           NA 0.0975        NA                NA
#> Nb_Comp_2 31.46903 33.40866           NA 0.0975        NA                NA
#> Nb_Comp_3 31.54404 33.96857           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
plot(cvtable(sum3))
#> 
#> CV Q2Chi2 criterion:
#> 0 
#> 0 
#> 
#> CV PreChi2 criterion:
#> 1 
#> 0 


data(aze_compl)
bbb <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,modele="pls",keepcoeffs=TRUE, verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>            [,1]        [,2]      [,3]       [,4]       [,5]         [,6]
#>  [1,] 0.2967168 -0.13006495 0.4244503 -0.1703251 0.26196401  0.043708007
#>  [2,] 0.2896904 -0.24144221 0.4654035 -0.2464788 0.26538330  0.171971504
#>  [3,] 0.1915698 -0.21842302 0.3791659 -0.1266468 0.44343504  0.139274504
#>  [4,] 0.3857919 -0.07736703 0.5467866 -0.1209728 0.09343678  0.065989000
#>  [5,] 0.4240014 -0.06216565 0.4129669 -0.2048635 0.21303281  0.077136550
#>  [6,] 0.3358596 -0.17211677 0.4105767 -0.2673842 0.21577097 -0.007922469
#>  [7,] 0.2822752 -0.10516666 0.4808965 -0.1494745 0.23466562  0.110696072
#>  [8,] 0.2527636 -0.10510133 0.4854560 -0.1000471 0.36979962  0.180110808
#>  [9,] 0.3993556 -0.16229696 0.4829590 -0.1711831 0.21189078  0.111610023
#> [10,] 0.2895455 -0.02869057 0.4658724 -0.3370568 0.37347644  0.107648292
#>              [,7]        [,8]       [,9]        [,10]       [,11]       [,12]
#>  [1,] -0.02220514 -0.08364986 -0.1979101  0.093472336 -0.06197899 0.056564347
#>  [2,]  0.05542231 -0.05243850 -0.2936932  0.076489521  0.01192569 0.154965271
#>  [3,] -0.08183618  0.02212690 -0.1642130 -0.004873956 -0.07652468 0.096745405
#>  [4,] -0.05663579  0.17570620 -0.1659681  0.039799153 -0.10128529 0.076611388
#>  [5,] -0.03846534  0.03037289 -0.1640829  0.013503127 -0.10236522 0.005930437
#>  [6,] -0.04313130  0.03739854 -0.1665734  0.165789571 -0.13665543 0.127189948
#>  [7,] -0.07223025  0.03784306 -0.1825016  0.016000105 -0.20028630 0.026568830
#>  [8,] -0.07565672 -0.14560626 -0.3282638  0.060146981 -0.16980593 0.006667436
#>  [9,] -0.10825685 -0.04553849 -0.1865241  0.018216191 -0.09792577 0.053831140
#> [10,] -0.06089697  0.07826849 -0.2590398 -0.002300607 -0.03250620 0.030586975
#>             [,13]      [,14]      [,15]       [,16]         [,17]      [,18]
#>  [1,] -0.07510736 0.23078667 0.10563173  0.18772625 -0.0300309287 0.21263873
#>  [2,] -0.23749636 0.04558970 0.08454484  0.14967967 -0.0004375664 0.39713371
#>  [3,] -0.12474268 0.04938246 0.15465556  0.03582883 -0.0059411494 0.26670994
#>  [4,] -0.15997105 0.08519115 0.18992075 -0.02233185  0.0599327774 0.25055260
#>  [5,] -0.09533200 0.07724877 0.16879089  0.02136612 -0.0049125484 0.19372702
#>  [6,] -0.13975090 0.11795042 0.08638693 -0.08176583  0.0465160346 0.20726721
#>  [7,] -0.13257979 0.06708428 0.07721613  0.14955894  0.0505285909 0.23974610
#>  [8,] -0.20749100 0.11427235 0.04686166  0.10262135  0.0530808304 0.09647482
#>  [9,] -0.10645950 0.04433262 0.17645282  0.09795228 -0.0541015653 0.23760656
#> [10,] -0.10994264 0.18823880 0.03820209 -0.05529992  0.0192444725 0.27622868
#>             [,19]        [,20]       [,21]       [,22]     [,23]       [,24]
#>  [1,] -0.07091179  0.109706049 -0.16681496  0.15284757 0.0754412 -0.06053938
#>  [2,] -0.03976640 -0.047429151 -0.19564985  0.11886939 0.1403371 -0.10750491
#>  [3,]  0.09536022  0.005360878 -0.16706574  0.04188904 0.2113393 -0.10795826
#>  [4,]  0.05422685  0.127388564 -0.16523735 -0.12147590 0.2623665 -0.14394757
#>  [5,] -0.02504900  0.079132188 -0.09007083  0.01233525 0.1854891 -0.07745486
#>  [6,]  0.19132045  0.130006132 -0.06078700  0.13899187 0.2478101 -0.17797556
#>  [7,] -0.02385385 -0.003605265 -0.13019927  0.11793327 0.1268651 -0.12357718
#>  [8,] -0.04146681  0.147465934 -0.18001147  0.10664554 0.2573214 -0.20949431
#>  [9,] -0.04400643  0.125181681 -0.06714280  0.04722346 0.1456566 -0.17888290
#> [10,]  0.03539030  0.048193158 -0.03333034  0.04484396 0.2358709 -0.25478533
#>             [,25]      [,26]     [,27]     [,28]       [,29]        [,30]
#>  [1,] -0.09798007 -0.2340953 0.1002317 0.1563176 -0.12281970 -0.002163329
#>  [2,] -0.08735315 -0.3198482 0.3133822 0.2198037 -0.12634550 -0.050496036
#>  [3,] -0.13482254 -0.2259121 0.1566395 0.2201068 -0.13156798 -0.026562696
#>  [4,] -0.23150811 -0.4056949 0.2175567 0.1802417 -0.07018451  0.019908433
#>  [5,] -0.18355385 -0.3447394 0.2071839 0.1330598 -0.11721313 -0.011874481
#>  [6,] -0.19855659 -0.3326028 0.2079497 0.1377566 -0.17212351 -0.061590527
#>  [7,] -0.23926395 -0.2702062 0.1546947 0.1956377 -0.10817234  0.063295607
#>  [8,] -0.19640686 -0.2375795 0.1940530 0.2020997  0.02779777  0.064559095
#>  [9,] -0.24837769 -0.1914205 0.1425370 0.2044805 -0.06714017  0.116370995
#> [10,] -0.20671878 -0.3224187 0.1734143 0.2649366 -0.07129189  0.062303086
#>             [,31]       [,32]      [,33]        [,34]
#>  [1,]  0.04193639  0.03039599 -0.3734408 -0.099984641
#>  [2,]  0.02320407 -0.11157988 -0.2941528  0.026787856
#>  [3,]  0.20541643 -0.06590305 -0.4589991 -0.043495605
#>  [4,] -0.02965628  0.09875478 -0.4339590 -0.097001286
#>  [5,]  0.11051232 -0.05876086 -0.3692653  0.017894402
#>  [6,]  0.22400399 -0.10653138 -0.3810163  0.022816898
#>  [7,]  0.18937429  0.01734272 -0.3555912  0.027072707
#>  [8,]  0.15928404 -0.07357715 -0.3343743  0.029830172
#>  [9,]  0.13015183 -0.11073135 -0.3820088 -0.007858931
#> [10,]  0.17073592 -0.13172571 -0.5408630  0.034492209
bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,modele="pls-glm-family",
family=binomial(probit),keepcoeffs=TRUE, verbose=FALSE)
bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,
modele="pls-glm-logistic",keepcoeffs=TRUE, verbose=FALSE)
summary(bbb,MClassed=TRUE)
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC MissClassed CV_MissClassed       Q2cum_Y LimQ2_Y       Q2_Y
#> Nb_Comp_0  154.6179          49             NA            NA      NA         NA
#> Nb_Comp_1  126.4083          27             47    -0.1316142  0.0975 -0.1316142
#> Nb_Comp_2  119.3375          25             50    -0.8611003  0.0975 -0.6446420
#> Nb_Comp_3  114.2313          27             49    -2.6590200  0.0975 -0.9660520
#> Nb_Comp_4  112.3463          23             49    -7.6014164  0.0975 -1.3507432
#> Nb_Comp_5  113.2362          22             51   -20.7208349  0.0975 -1.5252626
#> Nb_Comp_6  114.7620          21             50   -55.0297514  0.0975 -1.5795395
#> Nb_Comp_7  116.5264          20             51  -144.1617949  0.0975 -1.5907985
#> Nb_Comp_8  118.4601          20             51  -376.2570500  0.0975 -1.5988729
#> Nb_Comp_9  120.4452          19             50  -982.6365094  0.0975 -1.6073376
#> Nb_Comp_10 122.4395          19             50 -2568.2986579  0.0975 -1.6120408
#>             PRESS_Y    RSS_Y      R2_Y  AIC.std  DoF.dof sigmahat.dof   AIC.dof
#> Nb_Comp_0        NA 25.91346        NA 298.1344  1.00000    0.5015845 0.2540061
#> Nb_Comp_1  29.32404 19.38086 0.2520929 269.9248 22.55372    0.4848429 0.2883114
#> Nb_Comp_2  31.87458 17.76209 0.3145613 262.8540 27.31542    0.4781670 0.2908950
#> Nb_Comp_3  34.92119 16.58896 0.3598323 257.7478 30.52370    0.4719550 0.2902572
#> Nb_Comp_4  38.99639 15.98071 0.3833049 255.8628 34.00000    0.4744263 0.3008285
#> Nb_Comp_5  40.35548 15.81104 0.3898523 256.7527 34.00000    0.4719012 0.2976347
#> Nb_Comp_6  40.78520 15.73910 0.3926285 258.2785 34.00000    0.4708264 0.2962804
#> Nb_Comp_7  40.77683 15.70350 0.3940024 260.0429 33.71066    0.4693382 0.2937976
#> Nb_Comp_8  40.81139 15.69348 0.3943888 261.9766 34.00000    0.4701436 0.2954217
#> Nb_Comp_9  40.91821 15.69123 0.3944758 263.9617 33.87284    0.4696894 0.2945815
#> Nb_Comp_10 40.98613 15.69037 0.3945088 265.9560 34.00000    0.4700970 0.2953632
#>              BIC.dof  GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive
#> Nb_Comp_0  0.2604032 -67.17645         1      0.5015845 0.2540061 0.2604032
#> Nb_Comp_1  0.4231184 -53.56607         2      0.4358996 0.1936625 0.2033251
#> Nb_Comp_2  0.4496983 -52.42272         3      0.4193593 0.1809352 0.1943501
#> Nb_Comp_3  0.4631316 -51.93343         4      0.4072955 0.1722700 0.1891422
#> Nb_Comp_4  0.4954133 -50.37079         5      0.4017727 0.1691819 0.1897041
#> Nb_Comp_5  0.4901536 -50.65724         6      0.4016679 0.1706451 0.1952588
#> Nb_Comp_6  0.4879234 -50.78005         7      0.4028135 0.1731800 0.2020601
#> Nb_Comp_7  0.4826103 -51.05525         8      0.4044479 0.1761610 0.2094352
#> Nb_Comp_8  0.4865092 -50.85833         9      0.4064413 0.1794902 0.2172936
#> Nb_Comp_9  0.4845867 -50.95616        10      0.4085682 0.1829787 0.2254232
#> Nb_Comp_10 0.4864128 -50.86368        11      0.4107477 0.1865584 0.2337468
#>            GMDL.naive
#> Nb_Comp_0   -67.17645
#> Nb_Comp_1   -79.67755
#> Nb_Comp_2   -81.93501
#> Nb_Comp_3   -83.31503
#> Nb_Comp_4   -83.23369
#> Nb_Comp_5   -81.93513
#> Nb_Comp_6   -80.42345
#> Nb_Comp_7   -78.87607
#> Nb_Comp_8   -77.31942
#> Nb_Comp_9   -75.80069
#> Nb_Comp_10  -74.33325
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"
summary(bbb2,MClassed=TRUE)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC MissClassed CV_MissClassed Q2Chisqcum_Y  limQ2
#> Nb_Comp_0  145.8283 148.4727          49             NA           NA     NA
#> Nb_Comp_1  118.1398 123.4285          28             NA           NA 0.0975
#> Nb_Comp_2  109.9553 117.8885          26             NA           NA 0.0975
#> Nb_Comp_3  105.1591 115.7366          22             NA           NA 0.0975
#> Nb_Comp_4  103.8382 117.0601          21             NA           NA 0.0975
#> Nb_Comp_5  104.7338 120.6001          21             NA           NA 0.0975
#> Nb_Comp_6  105.6770 124.1878          21             NA           NA 0.0975
#> Nb_Comp_7  107.2828 128.4380          20             NA           NA 0.0975
#> Nb_Comp_8  109.0172 132.8167          22             NA           NA 0.0975
#> Nb_Comp_9  110.9354 137.3793          21             NA           NA 0.0975
#> Nb_Comp_10 112.9021 141.9904          20             NA           NA 0.0975
#>            Q2Chisq_Y PREChi2_Pearson_Y Chi2_Pearson_Y    RSS_Y      R2_Y
#> Nb_Comp_0         NA                NA      104.00000 25.91346        NA
#> Nb_Comp_1         NA                NA      100.53823 19.32272 0.2543365
#> Nb_Comp_2         NA                NA       99.17955 17.33735 0.3309519
#> Nb_Comp_3         NA                NA      123.37836 15.58198 0.3986915
#> Nb_Comp_4         NA                NA      114.77551 15.14046 0.4157299
#> Nb_Comp_5         NA                NA      105.35382 15.08411 0.4179043
#> Nb_Comp_6         NA                NA       98.87767 14.93200 0.4237744
#> Nb_Comp_7         NA                NA       97.04072 14.87506 0.4259715
#> Nb_Comp_8         NA                NA       98.90110 14.84925 0.4269676
#> Nb_Comp_9         NA                NA      100.35563 14.84317 0.4272022
#> Nb_Comp_10        NA                NA      102.85214 14.79133 0.4292027
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
kfolds2coeff(bbb2)
#>            [,1]       [,2]     [,3]       [,4]     [,5]        [,6]
#>  [1,] -2.880814 -1.6610263 3.231294 -1.7859007 2.157633  0.09432333
#>  [2,] -1.837237 -0.5422735 5.034633 -1.6529388 2.197167 -0.50874173
#>  [3,] -3.131430 -1.0823842 4.640383 -1.9632623 3.088607  0.68644255
#>  [4,] -3.354528 -1.3353782 3.834992 -2.0786740 2.428283  0.33873018
#>  [5,] -3.245508 -1.7768940 4.552638 -0.9719501 3.327824  0.76224721
#>  [6,] -2.308158 -1.8922548 3.744317 -1.8196755 2.538591  1.46435033
#>  [7,] -1.650276 -1.2782098 2.542370 -2.9520691 3.130694  1.36125422
#>  [8,] -1.728052 -0.6985864 3.971841 -1.4811625 0.654456  0.43064899
#>  [9,] -3.277268 -1.6679437 4.135433 -2.5336136 3.019127  0.82660980
#> [10,] -2.377851 -0.9410133 3.222897 -1.8056744 2.102928  0.48978524
#>               [,7]        [,8]       [,9]      [,10]      [,11]     [,12]
#>  [1,] -1.997484018  1.22591974 -0.7192582 0.42028283 -0.2706868 0.7392970
#>  [2,]  0.155967036 -0.11715847 -1.9112048 0.45891383 -0.8413667 0.1754463
#>  [3,]  0.340548203 -0.45334549 -2.5276452 0.08605549 -0.4430792 1.0583928
#>  [4,]  0.276769867 -1.37546144 -0.6611894 0.79977621 -0.7732384 0.4681988
#>  [5,] -0.088675517 -1.46521504 -2.5284904 0.06456009  0.1442771 1.4479132
#>  [6,] -0.111123814  0.03287493 -1.4500170 0.43603196 -1.2748238 1.2615131
#>  [7,] -0.296940964  0.78038251 -0.8419953 0.33314661 -1.6486065 0.6118850
#>  [8,]  0.008807352 -0.21504158 -0.2416204 0.50397000 -1.1884159 0.8579691
#>  [9,]  1.360703818 -1.14329026 -2.3391000 0.06476553 -1.0371017 0.6106706
#> [10,] -0.107797270 -0.25918784 -1.9158962 0.99421270 -0.7206721 0.4275451
#>            [,13]       [,14]     [,15]     [,16]     [,17]    [,18]       [,19]
#>  [1,] -0.4278569  1.84977208 0.7213886 0.4711776 0.8903936 1.478132  0.33921594
#>  [2,] -0.9246932  0.71166529 1.3680333 0.2241439 0.8338250 1.565577  0.13778568
#>  [3,] -0.5260922  1.37706211 0.2560763 0.7542326 0.6278519 1.423074 -0.07960995
#>  [4,] -0.4037473  1.10192476 0.4049219 0.9119733 0.9287992 2.226056 -1.57304555
#>  [5,] -2.1712519 -0.18420157 1.6775261 0.1399884 0.1444457 2.161048 -0.46876619
#>  [6,] -1.6532341  0.56014290 0.8652359 0.1954831 0.6710559 2.200706  0.13901295
#>  [7,] -0.7677609 -0.02597986 1.3072229 1.2606164 0.1696415 3.118533 -0.45016603
#>  [8,] -0.9129750  0.11012833 0.7019751 0.3357959 0.8365057 1.511821 -0.14507655
#>  [9,] -1.2165418  0.14458934 0.8904357 1.0063786 1.1406017 1.618479  1.02449582
#> [10,] -1.4244882  0.82272847 0.4990603 0.3643162 0.6079886 1.856440  0.42641384
#>           [,20]       [,21]       [,22]     [,23]      [,24]      [,25]
#>  [1,] 0.2575058  0.01726207  0.50918284 1.4893868 -1.8358027 -1.8616799
#>  [2,] 0.3253417 -0.55049975 -0.25476628 2.0204762 -1.5714178 -2.4341169
#>  [3,] 0.7697634 -0.47712966  0.68723657 1.2078782 -1.6740166 -1.9113482
#>  [4,] 0.7842468 -1.01466470  0.77432801 0.6910830 -0.8849110 -1.5607995
#>  [5,] 0.7868114 -1.25996659 -0.54563711 2.1526112 -0.7454372  0.1053691
#>  [6,] 0.4660470 -1.54183793  0.05692282 1.7389029 -0.6756828 -1.6862276
#>  [7,] 0.3718957 -1.23243175  1.07443205 0.8696556 -2.5510579 -2.0783182
#>  [8,] 1.2411450 -1.23903040  0.63379607 1.2300652 -1.2413804 -1.9471408
#>  [9,] 0.2915076 -2.54630562  0.48117693 1.7271066 -0.6554333 -0.9174380
#> [10,] 0.5668408 -0.73797468  0.79390221 1.3048951 -0.7869658 -1.4018839
#>           [,26]    [,27]    [,28]      [,29]       [,30]     [,31]       [,32]
#>  [1,] -2.122028 2.140007 1.968188 -2.0058317 -0.29732695 1.6882048 -0.26707951
#>  [2,] -2.098528 1.306625 1.209686 -0.7983181  0.62436327 0.9427557  0.16043806
#>  [3,] -2.397412 1.536631 2.218406 -0.9821636  0.79287028 1.6066021 -1.03753852
#>  [4,] -1.149318 1.370820 1.302808 -0.4334963  0.02000343 1.6506259  0.32427436
#>  [5,] -2.655590 1.670839 1.697827 -0.6602916  0.30060207 0.5024129 -0.64215111
#>  [6,] -2.502821 1.489366 1.737182 -0.4690624  0.29151549 1.7562107  0.17878560
#>  [7,] -2.762351 2.333225 1.527855  0.3855136  1.38226207 1.5746615 -1.29099677
#>  [8,] -2.169963 1.910589 1.780115 -0.1823360 -0.03561015 1.2138244 -0.06049045
#>  [9,] -1.833474 2.366218 1.711345 -0.6389232  1.27316103 1.8570508 -0.07265236
#> [10,] -2.217370 1.053345 1.150473 -0.2554370  0.44501671 1.3638434 -0.11398587
#>           [,33]       [,34]
#>  [1,] -3.404842  0.66893826
#>  [2,] -3.175950 -0.24471695
#>  [3,] -3.536399 -0.02685516
#>  [4,] -3.457776  0.66964879
#>  [5,] -2.170043  0.54111046
#>  [6,] -3.559590 -0.41736380
#>  [7,] -3.804775  0.22613603
#>  [8,] -3.389891 -0.25626557
#>  [9,] -3.807447  0.24450051
#> [10,] -2.120071 -0.36804768

kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0  145.8283 148.4727           NA     NA        NA                NA
#> Nb_Comp_1  118.1398 123.4285           NA 0.0975        NA                NA
#> Nb_Comp_2  109.9553 117.8885           NA 0.0975        NA                NA
#> Nb_Comp_3  105.1591 115.7366           NA 0.0975        NA                NA
#> Nb_Comp_4  103.8382 117.0601           NA 0.0975        NA                NA
#> Nb_Comp_5  104.7338 120.6001           NA 0.0975        NA                NA
#> Nb_Comp_6  105.6770 124.1878           NA 0.0975        NA                NA
#> Nb_Comp_7  107.2828 128.4380           NA 0.0975        NA                NA
#> Nb_Comp_8  109.0172 132.8167           NA 0.0975        NA                NA
#> Nb_Comp_9  110.9354 137.3793           NA 0.0975        NA                NA
#> Nb_Comp_10 112.9021 141.9904           NA 0.0975        NA                NA
#>            Chi2_Pearson_Y    RSS_Y      R2_Y
#> Nb_Comp_0       104.00000 25.91346        NA
#> Nb_Comp_1       100.53823 19.32272 0.2543365
#> Nb_Comp_2        99.17955 17.33735 0.3309519
#> Nb_Comp_3       123.37836 15.58198 0.3986915
#> Nb_Comp_4       114.77551 15.14046 0.4157299
#> Nb_Comp_5       105.35382 15.08411 0.4179043
#> Nb_Comp_6        98.87767 14.93200 0.4237744
#> Nb_Comp_7        97.04072 14.87506 0.4259715
#> Nb_Comp_8        98.90110 14.84925 0.4269676
#> Nb_Comp_9       100.35563 14.84317 0.4272022
#> Nb_Comp_10      102.85214 14.79133 0.4292027
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
rm(list=c("bbb","bbb2"))



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,modele="pls-glm-family",
family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,
modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>           [,1]         [,2]        [,3]      [,4]      [,5]      [,6]
#>  [1,] 16.88984 -0.006649928 -0.09824786 0.3363855 -2.214567 0.4159321
#>  [2,] 11.98800 -0.005092515 -0.07197508 0.1843301 -1.570331 0.3097856
#>  [3,] 10.33640 -0.005929898 -0.04154327 0.0237779 -1.607505 0.3175666
#>  [4,] 11.45685 -0.004411421 -0.05812383 0.1824373 -1.679110 0.3215815
#>  [5,] 13.09302 -0.004481021 -0.09047226 0.2081762 -2.341165 0.4426191
#>  [6,] 10.22433 -0.003429086 -0.05768415 0.1086879 -1.014885 0.2056477
#>  [7,] 12.75433 -0.004951018 -0.07474571 0.2162661 -1.649608 0.3307475
#>  [8,] 13.13244 -0.005813069 -0.06263790 0.2590391 -1.266201 0.2807201
#>  [9,] 10.64302 -0.004168722 -0.07657230 0.1395666 -1.426423 0.2982473
#> [10,] 15.02271 -0.005935281 -0.08119719 0.2504480 -1.930864 0.3405278
#>              [,7]        [,8]       [,9]      [,10]        [,11]
#>  [1,] -4.27141207  0.58217503 0.37162725 -0.9816309 -0.579296944
#>  [2,] -2.24528369  0.16957374 0.29069758 -1.0190103 -0.196499262
#>  [3,] -0.03758572  0.46330982 0.18723152 -1.1430583 -0.083563069
#>  [4,] -1.86830180  0.17097922 0.17589886 -1.2692590 -0.052901073
#>  [5,] -2.36307426  0.35268345 0.37128287 -1.3338238 -0.501250589
#>  [6,] -0.51427080 -0.66411407 0.08654318 -1.6578232 -0.056631993
#>  [7,] -2.75222087  0.44458319 0.12132515 -0.4196716 -0.582238180
#>  [8,] -2.96329854  0.12802673 0.23897351 -1.6034877  0.209432344
#>  [9,] -1.24219208 -0.06630621 0.12359687 -1.0839730  0.009831884
#> [10,] -2.81795223  0.29848599 0.38505722 -1.4699334 -0.454488794
boxplot(kfolds2coeff(bbb)[,1])


kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
summary(bbb)
#> ____************************************************____
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0  76.61170 78.10821           NA     NA        NA                NA
#> Nb_Comp_1  65.70029 68.69331           NA 0.0975        NA                NA
#> Nb_Comp_2  62.49440 66.98392           NA 0.0975        NA                NA
#> Nb_Comp_3  62.47987 68.46590           NA 0.0975        NA                NA
#> Nb_Comp_4  64.21704 71.69958           NA 0.0975        NA                NA
#> Nb_Comp_5  65.81654 74.79559           NA 0.0975        NA                NA
#> Nb_Comp_6  66.48888 76.96443           NA 0.0975        NA                NA
#> Nb_Comp_7  68.40234 80.37440           NA 0.0975        NA                NA
#> Nb_Comp_8  70.39399 83.86256           NA 0.0975        NA                NA
#> Nb_Comp_9  72.37642 87.34149           NA 0.0975        NA                NA
#> Nb_Comp_10 74.37612 90.83770           NA 0.0975        NA                NA
#>            Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0        33.75000 24.545455        NA
#> Nb_Comp_1        23.85891 12.599337 0.4866937
#> Nb_Comp_2        17.29992  9.056074 0.6310488
#> Nb_Comp_3        15.50937  8.232069 0.6646194
#> Nb_Comp_4        15.23934  8.125808 0.6689485
#> Nb_Comp_5        15.26275  7.862134 0.6796909
#> Nb_Comp_6        17.74629  6.203270 0.7472742
#> Nb_Comp_7        18.04460  5.879880 0.7604493
#> Nb_Comp_8        18.17881  5.827065 0.7626011
#> Nb_Comp_9        18.34925  5.837300 0.7621841
#> Nb_Comp_10       18.39332  5.832437 0.7623822
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm(ypine,Xpine,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                 AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0  82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1  63.61896  0.38248575  0.0975  0.38248575 12.844390 11.074659
#> Nb_Comp_2  58.47638  0.34836456  0.0975 -0.05525570 11.686597  8.919303
#> Nb_Comp_3  56.55421  0.23688359  0.0975 -0.17107874 10.445206  7.919786
#> Nb_Comp_4  54.35053  0.06999681  0.0975 -0.21869112  9.651773  6.972542
#> Nb_Comp_5  55.99834 -0.07691053  0.0975 -0.15796434  8.073955  6.898523
#> Nb_Comp_6  57.69592 -0.19968885  0.0975 -0.11400977  7.685022  6.835594
#> Nb_Comp_7  59.37953 -0.27722139  0.0975 -0.06462721  7.277359  6.770369
#> Nb_Comp_8  61.21213 -0.30602578  0.0975 -0.02255238  6.923057  6.736112
#> Nb_Comp_9  63.18426 -0.39920228  0.0975 -0.07134354  7.216690  6.730426
#> Nb_Comp_10 65.15982 -0.43743644  0.0975 -0.02732569  6.914340  6.725443
#>                 R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0         NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1  0.4675684 0.4675684   17.03781     19.76046  0.38248575 0.0975
#> Nb_Comp_2  0.5711905 0.5711905   13.72190     17.97925 -0.05525570 0.0975
#> Nb_Comp_3  0.6192438 0.6192438   12.18420     16.06943 -0.17107874 0.0975
#> Nb_Comp_4  0.6647841 0.6647841   10.72691     14.84877 -0.21869112 0.0975
#> Nb_Comp_5  0.6683426 0.6683426   10.61304     12.42138 -0.15796434 0.0975
#> Nb_Comp_6  0.6713681 0.6713681   10.51622     11.82303 -0.11400977 0.0975
#> Nb_Comp_7  0.6745039 0.6745039   10.41588     11.19586 -0.06462721 0.0975
#> Nb_Comp_8  0.6761508 0.6761508   10.36317     10.65078 -0.02255238 0.0975
#> Nb_Comp_9  0.6764242 0.6764242   10.35443     11.10252 -0.07134354 0.0975
#> Nb_Comp_10 0.6766638 0.6766638   10.34676     10.63737 -0.02732569 0.0975
#>            Q2cum_residY  AIC.std   DoF.dof sigmahat.dof   AIC.dof   BIC.dof
#> Nb_Comp_0            NA 96.63448  1.000000    0.8062287 0.6697018 0.6991787
#> Nb_Comp_1    0.38248575 77.83455  3.176360    0.5994089 0.4047616 0.4565153
#> Nb_Comp_2    0.34836456 72.69198  7.133559    0.5761829 0.4138120 0.5212090
#> Nb_Comp_3    0.23688359 70.76981  8.778329    0.5603634 0.4070516 0.5320535
#> Nb_Comp_4    0.06999681 68.56612  8.427874    0.5221703 0.3505594 0.4547689
#> Nb_Comp_5   -0.07691053 70.21393  9.308247    0.5285695 0.3666578 0.4845912
#> Nb_Comp_6   -0.19968885 71.91152  9.291931    0.5259794 0.3629363 0.4795121
#> Nb_Comp_7   -0.27722139 73.59512  9.756305    0.5284535 0.3702885 0.4938445
#> Nb_Comp_8   -0.30602578 75.42772 10.363948    0.5338475 0.3831339 0.5170783
#> Nb_Comp_9   -0.39920228 77.39986 10.732145    0.5378276 0.3920956 0.5328746
#> Nb_Comp_10  -0.43743644 79.37542 11.000000    0.5407500 0.3987417 0.5446065
#>             GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0  -3.605128         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1  -9.875081         2      0.5977015 0.3788984 0.4112998 -11.451340
#> Nb_Comp_2  -6.985517         3      0.5452615 0.3243383 0.3647862 -12.822703
#> Nb_Comp_3  -6.260610         4      0.5225859 0.3061986 0.3557368 -12.756838
#> Nb_Comp_4  -8.152986         5      0.4990184 0.2867496 0.3432131 -12.811575
#> Nb_Comp_5  -7.111583         6      0.5054709 0.3019556 0.3714754 -11.329638
#> Nb_Comp_6  -7.233043         7      0.5127450 0.3186757 0.4021333  -9.918688
#> Nb_Comp_7  -6.742195         8      0.5203986 0.3364668 0.4347156  -8.592770
#> Nb_Comp_8  -6.038372         9      0.5297842 0.3572181 0.4717708  -7.287834
#> Nb_Comp_9  -5.600238        10      0.5409503 0.3813021 0.5140048  -6.008747
#> Nb_Comp_10 -5.288422        11      0.5529032 0.4076026 0.5600977  -4.799453

data(pineNAX21)
bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
#>           [,1]         [,2]        [,3]         [,4]      [,5]      [,6]
#>  [1,] 12.95747 -0.005687251 -0.05628507  0.211369847 -1.448653 0.3441703
#>  [2,] 13.13803 -0.004172823 -0.07614726  0.230134186 -1.294829 0.2699066
#>  [3,] 18.03596 -0.003017714 -0.07587117  0.254074471 -3.039307 0.4537299
#>  [4,] 14.52213 -0.005066753 -0.07086053  0.215761258 -1.417661 0.2934312
#>  [5,] 19.21123 -0.007247620 -0.10117848  0.355481311 -2.172003 0.4227492
#>  [6,] 15.70833 -0.006167567 -0.06068698  0.249613484 -1.179408 0.2453340
#>  [7,] 11.58985 -0.003979485 -0.06657235  0.135803633 -1.215459 0.2466046
#>  [8,] 12.06044 -0.004590364 -0.07143708  0.194588892 -1.462225 0.2940832
#>  [9,] 10.41200 -0.005600532 -0.02390081 -0.001572928 -1.633021 0.3154482
#> [10,] 14.87391 -0.005071215 -0.07274074  0.239494159 -1.847655 0.3176912
#>             [,7]        [,8]       [,9]      [,10]       [,11]
#>  [1,] -2.2902871  0.31522911 0.11516292 -1.2737166  0.10092851
#>  [2,] -3.2052284  0.04533899 0.11766715 -0.4865288 -0.38456157
#>  [3,] -2.3541610 -0.57438947 0.98171345 -3.9504219 -1.51689439
#>  [4,] -2.3909020 -0.31690760 0.13919108 -1.0371601 -0.26455007
#>  [5,] -4.5423524  0.72099755 0.35837024 -1.0011027 -0.90262771
#>  [6,] -2.3593653  0.13850781 0.11398422 -1.6103403 -0.49817425
#>  [7,] -1.1242876 -0.46062310 0.08910399 -1.2782880  0.14967653
#>  [8,] -1.9968854 -0.04840744 0.16279077 -1.2302761 -0.08283930
#>  [9,]  0.2267452  0.36595690 0.16591772 -1.1682840  0.07720835
#> [10,] -2.6898706 -0.19041266 0.29060144 -1.2908633  0.08834721
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 76.61170 78.10821           NA     NA        NA                NA
#> Nb_Comp_1 65.74449 68.73751           NA 0.0975        NA                NA
#> Nb_Comp_2 62.35674 66.84626           NA 0.0975        NA                NA
#> Nb_Comp_3 62.39804 68.38407           NA 0.0975        NA                NA
#> Nb_Comp_4 64.08113 71.56366           NA 0.0975        NA                NA
#> Nb_Comp_5 65.63784 74.61689           NA 0.0975        NA                NA
#> Nb_Comp_6 67.18468 77.66024           NA 0.0975        NA                NA
#> Nb_Comp_7 68.61004 80.58210           NA 0.0975        NA                NA
#> Nb_Comp_8 70.54487 84.01344           NA 0.0975        NA                NA
#> Nb_Comp_9 72.37296 87.33803           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0       33.75000 24.545455        NA
#> Nb_Comp_1       23.89105 12.654950 0.4844280
#> Nb_Comp_2       17.31172  8.871122 0.6385839
#> Nb_Comp_3       15.51670  8.203709 0.6657748
#> Nb_Comp_4       15.31216  7.959332 0.6757309
#> Nb_Comp_5       15.51159  7.724832 0.6852846
#> Nb_Comp_6       16.30549  6.814620 0.7223673
#> Nb_Comp_7       17.52007  6.284737 0.7439552
#> Nb_Comp_8       17.75766  6.160827 0.7490034
#> Nb_Comp_9       18.30206  5.831059 0.7624383
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"

data(XpineNAX21)
PLS_lm(ypine,XpineNAX21,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> ____There are some NAs in X but not in Y____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE] < 10^{-12}
#> Warning only 9 components could thus be extracted
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#>                AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0 82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1 63.69250  0.35639805  0.0975  0.35639805 13.387018 11.099368
#> Nb_Comp_2 58.35228  0.28395028  0.0975 -0.11256611 12.348781  8.885823
#> Nb_Comp_3 56.36553  0.07664889  0.0975 -0.28950699 11.458331  7.874634
#> Nb_Comp_4 54.02416 -0.70355579  0.0975 -0.84497074 14.528469  6.903925
#> Nb_Comp_5 55.80450 -0.94905654  0.0975 -0.14411078  7.898855  6.858120
#> Nb_Comp_6 57.45753 -1.27568315  0.0975 -0.16758190  8.007417  6.786392
#> Nb_Comp_7 58.73951 -1.63309014  0.0975 -0.15705481  7.852227  6.640327
#> Nb_Comp_8 60.61227 -1.67907859  0.0975 -0.01746558  6.756304  6.614773
#> Nb_Comp_9 62.25948 -2.15165796  0.0975 -0.17639623  7.781594  6.544432
#>                R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0        NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1 0.4663804 0.4663804   17.07583     20.59526  0.35639805 0.0975
#> Nb_Comp_2 0.5728001 0.5728001   13.67040     18.99799 -0.11256611 0.0975
#> Nb_Comp_3 0.6214146 0.6214146   12.11473     17.62807 -0.28950699 0.0975
#> Nb_Comp_4 0.6680830 0.6680830   10.62135     22.35133 -0.84497074 0.0975
#> Nb_Comp_5 0.6702851 0.6702851   10.55088     12.15200 -0.14411078 0.0975
#> Nb_Comp_6 0.6737336 0.6737336   10.44053     12.31901 -0.16758190 0.0975
#> Nb_Comp_7 0.6807558 0.6807558   10.21581     12.08026 -0.15705481 0.0975
#> Nb_Comp_8 0.6819844 0.6819844   10.17650     10.39424 -0.01746558 0.0975
#> Nb_Comp_9 0.6853661 0.6853661   10.06828     11.97160 -0.17639623 0.0975
#>           Q2cum_residY  AIC.std DoF.dof sigmahat.dof AIC.dof BIC.dof GMDL.dof
#> Nb_Comp_0           NA 96.63448      NA           NA      NA      NA       NA
#> Nb_Comp_1   0.35639805 77.90810      NA           NA      NA      NA       NA
#> Nb_Comp_2   0.28395028 72.56787      NA           NA      NA      NA       NA
#> Nb_Comp_3   0.07664889 70.58113      NA           NA      NA      NA       NA
#> Nb_Comp_4  -0.70355579 68.23976      NA           NA      NA      NA       NA
#> Nb_Comp_5  -0.94905654 70.02009      NA           NA      NA      NA       NA
#> Nb_Comp_6  -1.27568315 71.67313      NA           NA      NA      NA       NA
#> Nb_Comp_7  -1.63309014 72.95511      NA           NA      NA      NA       NA
#> Nb_Comp_8  -1.67907859 74.82787      NA           NA      NA      NA       NA
#> Nb_Comp_9  -2.15165796 76.47507      NA           NA      NA      NA       NA
#>           DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1         2      0.5983679 0.3797438 0.4122175 -11.413749
#> Nb_Comp_2         3      0.5442372 0.3231208 0.3634169 -12.847656
#> Nb_Comp_3         4      0.5210941 0.3044529 0.3537087 -12.776843
#> Nb_Comp_4         5      0.4965569 0.2839276 0.3398355 -12.891035
#> Nb_Comp_5         6      0.5039885 0.3001871 0.3692997 -11.349498
#> Nb_Comp_6         7      0.5108963 0.3163819 0.3992388  -9.922119
#> Nb_Comp_7         8      0.5153766 0.3300041 0.4263658  -8.696873
#> Nb_Comp_8         9      0.5249910 0.3507834 0.4632727  -7.337679
#> Nb_Comp_9        10      0.5334234 0.3707649 0.4998004  -6.033403
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))
#> Warning: object 'XpineNAX21' not found



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",
family=Gamma,K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#> Warning: NaNs produced
#> Warning: NaNs produced
bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-Gamma",
K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>             [,1]        [,2]         [,3]       [,4]     [,5]       [,6]
#>  [1,] -10.836184 0.005831805  0.027778503 -0.1009921 1.471098 -0.3392421
#>  [2,] -17.179374 0.007606037  0.042714578 -0.3000448 1.981556 -0.2260388
#>  [3,] -11.265041 0.005206208  0.032087344 -0.2151383 1.316489 -0.3311818
#>  [4,] -12.582702 0.006398389  0.055068952 -0.1745120 2.039713 -0.4257204
#>  [5,] -11.658100 0.005729181  0.038509538 -0.1302073 1.656389 -0.3339884
#>  [6,] -12.070018 0.005420353  0.033023922 -0.2201332 1.467317 -0.2912866
#>  [7,] -11.975554 0.005816951  0.044743545 -0.1772773 1.731737 -0.3438069
#>  [8,]  -8.971586 0.007748655 -0.001844151  0.1855935 1.862339 -0.3806288
#>  [9,] -10.485088 0.005142101  0.040588404 -0.1293531 1.726443 -0.3553153
#> [10,] -11.232821 0.004305256  0.102967565 -0.1114788 2.717317 -0.5413746
#>            [,7]        [,8]         [,9]     [,10]     [,11]
#>  [1,]  1.308277 -0.13032589  0.006163703 0.5243983 0.6496267
#>  [2,]  3.809216 -0.01678435 -0.736814325 2.3413291 0.3561798
#>  [3,]  2.324821  0.07527371  0.003423844 1.2536491 0.2678706
#>  [4,]  2.247434 -0.37109752 -0.171210106 0.6803898 0.4791410
#>  [5,]  1.875297 -0.01975007 -0.249375619 0.9476588 0.6104855
#>  [6,]  2.819522 -0.32812907 -0.193388977 1.0284245 0.8039217
#>  [7,]  2.417400 -0.43933893 -0.257596200 0.9115575 0.7467749
#>  [8,] -1.791605 -0.41413601 -0.198103748 0.4482172 0.6519430
#>  [9,]  1.962835 -0.50552631 -0.163326356 0.4916565 0.8435801
#> [10,]  1.084328 -0.43093998 -0.220912618 0.9994902 0.5148541
boxplot(kfolds2coeff(bbb)[,1])


kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
summary(bbb)
#> ____************************************************____
#> 
#> Family: Gamma 
#> Link function: inverse 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0  56.60919 59.60220           NA     NA        NA                NA
#> Nb_Comp_1  39.01090 43.50042           NA 0.0975        NA                NA
#> Nb_Comp_2  37.30801 43.29404           NA 0.0975        NA                NA
#> Nb_Comp_3  36.87524 44.35777           NA 0.0975        NA                NA
#> Nb_Comp_4  36.55795 45.53700           NA 0.0975        NA                NA
#> Nb_Comp_5  37.13611 47.61167           NA 0.0975        NA                NA
#> Nb_Comp_6  38.27656 50.24862           NA 0.0975        NA                NA
#> Nb_Comp_7  39.39377 52.86234           NA 0.0975        NA                NA
#> Nb_Comp_8  40.96122 55.92630           NA 0.0975        NA                NA
#> Nb_Comp_9  42.90816 59.36974           NA 0.0975        NA                NA
#> Nb_Comp_10 44.90815 62.86625           NA 0.0975        NA                NA
#>            Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0        31.60805 20.800152        NA
#> Nb_Comp_1        17.31431 11.804594 0.4324756
#> Nb_Comp_2        17.01037  6.357437 0.6943562
#> Nb_Comp_3        15.83422  5.699662 0.7259798
#> Nb_Comp_4        13.52676  7.679741 0.6307844
#> Nb_Comp_5        13.60962  6.099077 0.7067773
#> Nb_Comp_6        13.91155  5.205052 0.7497590
#> Nb_Comp_7        14.94390  4.650377 0.7764258
#> Nb_Comp_8        15.25537  4.321314 0.7922461
#> Nb_Comp_9        15.15577  4.307757 0.7928978
#> Nb_Comp_10       15.15490  4.307391 0.7929154
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm(ypine,Xpine,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                 AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0  82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1  63.61896  0.38248575  0.0975  0.38248575 12.844390 11.074659
#> Nb_Comp_2  58.47638  0.34836456  0.0975 -0.05525570 11.686597  8.919303
#> Nb_Comp_3  56.55421  0.23688359  0.0975 -0.17107874 10.445206  7.919786
#> Nb_Comp_4  54.35053  0.06999681  0.0975 -0.21869112  9.651773  6.972542
#> Nb_Comp_5  55.99834 -0.07691053  0.0975 -0.15796434  8.073955  6.898523
#> Nb_Comp_6  57.69592 -0.19968885  0.0975 -0.11400977  7.685022  6.835594
#> Nb_Comp_7  59.37953 -0.27722139  0.0975 -0.06462721  7.277359  6.770369
#> Nb_Comp_8  61.21213 -0.30602578  0.0975 -0.02255238  6.923057  6.736112
#> Nb_Comp_9  63.18426 -0.39920228  0.0975 -0.07134354  7.216690  6.730426
#> Nb_Comp_10 65.15982 -0.43743644  0.0975 -0.02732569  6.914340  6.725443
#>                 R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0         NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1  0.4675684 0.4675684   17.03781     19.76046  0.38248575 0.0975
#> Nb_Comp_2  0.5711905 0.5711905   13.72190     17.97925 -0.05525570 0.0975
#> Nb_Comp_3  0.6192438 0.6192438   12.18420     16.06943 -0.17107874 0.0975
#> Nb_Comp_4  0.6647841 0.6647841   10.72691     14.84877 -0.21869112 0.0975
#> Nb_Comp_5  0.6683426 0.6683426   10.61304     12.42138 -0.15796434 0.0975
#> Nb_Comp_6  0.6713681 0.6713681   10.51622     11.82303 -0.11400977 0.0975
#> Nb_Comp_7  0.6745039 0.6745039   10.41588     11.19586 -0.06462721 0.0975
#> Nb_Comp_8  0.6761508 0.6761508   10.36317     10.65078 -0.02255238 0.0975
#> Nb_Comp_9  0.6764242 0.6764242   10.35443     11.10252 -0.07134354 0.0975
#> Nb_Comp_10 0.6766638 0.6766638   10.34676     10.63737 -0.02732569 0.0975
#>            Q2cum_residY  AIC.std   DoF.dof sigmahat.dof   AIC.dof   BIC.dof
#> Nb_Comp_0            NA 96.63448  1.000000    0.8062287 0.6697018 0.6991787
#> Nb_Comp_1    0.38248575 77.83455  3.176360    0.5994089 0.4047616 0.4565153
#> Nb_Comp_2    0.34836456 72.69198  7.133559    0.5761829 0.4138120 0.5212090
#> Nb_Comp_3    0.23688359 70.76981  8.778329    0.5603634 0.4070516 0.5320535
#> Nb_Comp_4    0.06999681 68.56612  8.427874    0.5221703 0.3505594 0.4547689
#> Nb_Comp_5   -0.07691053 70.21393  9.308247    0.5285695 0.3666578 0.4845912
#> Nb_Comp_6   -0.19968885 71.91152  9.291931    0.5259794 0.3629363 0.4795121
#> Nb_Comp_7   -0.27722139 73.59512  9.756305    0.5284535 0.3702885 0.4938445
#> Nb_Comp_8   -0.30602578 75.42772 10.363948    0.5338475 0.3831339 0.5170783
#> Nb_Comp_9   -0.39920228 77.39986 10.732145    0.5378276 0.3920956 0.5328746
#> Nb_Comp_10  -0.43743644 79.37542 11.000000    0.5407500 0.3987417 0.5446065
#>             GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0  -3.605128         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1  -9.875081         2      0.5977015 0.3788984 0.4112998 -11.451340
#> Nb_Comp_2  -6.985517         3      0.5452615 0.3243383 0.3647862 -12.822703
#> Nb_Comp_3  -6.260610         4      0.5225859 0.3061986 0.3557368 -12.756838
#> Nb_Comp_4  -8.152986         5      0.4990184 0.2867496 0.3432131 -12.811575
#> Nb_Comp_5  -7.111583         6      0.5054709 0.3019556 0.3714754 -11.329638
#> Nb_Comp_6  -7.233043         7      0.5127450 0.3186757 0.4021333  -9.918688
#> Nb_Comp_7  -6.742195         8      0.5203986 0.3364668 0.4347156  -8.592770
#> Nb_Comp_8  -6.038372         9      0.5297842 0.3572181 0.4717708  -7.287834
#> Nb_Comp_9  -5.600238        10      0.5409503 0.3813021 0.5140048  -6.008747
#> Nb_Comp_10 -5.288422        11      0.5529032 0.4076026 0.5600977  -4.799453

data(pineNAX21)
bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
modele="pls-glm-family",family=Gamma(),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#> Warning: NaNs produced

#For Jackknife computations
kfolds2coeff(bbb2)
#>            [,1]        [,2]         [,3]        [,4]     [,5]       [,6]
#>  [1,] -13.59603 0.005527492  0.051731114 -0.17037188 1.842562 -0.3673061
#>  [2,] -12.10434 0.005601656  0.048250067 -0.16431469 1.575107 -0.3115874
#>  [3,] -17.25415 0.006626677  0.033302172 -0.23958076 1.244828 -0.2182543
#>  [4,] -13.90314 0.005390505  0.023566808 -0.18238338 1.381805 -0.2735833
#>  [5,] -13.22736 0.008306338 -0.001237529  0.09849343 1.870009 -0.3781797
#>  [6,] -15.91502 0.006648752  0.070207533 -0.17096513 2.388587 -0.5157931
#>  [7,] -14.43444 0.006069759  0.021492664 -0.17526927 1.398683 -0.3476223
#>  [8,] -13.75847 0.005783179  0.040195585 -0.13173309 1.755915 -0.3805948
#>  [9,] -15.92496 0.004867993  0.075448031 -0.27376532 2.722445 -0.5239134
#> [10,] -13.13506 0.005891946  0.022408429 -0.06523313 1.893981 -0.3958723
#>             [,7]        [,8]        [,9]     [,10]     [,11]
#>  [1,]  2.2025699 -0.40354813 -0.23686330 0.7630927 0.7114202
#>  [2,]  2.3268517 -0.02724398 -0.26178903 0.7205300 0.8077048
#>  [3,]  3.1899718  0.23072864 -0.32871996 1.5660330 0.7430596
#>  [4,]  2.0651784 -0.11531317 -0.14608822 1.4613947 0.5097030
#>  [5,] -0.6325947 -0.21364620 -0.13312074 0.3390221 0.2535890
#>  [6,]  2.2514518 -0.50219256 -0.15840483 0.3946695 0.8559487
#>  [7,]  1.9363240  0.18489078 -0.03431276 1.3825344 0.3195522
#>  [8,]  1.8085802 -0.20865945 -0.11855281 0.5217650 0.8232937
#>  [9,]  4.0002385 -0.71132728 -0.42473313 0.7363125 1.2460135
#> [10,]  0.8152612 -0.31129185 -0.08833544 0.6799973 0.8404652
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Family: Gamma 
#> Link function: inverse 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 56.60919 59.60220           NA     NA        NA                NA
#> Nb_Comp_1 39.08940 43.57892           NA 0.0975        NA                NA
#> Nb_Comp_2 37.36154 43.34757           NA 0.0975        NA                NA
#> Nb_Comp_3 36.81173 44.29427           NA 0.0975        NA                NA
#> Nb_Comp_4 36.53654 45.51559           NA 0.0975        NA                NA
#> Nb_Comp_5 37.24312 47.71867           NA 0.0975        NA                NA
#> Nb_Comp_6 38.18649 50.15855           NA 0.0975        NA                NA
#> Nb_Comp_7 39.35575 52.82432           NA 0.0975        NA                NA
#> Nb_Comp_8 40.86209 55.82716           NA 0.0975        NA                NA
#> Nb_Comp_9 42.80511 59.26669           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0       31.60805 20.800152        NA
#> Nb_Comp_1       17.30890 12.031518 0.4215659
#> Nb_Comp_2       17.10360  6.183372 0.7027247
#> Nb_Comp_3       15.78579  5.756462 0.7232490
#> Nb_Comp_4       13.49013  7.630460 0.6331536
#> Nb_Comp_5       13.56918  6.303455 0.6969515
#> Nb_Comp_6       14.02295  5.274716 0.7464097
#> Nb_Comp_7       15.05896  4.867806 0.7659726
#> Nb_Comp_8       15.28052  4.317488 0.7924300
#> Nb_Comp_9       15.19429  4.298593 0.7933384
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
PLS_lm(ypine,XpineNAX21,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> ____There are some NAs in X but not in Y____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE] < 10^{-12}
#> Warning only 9 components could thus be extracted
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#>                AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0 82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1 63.69250  0.35639805  0.0975  0.35639805 13.387018 11.099368
#> Nb_Comp_2 58.35228  0.28395028  0.0975 -0.11256611 12.348781  8.885823
#> Nb_Comp_3 56.36553  0.07664889  0.0975 -0.28950699 11.458331  7.874634
#> Nb_Comp_4 54.02416 -0.70355579  0.0975 -0.84497074 14.528469  6.903925
#> Nb_Comp_5 55.80450 -0.94905654  0.0975 -0.14411078  7.898855  6.858120
#> Nb_Comp_6 57.45753 -1.27568315  0.0975 -0.16758190  8.007417  6.786392
#> Nb_Comp_7 58.73951 -1.63309014  0.0975 -0.15705481  7.852227  6.640327
#> Nb_Comp_8 60.61227 -1.67907859  0.0975 -0.01746558  6.756304  6.614773
#> Nb_Comp_9 62.25948 -2.15165796  0.0975 -0.17639623  7.781594  6.544432
#>                R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0        NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1 0.4663804 0.4663804   17.07583     20.59526  0.35639805 0.0975
#> Nb_Comp_2 0.5728001 0.5728001   13.67040     18.99799 -0.11256611 0.0975
#> Nb_Comp_3 0.6214146 0.6214146   12.11473     17.62807 -0.28950699 0.0975
#> Nb_Comp_4 0.6680830 0.6680830   10.62135     22.35133 -0.84497074 0.0975
#> Nb_Comp_5 0.6702851 0.6702851   10.55088     12.15200 -0.14411078 0.0975
#> Nb_Comp_6 0.6737336 0.6737336   10.44053     12.31901 -0.16758190 0.0975
#> Nb_Comp_7 0.6807558 0.6807558   10.21581     12.08026 -0.15705481 0.0975
#> Nb_Comp_8 0.6819844 0.6819844   10.17650     10.39424 -0.01746558 0.0975
#> Nb_Comp_9 0.6853661 0.6853661   10.06828     11.97160 -0.17639623 0.0975
#>           Q2cum_residY  AIC.std DoF.dof sigmahat.dof AIC.dof BIC.dof GMDL.dof
#> Nb_Comp_0           NA 96.63448      NA           NA      NA      NA       NA
#> Nb_Comp_1   0.35639805 77.90810      NA           NA      NA      NA       NA
#> Nb_Comp_2   0.28395028 72.56787      NA           NA      NA      NA       NA
#> Nb_Comp_3   0.07664889 70.58113      NA           NA      NA      NA       NA
#> Nb_Comp_4  -0.70355579 68.23976      NA           NA      NA      NA       NA
#> Nb_Comp_5  -0.94905654 70.02009      NA           NA      NA      NA       NA
#> Nb_Comp_6  -1.27568315 71.67313      NA           NA      NA      NA       NA
#> Nb_Comp_7  -1.63309014 72.95511      NA           NA      NA      NA       NA
#> Nb_Comp_8  -1.67907859 74.82787      NA           NA      NA      NA       NA
#> Nb_Comp_9  -2.15165796 76.47507      NA           NA      NA      NA       NA
#>           DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1         2      0.5983679 0.3797438 0.4122175 -11.413749
#> Nb_Comp_2         3      0.5442372 0.3231208 0.3634169 -12.847656
#> Nb_Comp_3         4      0.5210941 0.3044529 0.3537087 -12.776843
#> Nb_Comp_4         5      0.4965569 0.2839276 0.3398355 -12.891035
#> Nb_Comp_5         6      0.5039885 0.3001871 0.3692997 -11.349498
#> Nb_Comp_6         7      0.5108963 0.3163819 0.3992388  -9.922119
#> Nb_Comp_7         8      0.5153766 0.3300041 0.4263658  -8.696873
#> Nb_Comp_8         9      0.5249910 0.3507834 0.4632727  -7.337679
#> Nb_Comp_9        10      0.5334234 0.3707649 0.4998004  -6.033403
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2"))



data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,modele="pls",verbose=FALSE)
summary(bbb)
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y        Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA          NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.8820701  0.0975  0.88207011 55.16721  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8703549  0.0975 -0.09934049 39.29316  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.7769175  0.0975 -0.72071678 19.04249   4.418081 0.9905556
#> Nb_Comp_4 33.76477 -0.4033063  0.0975 -5.29052595 27.79206   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -8.0844064  0.0975 -5.47357360 27.89615   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975          NA       NA   3.496074 0.9925265
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000002    0.7633343  0.9711321  1.1359501 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"

cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=12,verbose=FALSE)
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 12
cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
family=inverse.gaussian,K=12,verbose=FALSE)
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 12
cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 8    NA   NA   NA
#> 
#> [[1]][[4]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 5    NA   NA   NA
#> 
#> [[1]][[5]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 6    NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 5    NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 6    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[2]][[5]]
#>    [,1] [,2] [,3]
#> 4    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[2]][[6]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> 
cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=inverse.gaussian(),
K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 1    NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 2    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 7    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 4    NA   NA   NA
#> 
#> [[2]][[2]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 7    NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 5    NA   NA   NA
#> 
#> 
cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 3    NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[1]][[3]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 7    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[4]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[2]][[5]]
#>    [,1] [,2] [,3]
#> 4    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 2    NA   NA   NA
#> 
#> 
cv.plsRglm(object=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
family=inverse.gaussian(link = "1/mu^2"),K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 6    NA   NA   NA
#> 
#> [[1]][[4]]
#>    [,1] [,2] [,3]
#> 7    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[1]][[6]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 4    NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 9    NA   NA   NA
#> 
#> [[2]][[2]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 7    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 6    NA   NA   NA
#> 
#> 

bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=10,
modele="pls-glm-inverse.gaussian",keepcoeffs=TRUE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
#>               [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,] -3.674228e-05  0.0029015201  0.0001494670 -0.0047774817  2.022648e-05
#> [2,]  3.009054e-05  0.0001134832  0.0001030883  0.0001929214  1.230990e-04
#> [3,]  9.011167e-03 -0.0069812070 -0.0088750830 -0.0122093894 -8.839722e-03
#> [4,] -3.897643e-03  0.0044100933  0.0040305466  0.0034972243  4.051040e-03
#> [5,]  6.260532e-04  0.0004636377 -0.0004860292 -0.0021876776 -5.246432e-04
#>               [,6]          [,7]          [,8]
#> [1,]  0.0001705475  8.021070e-05  1.546294e-03
#> [2,]  0.0001003421  6.649647e-05  5.961471e-05
#> [3,] -0.0088865464 -8.924653e-03 -8.735091e-03
#> [4,]  0.0040099909  3.998541e-03  3.983696e-03
#> [5,] -0.0005175895 -5.523313e-04 -6.405018e-05
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> 
#> Family: inverse.gaussian 
#> Link function: 1/mu^2 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 81.67928 82.64909           NA     NA        NA                NA
#> Nb_Comp_1 49.90521 51.35993           NA 0.0975        NA                NA
#> Nb_Comp_2 31.06918 33.00881           NA 0.0975        NA                NA
#> Nb_Comp_3 28.40632 30.83085           NA 0.0975        NA                NA
#> Nb_Comp_4 27.08522 29.99466           NA 0.0975        NA                NA
#> Nb_Comp_5 28.46056 31.85490           NA 0.0975        NA                NA
#> Nb_Comp_6 29.68366 33.56292           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0   6.729783e-04 467.796667        NA
#> Nb_Comp_1   3.957680e-05  32.478677 0.9305710
#> Nb_Comp_2   7.009452e-06   6.020269 0.9871306
#> Nb_Comp_3   4.727777e-06   3.795855 0.9918857
#> Nb_Comp_4   3.584346e-06   2.699884 0.9942285
#> Nb_Comp_5   3.408069e-06   2.598572 0.9944451
#> Nb_Comp_6   3.195402e-06   2.492371 0.9946721
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm(yCornell,XCornell,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                AIC   Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205        NA      NA          NA        NA 467.796667        NA
#> Nb_Comp_1 53.15173 0.8966556  0.0975  0.89665563 48.344150  35.742486 0.9235940
#> Nb_Comp_2 41.08283 0.9175426  0.0975  0.20210989 28.518576  11.066606 0.9763431
#> Nb_Comp_3 32.06411 0.9399676  0.0975  0.27195907  8.056942   4.418081 0.9905556
#> Nb_Comp_4 33.76477 0.9197009  0.0975 -0.33759604  5.909608   4.309235 0.9907882
#> Nb_Comp_5 33.34373 0.9281373  0.0975  0.10506161  3.856500   3.521924 0.9924713
#> Nb_Comp_6 35.25533 0.9232562  0.0975 -0.06792167  3.761138   3.496074 0.9925265
#>           R2_residY  RSS_residY PRESS_residY   Q2_residY  LimQ2 Q2cum_residY
#> Nb_Comp_0        NA 11.00000000           NA          NA     NA           NA
#> Nb_Comp_1 0.9235940  0.84046633   1.13678803  0.89665563 0.0975    0.8966556
#> Nb_Comp_2 0.9763431  0.26022559   0.67059977  0.20210989 0.0975    0.9175426
#> Nb_Comp_3 0.9905556  0.10388893   0.18945488  0.27195907 0.0975    0.9399676
#> Nb_Comp_4 0.9907882  0.10132947   0.13896142 -0.33759604 0.0975    0.9197009
#> Nb_Comp_5 0.9924713  0.08281624   0.09068364  0.10506161 0.0975    0.9281373
#> Nb_Comp_6 0.9925265  0.08220840   0.08844125 -0.06792167 0.0975    0.9232562
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000002    0.7633343  0.9711321  1.1359501 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
rm(list=c("XCornell","yCornell","bbb","bbb2"))
# }
data(Cornell)
bbb <- cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,modele="pls")
#> 
#> Model: pls 
#> 
#> NK: 1 
#> Number of groups : 5 
#> 1 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
#> 2 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> Warning :  < 10^{-12}
#> Warning only 5 components could thus be extracted
#> ****________________________________________________****
#> 
#> 3 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
#> 4 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
#> 5 
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning :  < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ****________________________________________________****
#> 
summary(bbb)
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y       Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA         NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.9008515  0.0975  0.9008515 46.38133  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8690684  0.0975 -0.3205604 47.20011  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.6302084  0.0975 -1.8243125 31.25555   4.418081 0.9905556
#> Nb_Comp_4 33.76477 -0.8107407  0.0975 -3.8966512 21.63380   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -9.5546332  0.0975 -4.8289039 25.11812   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975         NA       NA   3.496074 0.9925265
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000002    0.7633343  0.9711321  1.1359501 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"

cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=12)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> NK: 1 
#> Leave One Out
#> Number of groups : 12 
#> 1 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 2 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 3 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 4 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 5 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 6 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 7 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 8 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 9 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 10 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 11 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> 12 
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 12

# \donttest{
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,
NK=2,random=TRUE,keepfolds=TRUE,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 4    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[1]][[3]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[1]][[5]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 6    NA   NA   NA
#> 
#> [[1]][[6]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 3    NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[2]][[2]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 7    NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[2]][[6]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> 

#Different ways of model specifications
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,
NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 5    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[1]][[3]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[1]][[4]]
#>    [,1] [,2] [,3]
#> 4    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[1]][[5]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 7    NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 9    NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 7    NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 3    NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[6]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> 
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian,
K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 7    NA   NA   NA
#> 
#> [[1]][[3]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[1]][[5]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 2    NA   NA   NA
#> 
#> [[1]][[6]]
#>    [,1] [,2] [,3]
#> 8    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[2]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[2]][[4]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[2]][[5]]
#>    [,1] [,2] [,3]
#> 9    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 11   NA   NA   NA
#> 
#> 
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),
K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[1]][[4]]
#>    [,1] [,2] [,3]
#> 3    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[2]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[2]][[3]]
#>    [,1] [,2] [,3]
#> 6    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 5    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 9    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> 
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),
K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 2    NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 3    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> [[1]][[6]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 6    NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 4    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[2]][[3]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 2    NA   NA   NA
#> 
#> [[2]][[4]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 7    NA   NA   NA
#> 
#> 

bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=10,
modele="pls-glm-gaussian",keepcoeffs=TRUE,verbose=FALSE)
bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",
family=gaussian(link=log),K=6,keepcoeffs=TRUE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
#>          [,1]        [,2]         [,3]        [,4]        [,5]         [,6]
#> [1,] 4.465364 -0.09104103 -0.011451364 -0.14597215 -0.03832459  0.069302213
#> [2,] 4.500596 -0.10304884 -0.036075886 -0.17189163 -0.06063713 -0.024692773
#> [3,] 4.471075 -0.08641298 -0.017500026 -0.14426712 -0.04722771  0.010149405
#> [4,] 4.456008 -0.04771799 -0.007331319 -0.07597615 -0.05095050 -0.005984889
#> [5,] 4.472561 -0.07675484 -0.023421939 -0.13048323 -0.05421782  0.062538452
#> [6,] 4.462804 -0.08690227 -0.008551162 -0.14521048 -0.04277721 -0.023997123
#>           [,7]       [,8]
#> [1,] 0.1526247 -0.1795973
#> [2,] 0.1388888 -0.3900981
#> [3,] 0.1631165 -0.2269663
#> [4,] 0.1852185 -0.2171953
#> [5,] 0.1556353 -0.1991052
#> [6,] 0.1758965 -0.1540098
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: log 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA        NA                NA
#> Nb_Comp_1 52.67938 54.13410           NA 0.0975        NA                NA
#> Nb_Comp_2 32.16524 34.10487           NA 0.0975        NA                NA
#> Nb_Comp_3 30.58789 33.01242           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      34.362913  34.362913 0.9265431
#> Nb_Comp_2       5.263520   5.263520 0.9887483
#> Nb_Comp_3       3.906676   3.906676 0.9916488
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(Y~.,data=Cornell,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning : 1 2 3 4 5 6 7 < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                AIC   Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205        NA      NA          NA        NA 467.796667        NA
#> Nb_Comp_1 53.15173 0.8966556  0.0975  0.89665563 48.344150  35.742486 0.9235940
#> Nb_Comp_2 41.08283 0.9175426  0.0975  0.20210989 28.518576  11.066606 0.9763431
#> Nb_Comp_3 32.06411 0.9399676  0.0975  0.27195907  8.056942   4.418081 0.9905556
#> Nb_Comp_4 33.76477 0.9197009  0.0975 -0.33759604  5.909608   4.309235 0.9907882
#> Nb_Comp_5 33.34373 0.9281373  0.0975  0.10506161  3.856500   3.521924 0.9924713
#> Nb_Comp_6 35.25533 0.9232562  0.0975 -0.06792167  3.761138   3.496074 0.9925265
#>           R2_residY  RSS_residY PRESS_residY   Q2_residY  LimQ2 Q2cum_residY
#> Nb_Comp_0        NA 11.00000000           NA          NA     NA           NA
#> Nb_Comp_1 0.9235940  0.84046633   1.13678803  0.89665563 0.0975    0.8966556
#> Nb_Comp_2 0.9763431  0.26022559   0.67059977  0.20210989 0.0975    0.9175426
#> Nb_Comp_3 0.9905556  0.10388893   0.18945488  0.27195907 0.0975    0.9399676
#> Nb_Comp_4 0.9907882  0.10132947   0.13896142 -0.33759604 0.0975    0.9197009
#> Nb_Comp_5 0.9924713  0.08281624   0.09068364  0.10506161 0.0975    0.9281373
#> Nb_Comp_6 0.9925265  0.08220840   0.08844125 -0.06792167 0.0975    0.9232562
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000001    0.7633342  0.9711319  1.1359499 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
rm(list=c("bbb","bbb2"))


data(pine)
bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",
family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",family=gaussian(),
K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>           [,1]         [,2]        [,3]        [,4]        [,5]        [,6]
#>  [1,] 8.827411 -0.003024637 -0.03212902  0.03318621 -0.61824497 0.116060970
#>  [2,] 8.381799 -0.002865575 -0.03081204  0.04562679 -0.66246828 0.145696363
#>  [3,] 8.958661 -0.003360521 -0.04519279  0.03007382 -0.39761466 0.110767333
#>  [4,] 7.413166 -0.002164980 -0.02798660 -0.02742069  0.03090166 0.004760803
#>  [5,] 8.039352 -0.002771071 -0.03364280  0.02339507 -0.43817893 0.093654003
#>  [6,] 8.462344 -0.002886290 -0.03918590  0.02001647 -0.44446226 0.094024062
#>  [7,] 9.313106 -0.003322881 -0.03783576  0.06761070 -0.49574386 0.108648773
#>  [8,] 7.212905 -0.002927104 -0.03104309  0.02585008 -0.47071133 0.102927075
#>  [9,] 7.709035 -0.002698856 -0.02905105  0.02620614 -0.46271909 0.120018938
#> [10,] 8.102813 -0.002955022 -0.03805164  0.02010760 -0.44712101 0.104991254
#>              [,7]       [,8]         [,9]      [,10]      [,11]
#>  [1,]  0.06231920 -0.1540749  0.042242967 -0.9398882 -0.4511634
#>  [2,] -0.09954001 -0.3812528  0.072218736 -0.9330554 -0.3424166
#>  [3,] -0.11821767 -0.2670741 -0.010633715 -0.6872869 -0.3434467
#>  [4,]  1.00678066 -0.9941689 -0.084677908 -1.1702028  0.0157466
#>  [5,]  0.10257199 -0.3614453  0.017161596 -0.8085061 -0.2766541
#>  [6,]  0.19999283 -0.4774163  0.034196834 -0.8977215 -0.2400707
#>  [7,] -0.50417691  0.1344707  0.034509256 -0.8107578 -0.7153230
#>  [8,]  0.07726327  0.2126788  0.033672355 -0.8534570 -0.2883699
#>  [9,]  0.14461054 -0.2346029  0.037680525 -1.0574493 -0.3637227
#> [10,]  0.15390012 -0.2765310 -0.008502297 -0.7278936 -0.2552828
boxplot(kfolds2coeff(bbb)[,1])


kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
summary(bbb)
#> ____************************************************____
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0  82.41888 85.41190           NA     NA        NA                NA
#> Nb_Comp_1  63.61896 68.10848           NA 0.0975        NA                NA
#> Nb_Comp_2  54.15489 60.14092           NA 0.0975        NA                NA
#> Nb_Comp_3  53.47303 60.95556           NA 0.0975        NA                NA
#> Nb_Comp_4  54.83398 63.81302           NA 0.0975        NA                NA
#> Nb_Comp_5  56.32757 66.80312           NA 0.0975        NA                NA
#> Nb_Comp_6  57.45220 69.42426           NA 0.0975        NA                NA
#> Nb_Comp_7  59.31417 72.78274           NA 0.0975        NA                NA
#> Nb_Comp_8  61.20356 76.16863           NA 0.0975        NA                NA
#> Nb_Comp_9  63.16270 79.62429           NA 0.0975        NA                NA
#> Nb_Comp_10 65.15982 83.11791           NA 0.0975        NA                NA
#>            Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0       20.800152 20.800152        NA
#> Nb_Comp_1       11.074659 11.074659 0.4675684
#> Nb_Comp_2        7.824528  7.824528 0.6238235
#> Nb_Comp_3        7.213793  7.213793 0.6531855
#> Nb_Comp_4        7.075441  7.075441 0.6598370
#> Nb_Comp_5        6.967693  6.967693 0.6650172
#> Nb_Comp_6        6.785296  6.785296 0.6737862
#> Nb_Comp_7        6.756973  6.756973 0.6751479
#> Nb_Comp_8        6.734363  6.734363 0.6762349
#> Nb_Comp_9        6.726030  6.726030 0.6766355
#> Nb_Comp_10       6.725443  6.725443 0.6766638
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(x11~.,data=pine,nt=10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                 AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0  82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1  63.61896  0.38248575  0.0975  0.38248575 12.844390 11.074659
#> Nb_Comp_2  58.47638  0.34836456  0.0975 -0.05525570 11.686597  8.919303
#> Nb_Comp_3  56.55421  0.23688359  0.0975 -0.17107874 10.445206  7.919786
#> Nb_Comp_4  54.35053  0.06999681  0.0975 -0.21869112  9.651773  6.972542
#> Nb_Comp_5  55.99834 -0.07691053  0.0975 -0.15796434  8.073955  6.898523
#> Nb_Comp_6  57.69592 -0.19968885  0.0975 -0.11400977  7.685022  6.835594
#> Nb_Comp_7  59.37953 -0.27722139  0.0975 -0.06462721  7.277359  6.770369
#> Nb_Comp_8  61.21213 -0.30602578  0.0975 -0.02255238  6.923057  6.736112
#> Nb_Comp_9  63.18426 -0.39920228  0.0975 -0.07134354  7.216690  6.730426
#> Nb_Comp_10 65.15982 -0.43743644  0.0975 -0.02732569  6.914340  6.725443
#>                 R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0         NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1  0.4675684 0.4675684   17.03781     19.76046  0.38248575 0.0975
#> Nb_Comp_2  0.5711905 0.5711905   13.72190     17.97925 -0.05525570 0.0975
#> Nb_Comp_3  0.6192438 0.6192438   12.18420     16.06943 -0.17107874 0.0975
#> Nb_Comp_4  0.6647841 0.6647841   10.72691     14.84877 -0.21869112 0.0975
#> Nb_Comp_5  0.6683426 0.6683426   10.61304     12.42138 -0.15796434 0.0975
#> Nb_Comp_6  0.6713681 0.6713681   10.51622     11.82303 -0.11400977 0.0975
#> Nb_Comp_7  0.6745039 0.6745039   10.41588     11.19586 -0.06462721 0.0975
#> Nb_Comp_8  0.6761508 0.6761508   10.36317     10.65078 -0.02255238 0.0975
#> Nb_Comp_9  0.6764242 0.6764242   10.35443     11.10252 -0.07134354 0.0975
#> Nb_Comp_10 0.6766638 0.6766638   10.34676     10.63737 -0.02732569 0.0975
#>            Q2cum_residY  AIC.std   DoF.dof sigmahat.dof   AIC.dof   BIC.dof
#> Nb_Comp_0            NA 96.63448  1.000000    0.8062287 0.6697018 0.6991787
#> Nb_Comp_1    0.38248575 77.83455  3.176360    0.5994089 0.4047616 0.4565153
#> Nb_Comp_2    0.34836456 72.69198  7.133559    0.5761829 0.4138120 0.5212090
#> Nb_Comp_3    0.23688359 70.76981  8.778329    0.5603634 0.4070516 0.5320535
#> Nb_Comp_4    0.06999681 68.56612  8.427874    0.5221703 0.3505594 0.4547689
#> Nb_Comp_5   -0.07691053 70.21393  9.308247    0.5285695 0.3666578 0.4845912
#> Nb_Comp_6   -0.19968885 71.91152  9.291931    0.5259794 0.3629363 0.4795121
#> Nb_Comp_7   -0.27722139 73.59512  9.756305    0.5284535 0.3702885 0.4938445
#> Nb_Comp_8   -0.30602578 75.42772 10.363948    0.5338475 0.3831339 0.5170783
#> Nb_Comp_9   -0.39920228 77.39986 10.732146    0.5378276 0.3920957 0.5328746
#> Nb_Comp_10  -0.43743644 79.37542 11.000000    0.5407500 0.3987417 0.5446065
#>             GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0  -3.605128         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1  -9.875081         2      0.5977015 0.3788984 0.4112998 -11.451340
#> Nb_Comp_2  -6.985517         3      0.5452615 0.3243383 0.3647862 -12.822703
#> Nb_Comp_3  -6.260610         4      0.5225859 0.3061986 0.3557368 -12.756838
#> Nb_Comp_4  -8.152986         5      0.4990184 0.2867496 0.3432131 -12.811575
#> Nb_Comp_5  -7.111583         6      0.5054709 0.3019556 0.3714754 -11.329638
#> Nb_Comp_6  -7.233043         7      0.5127450 0.3186757 0.4021333  -9.918688
#> Nb_Comp_7  -6.742195         8      0.5203986 0.3364668 0.4347156  -8.592770
#> Nb_Comp_8  -6.038372         9      0.5297842 0.3572181 0.4717708  -7.287834
#> Nb_Comp_9  -5.600237        10      0.5409503 0.3813021 0.5140048  -6.008747
#> Nb_Comp_10 -5.288422        11      0.5529032 0.4076026 0.5600977  -4.799453

data(pineNAX21)
bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
modele="pls-glm-family",family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
modele="pls-glm-gaussian",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
#>           [,1]          [,2]        [,3]         [,4]       [,5]       [,6]
#>  [1,] 3.320094  0.0003034149  0.02340033 -0.222624357  2.1105469 -0.3999920
#>  [2,] 8.860816 -0.0043033126 -0.03665880  0.093040486 -0.7851291  0.1805536
#>  [3,] 7.615754 -0.0034017330 -0.03188923  0.046195833 -0.4889310  0.1347685
#>  [4,] 7.171507 -0.0015893703 -0.03225535  0.043853427 -0.6011096  0.1027683
#>  [5,] 7.816351 -0.0028312877 -0.03366530  0.021577097 -0.5411489  0.1152638
#>  [6,] 6.165179 -0.0031097763 -0.03168662 -0.008606572 -0.5312360  0.1469865
#>  [7,] 7.815241 -0.0030006293 -0.03160637  0.064130903 -0.6043211  0.1182370
#>  [8,] 7.891825 -0.0030669630 -0.03669431  0.056292297 -0.4957551  0.1131661
#>  [9,] 9.410062 -0.0043633297 -0.03393477  0.151454704 -0.4869584  0.1187660
#> [10,] 7.488769 -0.0028843195 -0.02879965  0.032286072 -0.5071476  0.1122509
#>              [,7]        [,8]         [,9]      [,10]      [,11]
#>  [1,]  4.91899618 -4.58295354 -0.328270297 -3.4760604  1.8028654
#>  [2,] -0.84650686  0.78498808  0.044412465 -0.5723738 -0.8096367
#>  [3,] -0.13694666  0.03819981  0.004212034 -0.9081770 -0.3631493
#>  [4,] -0.26399426 -0.47486794  0.138358585 -1.0400658 -0.7125486
#>  [5,]  0.18724953 -0.02069456  0.007544203 -0.8153538 -0.3566206
#>  [6,]  0.58572883  0.41305648 -0.134905394 -0.3445386 -0.5380991
#>  [7,] -0.53232551  0.02030175  0.082162364 -0.8060750 -0.5220032
#>  [8,] -0.35838855 -0.09109065  0.057031529 -0.9074141 -0.4935589
#>  [9,] -1.76116016  1.48609833 -0.011080337 -0.1975556 -1.6475023
#> [10,] -0.02893316 -0.15482830  0.023685857 -0.7421772 -0.6014031
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.41888 85.41190           NA     NA        NA                NA
#> Nb_Comp_1 63.90814 68.39766           NA 0.0975        NA                NA
#> Nb_Comp_2 54.06295 60.04898           NA 0.0975        NA                NA
#> Nb_Comp_3 53.77276 61.25530           NA 0.0975        NA                NA
#> Nb_Comp_4 55.18223 64.16127           NA 0.0975        NA                NA
#> Nb_Comp_5 56.53963 67.01518           NA 0.0975        NA                NA
#> Nb_Comp_6 57.73540 69.70746           NA 0.0975        NA                NA
#> Nb_Comp_7 59.46634 72.93491           NA 0.0975        NA                NA
#> Nb_Comp_8 60.79943 75.76451           NA 0.0975        NA                NA
#> Nb_Comp_9 62.14147 78.60305           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0      20.800152 20.800152        NA
#> Nb_Comp_1      11.172133 11.172133 0.4628821
#> Nb_Comp_2       7.802760  7.802760 0.6248700
#> Nb_Comp_3       7.279614  7.279614 0.6500211
#> Nb_Comp_4       7.150504  7.150504 0.6562283
#> Nb_Comp_5       7.012612  7.012612 0.6628577
#> Nb_Comp_6       6.843775  6.843775 0.6709747
#> Nb_Comp_7       6.788203  6.788203 0.6736465
#> Nb_Comp_8       6.652395  6.652395 0.6801757
#> Nb_Comp_9       6.521071  6.521071 0.6864893
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(x11~.,data=pineNAX21,nt=10,typeVC="standard")$InfCrit
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> ____There are some NAs in X but not in Y____
#> ____TypeVC____ standard ____
#> ____TypeVC____ standard ____unknown____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE] < 10^{-12}
#> Warning only 9 components could thus be extracted
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#>                AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0 82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1 63.69250  0.35639805  0.0975  0.35639805 13.387018 11.099368
#> Nb_Comp_2 58.35228  0.28395028  0.0975 -0.11256611 12.348781  8.885823
#> Nb_Comp_3 56.36553  0.07664889  0.0975 -0.28950699 11.458331  7.874634
#> Nb_Comp_4 54.02416 -0.70355579  0.0975 -0.84497074 14.528469  6.903925
#> Nb_Comp_5 55.80450 -0.94905654  0.0975 -0.14411078  7.898855  6.858120
#> Nb_Comp_6 57.45753 -1.27568315  0.0975 -0.16758190  8.007417  6.786392
#> Nb_Comp_7 58.73951 -1.63309014  0.0975 -0.15705481  7.852227  6.640327
#> Nb_Comp_8 60.61227 -1.67907859  0.0975 -0.01746558  6.756304  6.614773
#> Nb_Comp_9 62.25948 -2.15165796  0.0975 -0.17639623  7.781594  6.544432
#>                R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0        NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1 0.4663804 0.4663804   17.07583     20.59526  0.35639805 0.0975
#> Nb_Comp_2 0.5728001 0.5728001   13.67040     18.99799 -0.11256611 0.0975
#> Nb_Comp_3 0.6214146 0.6214146   12.11473     17.62807 -0.28950699 0.0975
#> Nb_Comp_4 0.6680830 0.6680830   10.62135     22.35133 -0.84497074 0.0975
#> Nb_Comp_5 0.6702851 0.6702851   10.55088     12.15200 -0.14411078 0.0975
#> Nb_Comp_6 0.6737336 0.6737336   10.44053     12.31901 -0.16758190 0.0975
#> Nb_Comp_7 0.6807558 0.6807558   10.21581     12.08026 -0.15705481 0.0975
#> Nb_Comp_8 0.6819844 0.6819844   10.17650     10.39424 -0.01746558 0.0975
#> Nb_Comp_9 0.6853661 0.6853661   10.06828     11.97160 -0.17639623 0.0975
#>           Q2cum_residY  AIC.std DoF.dof sigmahat.dof AIC.dof BIC.dof GMDL.dof
#> Nb_Comp_0           NA 96.63448      NA           NA      NA      NA       NA
#> Nb_Comp_1   0.35639805 77.90810      NA           NA      NA      NA       NA
#> Nb_Comp_2   0.28395028 72.56787      NA           NA      NA      NA       NA
#> Nb_Comp_3   0.07664889 70.58113      NA           NA      NA      NA       NA
#> Nb_Comp_4  -0.70355579 68.23976      NA           NA      NA      NA       NA
#> Nb_Comp_5  -0.94905654 70.02009      NA           NA      NA      NA       NA
#> Nb_Comp_6  -1.27568315 71.67313      NA           NA      NA      NA       NA
#> Nb_Comp_7  -1.63309014 72.95511      NA           NA      NA      NA       NA
#> Nb_Comp_8  -1.67907859 74.82787      NA           NA      NA      NA       NA
#> Nb_Comp_9  -2.15165796 76.47507      NA           NA      NA      NA       NA
#>           DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1         2      0.5983679 0.3797438 0.4122175 -11.413749
#> Nb_Comp_2         3      0.5442372 0.3231208 0.3634169 -12.847656
#> Nb_Comp_3         4      0.5210941 0.3044529 0.3537087 -12.776843
#> Nb_Comp_4         5      0.4965569 0.2839276 0.3398355 -12.891035
#> Nb_Comp_5         6      0.5039885 0.3001871 0.3692997 -11.349498
#> Nb_Comp_6         7      0.5108963 0.3163819 0.3992388  -9.922119
#> Nb_Comp_7         8      0.5153766 0.3300041 0.4263658  -8.696873
#> Nb_Comp_8         9      0.5249910 0.3507834 0.4632727  -7.337679
#> Nb_Comp_9        10      0.5334234 0.3707649 0.4998004  -6.033403
rm(list=c("bbb","bbb2"))


data(aze_compl)
bbb <- cv.plsRglm(y~.,data=aze_compl,nt=10,K=10,modele="pls",
keepcoeffs=TRUE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>            [,1]        [,2]      [,3]       [,4]      [,5]       [,6]
#>  [1,] 0.4720148 -0.18921733 0.3117505 -0.2479646 0.2701179 0.11141452
#>  [2,] 0.1711700 -0.16636435 0.5022693 -0.1793904 0.4002111 0.05562376
#>  [3,] 0.4535198 -0.11975322 0.4236443 -0.1775557 0.2044521 0.06628416
#>  [4,] 0.2564124 -0.27608027 0.4933865 -0.1385890 0.3804256 0.10334714
#>  [5,] 0.2385741 -0.06422092 0.4577761 -0.1904053 0.2261762 0.13472776
#>  [6,] 0.2164312 -0.15323727 0.5928721 -0.2073372 0.2285467 0.12942366
#>  [7,] 0.3287707 -0.02781054 0.3349072 -0.1713153 0.2506213 0.02721249
#>  [8,] 0.2665835 -0.14027032 0.5268426 -0.2403808 0.3338911 0.11511991
#>  [9,] 0.2832101 -0.13135471 0.3974422 -0.1411602 0.2226834 0.17670768
#> [10,] 0.4421650 -0.12126484 0.4628998 -0.1914751 0.1520873 0.10008102
#>              [,7]        [,8]       [,9]        [,10]        [,11]        [,12]
#>  [1,] -0.10691664  0.16014049 -0.2125185  0.055643198 -0.008628953  0.087185137
#>  [2,] -0.04332097 -0.10134240 -0.2144580  0.067479370 -0.055989742 -0.050765888
#>  [3,] -0.08096434 -0.02156994 -0.1851704  0.087831997 -0.065191410 -0.001392207
#>  [4,] -0.03374670 -0.04031349 -0.3915468  0.061955068 -0.230966792  0.118821250
#>  [5,]  0.01112923 -0.08867503 -0.1151823  0.060439076 -0.074725936  0.038721524
#>  [6,] -0.05148663 -0.09335017 -0.2626375  0.009312704 -0.004776961  0.129024818
#>  [7,] -0.03576488  0.03662383 -0.1751392  0.052725198 -0.145109050  0.029481879
#>  [8,] -0.04541675  0.09679495 -0.1998142  0.024789710 -0.083069169  0.107257474
#>  [9,] -0.05901319  0.07777635 -0.1529209 -0.008701909 -0.152317948  0.114393720
#> [10,] -0.05375944  0.06166435 -0.1976978  0.051212795 -0.158988869  0.023802369
#>             [,13]      [,14]       [,15]       [,16]        [,17]     [,18]
#>  [1,] -0.13730479 0.08730486 0.090055403 0.101725293  0.046686558 0.2598979
#>  [2,] -0.18726936 0.16361149 0.003950767 0.127103677  0.017956953 0.2449673
#>  [3,] -0.11066851 0.13072128 0.096725673 0.045505228 -0.094775043 0.2129385
#>  [4,] -0.09517085 0.11354397 0.024545665 0.152723668  0.120022276 0.1179847
#>  [5,] -0.09561318 0.17510759 0.166450281 0.046255708  0.021027950 0.2088890
#>  [6,] -0.18359290 0.05455750 0.231373661 0.010849203 -0.004709846 0.3075805
#>  [7,] -0.15861236 0.10732489 0.165768812 0.070528269  0.040219947 0.2191750
#>  [8,] -0.11745752 0.12749012 0.123739099 0.061028238 -0.043563708 0.2767255
#>  [9,] -0.06029285 0.02481337 0.118837306 0.005175933 -0.052303082 0.3126029
#> [10,] -0.16534441 0.07022517 0.106346292 0.036746957  0.053761965 0.1826134
#>             [,19]      [,20]       [,21]        [,22]      [,23]       [,24]
#>  [1,]  0.05489183 0.04105146 -0.09982589 -0.008548760 0.14009828 -0.14218442
#>  [2,]  0.01526206 0.07296655 -0.07726314  0.069102798 0.18090255 -0.20634371
#>  [3,]  0.04030061 0.12321151 -0.12591167  0.094725156 0.18303355 -0.19111314
#>  [4,]  0.02253328 0.14467030 -0.07756995  0.159450058 0.26964128 -0.15335224
#>  [5,] -0.05441518 0.07470739 -0.16463246  0.098655109 0.07089352 -0.08800989
#>  [6,] -0.07826768 0.14844510 -0.05355764  0.032822849 0.20318608 -0.12872370
#>  [7,] -0.09326496 0.03187084 -0.13810902  0.059039405 0.19877376 -0.10801368
#>  [8,]  0.05227711 0.05856077 -0.09481483 -0.008473124 0.29002600 -0.18715780
#>  [9,]  0.08077597 0.06733765 -0.18276062  0.097751864 0.07240900 -0.04572254
#> [10,]  0.04618316 0.03247623 -0.19632683  0.102724284 0.18019885 -0.09224520
#>            [,25]      [,26]     [,27]      [,28]       [,29]        [,30]
#>  [1,] -0.1682622 -0.2757798 0.1564560 0.15916213 -0.15533252 -0.007424773
#>  [2,] -0.2177099 -0.2083308 0.1880946 0.16809465 -0.02933524  0.045575536
#>  [3,] -0.1764635 -0.3151767 0.1791981 0.19351771 -0.11894254  0.057250941
#>  [4,] -0.1937295 -0.2422658 0.1683737 0.27035747 -0.23986454 -0.028881509
#>  [5,] -0.1095170 -0.2472249 0.1745536 0.14954646 -0.04875785 -0.019105746
#>  [6,] -0.1082968 -0.3431003 0.1806241 0.18953853 -0.12064895 -0.001275921
#>  [7,] -0.1930979 -0.2927788 0.1525689 0.18124672 -0.05889064 -0.078729524
#>  [8,] -0.2731259 -0.2940052 0.2229317 0.09989895 -0.09095350  0.075912339
#>  [9,] -0.2139054 -0.2908271 0.1879460 0.26790441 -0.10818269  0.027191788
#> [10,] -0.1461844 -0.3013340 0.2342328 0.24324540 -0.06040999  0.012626419
#>              [,31]        [,32]      [,33]        [,34]
#>  [1,]  0.072724076 -0.024300839 -0.4241927 -0.041183813
#>  [2,]  0.265655122 -0.193880129 -0.3599459  0.059396124
#>  [3,]  0.117953111 -0.024440221 -0.4025196 -0.054336434
#>  [4,]  0.141764207 -0.082363842 -0.2219423 -0.060927469
#>  [5,]  0.129499168 -0.073962128 -0.4984025  0.035810900
#>  [6,]  0.006096966 -0.079120257 -0.2584478  0.001626410
#>  [7,]  0.089672323 -0.007992422 -0.1951715 -0.010292474
#>  [8,] -0.011769041 -0.009759173 -0.4453213 -0.004802358
#>  [9,]  0.192037167 -0.010225424 -0.5112457 -0.047196818
#> [10,]  0.148758244  0.029029514 -0.5018406 -0.034286216
bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=3,K=10,
modele="pls-glm-family",family=binomial(probit),keepcoeffs=TRUE,verbose=FALSE)
bbb2 <- cv.plsRglm(y~.,data=aze_compl,nt=3,K=10,
modele="pls-glm-logistic",keepcoeffs=TRUE,verbose=FALSE)
summary(bbb,MClassed=TRUE)
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC MissClassed CV_MissClassed       Q2cum_Y LimQ2_Y       Q2_Y
#> Nb_Comp_0  154.6179          49             NA            NA      NA         NA
#> Nb_Comp_1  126.4083          27             49    -0.1938111  0.0975 -0.1938111
#> Nb_Comp_2  119.3375          25             43    -0.8504107  0.0975 -0.5500030
#> Nb_Comp_3  114.2313          27             46    -2.4900820  0.0975 -0.8861121
#> Nb_Comp_4  112.3463          23             48    -6.8145995  0.0975 -1.2390877
#> Nb_Comp_5  113.2362          22             44   -17.2537261  0.0975 -1.3358492
#> Nb_Comp_6  114.7620          21             45   -42.6099196  0.0975 -1.3890969
#> Nb_Comp_7  116.5264          20             46  -102.8477774  0.0975 -1.3812880
#> Nb_Comp_8  118.4601          20             46  -245.6570159  0.0975 -1.3751786
#> Nb_Comp_9  120.4452          19             45  -584.6393179  0.0975 -1.3743063
#> Nb_Comp_10 122.4395          19             45 -1390.1500229  0.0975 -1.3754382
#>             PRESS_Y    RSS_Y      R2_Y  AIC.std  DoF.dof sigmahat.dof   AIC.dof
#> Nb_Comp_0        NA 25.91346        NA 298.1344  1.00000    0.5015845 0.2540061
#> Nb_Comp_1  30.93578 19.38086 0.2520929 269.9248 22.55372    0.4848429 0.2883114
#> Nb_Comp_2  30.04039 17.76209 0.3145613 262.8540 27.31542    0.4781670 0.2908950
#> Nb_Comp_3  33.50129 16.58896 0.3598323 257.7478 30.52370    0.4719550 0.2902572
#> Nb_Comp_4  37.14414 15.98071 0.3833049 255.8628 34.00000    0.4744263 0.3008285
#> Nb_Comp_5  37.32852 15.81104 0.3898523 256.7527 34.00000    0.4719012 0.2976347
#> Nb_Comp_6  37.77411 15.73910 0.3926285 258.2785 34.00000    0.4708264 0.2962804
#> Nb_Comp_7  37.47933 15.70350 0.3940024 260.0429 33.71066    0.4693382 0.2937976
#> Nb_Comp_8  37.29861 15.69348 0.3943888 261.9766 34.00000    0.4701436 0.2954217
#> Nb_Comp_9  37.26113 15.69123 0.3944758 263.9617 33.87284    0.4696894 0.2945815
#> Nb_Comp_10 37.27354 15.69037 0.3945088 265.9560 34.00000    0.4700970 0.2953632
#>              BIC.dof  GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive
#> Nb_Comp_0  0.2604032 -67.17645         1      0.5015845 0.2540061 0.2604032
#> Nb_Comp_1  0.4231184 -53.56607         2      0.4358996 0.1936625 0.2033251
#> Nb_Comp_2  0.4496983 -52.42272         3      0.4193593 0.1809352 0.1943501
#> Nb_Comp_3  0.4631316 -51.93343         4      0.4072955 0.1722700 0.1891422
#> Nb_Comp_4  0.4954133 -50.37079         5      0.4017727 0.1691819 0.1897041
#> Nb_Comp_5  0.4901536 -50.65724         6      0.4016679 0.1706451 0.1952588
#> Nb_Comp_6  0.4879234 -50.78005         7      0.4028135 0.1731800 0.2020601
#> Nb_Comp_7  0.4826103 -51.05525         8      0.4044479 0.1761610 0.2094352
#> Nb_Comp_8  0.4865092 -50.85833         9      0.4064413 0.1794902 0.2172936
#> Nb_Comp_9  0.4845867 -50.95616        10      0.4085682 0.1829787 0.2254232
#> Nb_Comp_10 0.4864128 -50.86368        11      0.4107477 0.1865584 0.2337468
#>            GMDL.naive
#> Nb_Comp_0   -67.17645
#> Nb_Comp_1   -79.67755
#> Nb_Comp_2   -81.93501
#> Nb_Comp_3   -83.31503
#> Nb_Comp_4   -83.23369
#> Nb_Comp_5   -81.93513
#> Nb_Comp_6   -80.42345
#> Nb_Comp_7   -78.87607
#> Nb_Comp_8   -77.31942
#> Nb_Comp_9   -75.80069
#> Nb_Comp_10  -74.33325
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"
summary(bbb2,MClassed=TRUE)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC MissClassed CV_MissClassed Q2Chisqcum_Y  limQ2
#> Nb_Comp_0 145.8283 148.4727          49             NA           NA     NA
#> Nb_Comp_1 118.1398 123.4285          28             NA           NA 0.0975
#> Nb_Comp_2 109.9553 117.8885          26             NA           NA 0.0975
#> Nb_Comp_3 105.1591 115.7366          22             NA           NA 0.0975
#>           Q2Chisq_Y PREChi2_Pearson_Y Chi2_Pearson_Y    RSS_Y      R2_Y
#> Nb_Comp_0        NA                NA      104.00000 25.91346        NA
#> Nb_Comp_1        NA                NA      100.53823 19.32272 0.2543365
#> Nb_Comp_2        NA                NA       99.17955 17.33735 0.3309519
#> Nb_Comp_3        NA                NA      123.37836 15.58198 0.3986915
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
kfolds2coeff(bbb2)
#>             [,1]        [,2]     [,3]       [,4]      [,5]        [,6]
#>  [1,] -0.5995726 -0.44472026 2.428467 -0.4996669 1.0112817  0.11544997
#>  [2,] -2.1266385 -0.59106118 2.578529 -0.3261298 1.4079666 -0.16507148
#>  [3,] -1.2138421 -1.13061925 3.296518 -0.7557706 1.0228252 -0.17880376
#>  [4,] -0.2422708  0.09325782 1.820126 -0.4709656 0.7629987 -0.14759472
#>  [5,] -0.4873652 -0.69065602 1.917047 -0.3012152 0.9766570  0.14468032
#>  [6,] -1.9912346 -0.05482392 2.556544 -0.4290094 1.4705423  0.22713821
#>  [7,] -1.1349518 -0.40726347 3.471090 -0.3938454 1.2294665 -0.10738233
#>  [8,] -1.2319146 -0.49553591 2.484675 -0.5471586 0.7059237  0.46257661
#>  [9,] -0.7400858 -0.43309053 1.850899 -0.2687466 0.6981307 -0.09826402
#> [10,] -0.9686851 -0.44031448 1.664519 -0.3596945 0.8024983  0.17024152
#>              [,7]        [,8]       [,9]       [,10]       [,11]        [,12]
#>  [1,] -0.35155749  0.27098144 -0.5643713 -0.13662260 -0.79103105  0.123831336
#>  [2,] -0.33875920 -0.21505556 -0.4924509  0.24099119 -0.52709571 -0.390424106
#>  [3,] -1.24122353  0.24140942 -0.8131511  0.39364901 -0.57020002  0.313976921
#>  [4,] -0.08958657  0.14840070 -0.4399510 -0.08706155 -0.77643530 -0.325078641
#>  [5,] -0.49416532  0.25744088 -0.7317188 -0.09530154 -0.65071265  0.178121483
#>  [6,] -0.83712540  0.39844834 -1.1228767  0.19861598 -0.98776848 -0.001330256
#>  [7,] -0.61232257  0.35391900 -0.8306103  0.69407364 -0.76966979 -0.473317800
#>  [8,] -0.44514582  0.09214555 -0.5011814 -0.13459763 -0.08211351  0.005477666
#>  [9,] -0.78752228 -0.02973437 -0.6036081 -0.04581873 -1.17048958  0.314363906
#> [10,] -0.70399320  0.71086730 -0.8258186  0.31612495 -0.62690471 -0.155449632
#>             [,13]     [,14]     [,15]     [,16]       [,17]     [,18]
#>  [1,] -0.90478069 0.9426560 0.9392174 0.4051812 -0.09400762 0.8587024
#>  [2,] -1.13364345 0.3755652 1.0299148 0.7411723  0.34159334 1.1516584
#>  [3,] -0.33411674 1.2884428 1.5446762 0.2367645 -0.15036339 1.2332013
#>  [4,] -0.15739680 0.8114167 1.3363416 0.3278741  0.28474462 1.1528850
#>  [5,] -0.83464052 0.9477550 0.8904600 0.5117776 -0.20822990 1.3105120
#>  [6,] -0.68377781 1.8491355 1.0079011 0.8513415  0.10314175 0.3888767
#>  [7,] -0.47614510 0.4417827 1.6202618 1.1324293 -0.58775337 1.2921074
#>  [8,]  0.09387891 0.1994558 1.0025869 0.2397655  0.13102599 1.6158143
#>  [9,]  0.06992107 0.8520815 1.2293568 0.1859066 -0.32419226 1.0448063
#> [10,] -0.55783037 0.4099192 1.0839515 0.7583290  0.22659004 0.6263444
#>             [,19]      [,20]      [,21]       [,22]     [,23]      [,24]
#>  [1,]  0.57236408  0.4614380 -0.6054459  0.53847881 0.8571912 -0.5698025
#>  [2,]  0.03679744  1.0768659 -0.9431263 -0.04749470 0.6459030 -0.2986934
#>  [3,] -0.19951780  1.3470703 -1.0062255 -0.07932805 1.0137349 -0.7226137
#>  [4,]  0.06485926 -0.1278838 -0.8178764  0.19150276 0.5506307 -0.9877573
#>  [5,]  0.34562437  0.4725870 -0.6030022  0.14711502 0.7539422 -0.3819852
#>  [6,] -0.33445395  0.8292700 -1.0593607  0.62429707 0.1838399 -0.3435331
#>  [7,]  0.26392703  1.0854399 -0.7062498 -0.28945809 1.1722975 -2.1293273
#>  [8,]  0.15152775  0.9261077 -0.9191057  0.34310697 0.5278247 -0.3552084
#>  [9,]  0.11838814  0.2778435 -0.7441412  0.20325502 0.8960014 -0.4272961
#> [10,]  0.26902944  0.8003778 -0.9202317  0.64246105 0.9081169 -0.7145389
#>            [,25]      [,26]     [,27]     [,28]      [,29]       [,30]
#>  [1,] -1.5753954 -1.1889646 0.5312176 0.6521445 -1.1081268  0.12881035
#>  [2,] -1.0306700 -1.2765836 0.6486454 1.1374112 -0.4380608  0.07143747
#>  [3,] -0.8429597 -1.5521500 0.2745060 0.8015847 -1.5620342 -0.47561243
#>  [4,] -0.7560257 -1.2734658 0.6606161 1.1740409 -0.4691185 -0.29698463
#>  [5,] -1.0094902 -1.5021999 0.4083606 0.3766212 -0.6872039 -0.49113345
#>  [6,] -1.5086798 -1.6641615 0.6192488 1.3529200 -1.3625585 -0.54561109
#>  [7,] -1.9619138 -1.6907278 0.8024144 0.7807298 -1.2808134  0.45032978
#>  [8,] -1.3538001 -0.9795329 0.7635753 0.8803881 -1.6185323 -0.51563951
#>  [9,] -1.4951981 -1.3452694 0.7483078 1.1681757 -0.8630678 -0.40056571
#> [10,] -0.7001570 -1.3105793 0.4265861 0.6063143 -1.2915884  0.36626078
#>           [,31]       [,32]     [,33]       [,34]
#>  [1,] 0.6769140  0.45327664 -2.274928 -0.05510173
#>  [2,] 1.6138060  0.21161020 -2.535652 -0.10920675
#>  [3,] 1.0329460  0.71972460 -3.001135  0.63506287
#>  [4,] 0.7265553 -0.08868924 -2.422893 -0.10321454
#>  [5,] 0.7890141  0.67393074 -2.184437 -0.12389639
#>  [6,] 1.8548737  0.49997956 -2.491934  0.37293002
#>  [7,] 0.9496569  0.97199028 -2.784751 -0.20664749
#>  [8,] 0.5995813 -0.25677041 -1.612093 -0.29992384
#>  [9,] 1.4783761  0.68995579 -1.948575 -0.40795007
#> [10,] 1.1216384  0.21703659 -2.242087 -0.33931344

kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 145.8283 148.4727           NA     NA        NA                NA
#> Nb_Comp_1 118.1398 123.4285           NA 0.0975        NA                NA
#> Nb_Comp_2 109.9553 117.8885           NA 0.0975        NA                NA
#> Nb_Comp_3 105.1591 115.7366           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y    RSS_Y      R2_Y
#> Nb_Comp_0      104.00000 25.91346        NA
#> Nb_Comp_1      100.53823 19.32272 0.2543365
#> Nb_Comp_2       99.17955 17.33735 0.3309519
#> Nb_Comp_3      123.37836 15.58198 0.3986915
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
rm(list=c("bbb","bbb2"))



data(pine)
bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,
modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb <- cv.plsRglm(round(x11)~.,data=pine,nt=10,
modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>           [,1]         [,2]        [,3]       [,4]      [,5]      [,6]
#>  [1,] 13.43440 -0.005020346 -0.07102968 0.21380483 -1.576453 0.3012073
#>  [2,] 13.71955 -0.004604089 -0.06125666 0.26993320 -1.832569 0.3159689
#>  [3,] 10.59279 -0.003987051 -0.05439997 0.16337416 -1.217370 0.2316912
#>  [4,] 11.89610 -0.005057295 -0.07836379 0.15835572 -1.495708 0.3082888
#>  [5,] 14.93021 -0.006382132 -0.06393526 0.28622870 -1.674858 0.3683381
#>  [6,] 14.25039 -0.007809612 -0.07682626 0.08385396 -2.353049 0.4831968
#>  [7,] 12.14547 -0.003915827 -0.09098797 0.19647456 -2.165601 0.4008132
#>  [8,] 11.34905 -0.004015198 -0.06397980 0.15729686 -1.559873 0.2735091
#>  [9,] 13.68196 -0.006068322 -0.06995102 0.20394086 -1.146496 0.2648040
#> [10,] 11.09734 -0.004425391 -0.07139576 0.14445299 -1.553145 0.3067691
#>            [,7]        [,8]        [,9]      [,10]       [,11]
#>  [1,] -2.436453 -0.09249864  0.27182416 -1.2796871 -0.37870030
#>  [2,] -3.563249  0.67852967  0.31162528 -0.9422022 -0.98576020
#>  [3,] -1.535269 -0.18266390  0.11196531 -1.4066954  0.27066190
#>  [4,] -1.725937  0.04360096  0.20761316 -1.0202513 -0.13200144
#>  [5,] -3.233393  0.31120118  0.18421458 -1.3188269 -0.34788498
#>  [6,] -1.251317  0.70607959  0.25570572 -0.5042740 -0.47483515
#>  [7,] -2.102087  0.11738924  0.37224207 -1.6058639 -0.09056626
#>  [8,] -1.463213 -0.26582812  0.31781719 -1.6796060 -0.08296238
#>  [9,] -1.868919  0.08659158 -0.01138339 -0.9786412 -0.37504059
#> [10,] -1.336106  0.01461444  0.17950231 -1.1905006 -0.04232857
boxplot(kfolds2coeff(bbb)[,1])


kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
summary(bbb)
#> ____************************************************____
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0  76.61170 78.10821           NA     NA        NA                NA
#> Nb_Comp_1  65.70029 68.69331           NA 0.0975        NA                NA
#> Nb_Comp_2  62.49440 66.98392           NA 0.0975        NA                NA
#> Nb_Comp_3  62.47987 68.46590           NA 0.0975        NA                NA
#> Nb_Comp_4  64.21704 71.69958           NA 0.0975        NA                NA
#> Nb_Comp_5  65.81654 74.79559           NA 0.0975        NA                NA
#> Nb_Comp_6  66.48888 76.96443           NA 0.0975        NA                NA
#> Nb_Comp_7  68.40234 80.37440           NA 0.0975        NA                NA
#> Nb_Comp_8  70.39399 83.86256           NA 0.0975        NA                NA
#> Nb_Comp_9  72.37642 87.34149           NA 0.0975        NA                NA
#> Nb_Comp_10 74.37612 90.83770           NA 0.0975        NA                NA
#>            Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0        33.75000 24.545455        NA
#> Nb_Comp_1        23.85891 12.599337 0.4866937
#> Nb_Comp_2        17.29992  9.056074 0.6310488
#> Nb_Comp_3        15.50937  8.232069 0.6646194
#> Nb_Comp_4        15.23934  8.125808 0.6689485
#> Nb_Comp_5        15.26275  7.862134 0.6796909
#> Nb_Comp_6        17.74629  6.203270 0.7472742
#> Nb_Comp_7        18.04460  5.879880 0.7604493
#> Nb_Comp_8        18.17881  5.827065 0.7626011
#> Nb_Comp_9        18.34925  5.837300 0.7621841
#> Nb_Comp_10       18.39332  5.832437 0.7623822
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(x11~.,data=pine,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                 AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0  82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1  63.61896  0.38248575  0.0975  0.38248575 12.844390 11.074659
#> Nb_Comp_2  58.47638  0.34836456  0.0975 -0.05525570 11.686597  8.919303
#> Nb_Comp_3  56.55421  0.23688359  0.0975 -0.17107874 10.445206  7.919786
#> Nb_Comp_4  54.35053  0.06999681  0.0975 -0.21869112  9.651773  6.972542
#> Nb_Comp_5  55.99834 -0.07691053  0.0975 -0.15796434  8.073955  6.898523
#> Nb_Comp_6  57.69592 -0.19968885  0.0975 -0.11400977  7.685022  6.835594
#> Nb_Comp_7  59.37953 -0.27722139  0.0975 -0.06462721  7.277359  6.770369
#> Nb_Comp_8  61.21213 -0.30602578  0.0975 -0.02255238  6.923057  6.736112
#> Nb_Comp_9  63.18426 -0.39920228  0.0975 -0.07134354  7.216690  6.730426
#> Nb_Comp_10 65.15982 -0.43743644  0.0975 -0.02732569  6.914340  6.725443
#>                 R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0         NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1  0.4675684 0.4675684   17.03781     19.76046  0.38248575 0.0975
#> Nb_Comp_2  0.5711905 0.5711905   13.72190     17.97925 -0.05525570 0.0975
#> Nb_Comp_3  0.6192438 0.6192438   12.18420     16.06943 -0.17107874 0.0975
#> Nb_Comp_4  0.6647841 0.6647841   10.72691     14.84877 -0.21869112 0.0975
#> Nb_Comp_5  0.6683426 0.6683426   10.61304     12.42138 -0.15796434 0.0975
#> Nb_Comp_6  0.6713681 0.6713681   10.51622     11.82303 -0.11400977 0.0975
#> Nb_Comp_7  0.6745039 0.6745039   10.41588     11.19586 -0.06462721 0.0975
#> Nb_Comp_8  0.6761508 0.6761508   10.36317     10.65078 -0.02255238 0.0975
#> Nb_Comp_9  0.6764242 0.6764242   10.35443     11.10252 -0.07134354 0.0975
#> Nb_Comp_10 0.6766638 0.6766638   10.34676     10.63737 -0.02732569 0.0975
#>            Q2cum_residY  AIC.std   DoF.dof sigmahat.dof   AIC.dof   BIC.dof
#> Nb_Comp_0            NA 96.63448  1.000000    0.8062287 0.6697018 0.6991787
#> Nb_Comp_1    0.38248575 77.83455  3.176360    0.5994089 0.4047616 0.4565153
#> Nb_Comp_2    0.34836456 72.69198  7.133559    0.5761829 0.4138120 0.5212090
#> Nb_Comp_3    0.23688359 70.76981  8.778329    0.5603634 0.4070516 0.5320535
#> Nb_Comp_4    0.06999681 68.56612  8.427874    0.5221703 0.3505594 0.4547689
#> Nb_Comp_5   -0.07691053 70.21393  9.308247    0.5285695 0.3666578 0.4845912
#> Nb_Comp_6   -0.19968885 71.91152  9.291931    0.5259794 0.3629363 0.4795121
#> Nb_Comp_7   -0.27722139 73.59512  9.756305    0.5284535 0.3702885 0.4938445
#> Nb_Comp_8   -0.30602578 75.42772 10.363948    0.5338475 0.3831339 0.5170783
#> Nb_Comp_9   -0.39920228 77.39986 10.732146    0.5378276 0.3920957 0.5328746
#> Nb_Comp_10  -0.43743644 79.37542 11.000000    0.5407500 0.3987417 0.5446065
#>             GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0  -3.605128         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1  -9.875081         2      0.5977015 0.3788984 0.4112998 -11.451340
#> Nb_Comp_2  -6.985517         3      0.5452615 0.3243383 0.3647862 -12.822703
#> Nb_Comp_3  -6.260610         4      0.5225859 0.3061986 0.3557368 -12.756838
#> Nb_Comp_4  -8.152986         5      0.4990184 0.2867496 0.3432131 -12.811575
#> Nb_Comp_5  -7.111583         6      0.5054709 0.3019556 0.3714754 -11.329638
#> Nb_Comp_6  -7.233043         7      0.5127450 0.3186757 0.4021333  -9.918688
#> Nb_Comp_7  -6.742195         8      0.5203986 0.3364668 0.4347156  -8.592770
#> Nb_Comp_8  -6.038372         9      0.5297842 0.3572181 0.4717708  -7.287834
#> Nb_Comp_9  -5.600237        10      0.5409503 0.3813021 0.5140048  -6.008747
#> Nb_Comp_10 -5.288422        11      0.5529032 0.4076026 0.5600977  -4.799453

data(pineNAX21)
bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb2 <- cv.plsRglm(round(x11)~.,data=pineNAX21,nt=10,
modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
#>           [,1]         [,2]        [,3]        [,4]       [,5]      [,6]
#>  [1,] 11.93974 -0.005451153 -0.01513575  0.07293697 -1.9802573 0.3297352
#>  [2,] 13.95860 -0.004698416 -0.06250622  0.24567067 -1.1614547 0.1829770
#>  [3,] 14.33562 -0.005968918 -0.05306840  0.26155441 -1.4106086 0.3524191
#>  [4,] -0.95380 -0.003889616 -0.02094643 -0.40660268 -0.2379363 0.4700880
#>  [5,] 16.91036 -0.003613692 -0.17422546  0.38123625 -3.7546680 0.7095015
#>  [6,] 16.50836 -0.006039170 -0.07570081  0.23222456 -1.6627985 0.2815828
#>  [7,] 13.32782 -0.005093032 -0.07009140  0.20915251 -1.4816230 0.3128028
#>  [8,] 13.54216 -0.004381571 -0.06509680  0.17877374 -1.3388351 0.2590691
#>  [9,] 12.33722 -0.004615619 -0.05868211  0.17233204 -1.5654872 0.3141138
#> [10,] 12.75584 -0.004969525 -0.07158048  0.19656792 -1.6231507 0.3139152
#>            [,7]        [,8]        [,9]      [,10]       [,11]
#>  [1,] -1.471488 -0.31569009  0.40041036 -0.7978037  0.24257568
#>  [2,] -2.879842 -0.15397870  0.34803702 -1.6800830 -0.18890315
#>  [3,] -2.971420  0.26736657  0.06220152 -1.1311851 -0.22968578
#>  [4,]  6.027726 -0.51308092 -0.94737365  0.5729701  0.74382441
#>  [5,] -6.543311  2.39186343  0.52986952  0.5726013 -1.78426859
#>  [6,] -2.264299 -0.19053645  0.34584736 -1.9090284 -0.17393189
#>  [7,] -2.348482  0.25693332  0.06439831 -0.6502272 -0.31032315
#>  [8,] -1.817733 -0.52079705  0.33565284 -1.9505093 -0.24126970
#>  [9,] -1.723525  0.13694904  0.12240060 -1.1251920 -0.08248523
#> [10,] -2.125138  0.03076682  0.23504469 -1.2578376 -0.21416307
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 76.61170 78.10821           NA     NA        NA                NA
#> Nb_Comp_1 65.74449 68.73751           NA 0.0975        NA                NA
#> Nb_Comp_2 62.35674 66.84626           NA 0.0975        NA                NA
#> Nb_Comp_3 62.39804 68.38407           NA 0.0975        NA                NA
#> Nb_Comp_4 64.08113 71.56366           NA 0.0975        NA                NA
#> Nb_Comp_5 65.63784 74.61689           NA 0.0975        NA                NA
#> Nb_Comp_6 67.18468 77.66024           NA 0.0975        NA                NA
#> Nb_Comp_7 68.61004 80.58210           NA 0.0975        NA                NA
#> Nb_Comp_8 70.54487 84.01344           NA 0.0975        NA                NA
#> Nb_Comp_9 72.37296 87.33803           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0       33.75000 24.545455        NA
#> Nb_Comp_1       23.89105 12.654950 0.4844280
#> Nb_Comp_2       17.31172  8.871122 0.6385839
#> Nb_Comp_3       15.51670  8.203709 0.6657748
#> Nb_Comp_4       15.31216  7.959332 0.6757309
#> Nb_Comp_5       15.51159  7.724832 0.6852846
#> Nb_Comp_6       16.30549  6.814620 0.7223673
#> Nb_Comp_7       17.52007  6.284737 0.7439552
#> Nb_Comp_8       17.75766  6.160827 0.7490034
#> Nb_Comp_9       18.30206  5.831059 0.7624383
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(x11~.,data=pineNAX21,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> ____There are some NAs in X but not in Y____
#> ____TypeVC____ standard ____
#> ____TypeVC____ standard ____unknown____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE] < 10^{-12}
#> Warning only 9 components could thus be extracted
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#>                AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0 82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1 63.69250  0.35639805  0.0975  0.35639805 13.387018 11.099368
#> Nb_Comp_2 58.35228  0.28395028  0.0975 -0.11256611 12.348781  8.885823
#> Nb_Comp_3 56.36553  0.07664889  0.0975 -0.28950699 11.458331  7.874634
#> Nb_Comp_4 54.02416 -0.70355579  0.0975 -0.84497074 14.528469  6.903925
#> Nb_Comp_5 55.80450 -0.94905654  0.0975 -0.14411078  7.898855  6.858120
#> Nb_Comp_6 57.45753 -1.27568315  0.0975 -0.16758190  8.007417  6.786392
#> Nb_Comp_7 58.73951 -1.63309014  0.0975 -0.15705481  7.852227  6.640327
#> Nb_Comp_8 60.61227 -1.67907859  0.0975 -0.01746558  6.756304  6.614773
#> Nb_Comp_9 62.25948 -2.15165796  0.0975 -0.17639623  7.781594  6.544432
#>                R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0        NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1 0.4663804 0.4663804   17.07583     20.59526  0.35639805 0.0975
#> Nb_Comp_2 0.5728001 0.5728001   13.67040     18.99799 -0.11256611 0.0975
#> Nb_Comp_3 0.6214146 0.6214146   12.11473     17.62807 -0.28950699 0.0975
#> Nb_Comp_4 0.6680830 0.6680830   10.62135     22.35133 -0.84497074 0.0975
#> Nb_Comp_5 0.6702851 0.6702851   10.55088     12.15200 -0.14411078 0.0975
#> Nb_Comp_6 0.6737336 0.6737336   10.44053     12.31901 -0.16758190 0.0975
#> Nb_Comp_7 0.6807558 0.6807558   10.21581     12.08026 -0.15705481 0.0975
#> Nb_Comp_8 0.6819844 0.6819844   10.17650     10.39424 -0.01746558 0.0975
#> Nb_Comp_9 0.6853661 0.6853661   10.06828     11.97160 -0.17639623 0.0975
#>           Q2cum_residY  AIC.std DoF.dof sigmahat.dof AIC.dof BIC.dof GMDL.dof
#> Nb_Comp_0           NA 96.63448      NA           NA      NA      NA       NA
#> Nb_Comp_1   0.35639805 77.90810      NA           NA      NA      NA       NA
#> Nb_Comp_2   0.28395028 72.56787      NA           NA      NA      NA       NA
#> Nb_Comp_3   0.07664889 70.58113      NA           NA      NA      NA       NA
#> Nb_Comp_4  -0.70355579 68.23976      NA           NA      NA      NA       NA
#> Nb_Comp_5  -0.94905654 70.02009      NA           NA      NA      NA       NA
#> Nb_Comp_6  -1.27568315 71.67313      NA           NA      NA      NA       NA
#> Nb_Comp_7  -1.63309014 72.95511      NA           NA      NA      NA       NA
#> Nb_Comp_8  -1.67907859 74.82787      NA           NA      NA      NA       NA
#> Nb_Comp_9  -2.15165796 76.47507      NA           NA      NA      NA       NA
#>           DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1         2      0.5983679 0.3797438 0.4122175 -11.413749
#> Nb_Comp_2         3      0.5442372 0.3231208 0.3634169 -12.847656
#> Nb_Comp_3         4      0.5210941 0.3044529 0.3537087 -12.776843
#> Nb_Comp_4         5      0.4965569 0.2839276 0.3398355 -12.891035
#> Nb_Comp_5         6      0.5039885 0.3001871 0.3692997 -11.349498
#> Nb_Comp_6         7      0.5108963 0.3163819 0.3992388  -9.922119
#> Nb_Comp_7         8      0.5153766 0.3300041 0.4263658  -8.696873
#> Nb_Comp_8         9      0.5249910 0.3507834 0.4632727  -7.337679
#> Nb_Comp_9        10      0.5334234 0.3707649 0.4998004  -6.033403
rm(list=c("bbb","bbb2"))



data(pine)
bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-family",
family=Gamma,K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb <- cv.plsRglm(x11~.,data=pine,nt=10,modele="pls-glm-Gamma",
K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
#>            [,1]        [,2]        [,3]       [,4]     [,5]       [,6]
#>  [1,] -13.95857 0.005492044 0.077891499 -0.2695253 2.767389 -0.5047581
#>  [2,] -14.31515 0.007928466 0.015030334 -0.2008144 1.669468 -0.4086851
#>  [3,] -10.91445 0.004996889 0.026411511 -0.1360011 1.737344 -0.3257053
#>  [4,] -11.45593 0.005668026 0.034349894 -0.1636443 1.509233 -0.3207658
#>  [5,] -12.63340 0.006648547 0.059439953 -0.1229428 2.401614 -0.4610994
#>  [6,] -13.39863 0.006341036 0.036360527 -0.2297059 1.547118 -0.3087343
#>  [7,] -11.28863 0.004808237 0.044881971 -0.2194449 1.326409 -0.2892386
#>  [8,] -13.42320 0.005774734 0.018779268 -0.1310495 1.013063 -0.2342447
#>  [9,] -10.52738 0.008452585 0.007726981  0.1829999 1.975441 -0.4190091
#> [10,] -10.66290 0.005346505 0.045656541 -0.1292047 1.829342 -0.3524160
#>            [,7]        [,8]        [,9]     [,10]     [,11]
#>  [1,]  3.892706 -0.97671926 -0.54252584 1.0714006 1.1478947
#>  [2,]  1.338291  0.14027394 -0.06187781 1.9735920 0.3743748
#>  [3,]  1.730714 -0.13559519 -0.15012613 0.8618811 0.5970222
#>  [4,]  2.212397 -0.36603986 -0.10691877 0.5502711 0.8627822
#>  [5,]  1.593489 -0.50499819 -0.35236772 0.8618332 0.6842768
#>  [6,]  2.876549 -0.34983124 -0.24479651 1.4340965 0.6853397
#>  [7,]  2.604138 -0.06249697 -0.07013608 0.9064415 0.5230287
#>  [8,]  1.601278  1.38080889 -0.06193123 1.2014655 0.5202962
#>  [9,] -1.982126 -0.34224934 -0.11175185 0.5840213 0.6573964
#> [10,]  1.911670 -0.58048487 -0.29458841 0.7000467 0.8239750
boxplot(kfolds2coeff(bbb)[,1])


kfolds2Chisqind(bbb)
#> [[1]]
#> [[1]][[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb)
#> [[1]]
#>  [1] NA NA NA NA NA NA NA NA NA NA
#> 
summary(bbb)
#> ____************************************************____
#> 
#> Family: Gamma 
#> Link function: inverse 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                 AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0  56.60919 59.60220           NA     NA        NA                NA
#> Nb_Comp_1  39.01090 43.50042           NA 0.0975        NA                NA
#> Nb_Comp_2  37.30801 43.29404           NA 0.0975        NA                NA
#> Nb_Comp_3  36.87524 44.35777           NA 0.0975        NA                NA
#> Nb_Comp_4  36.55795 45.53700           NA 0.0975        NA                NA
#> Nb_Comp_5  37.13611 47.61167           NA 0.0975        NA                NA
#> Nb_Comp_6  38.27656 50.24862           NA 0.0975        NA                NA
#> Nb_Comp_7  39.39377 52.86234           NA 0.0975        NA                NA
#> Nb_Comp_8  40.96122 55.92630           NA 0.0975        NA                NA
#> Nb_Comp_9  42.90816 59.36974           NA 0.0975        NA                NA
#> Nb_Comp_10 44.90815 62.86625           NA 0.0975        NA                NA
#>            Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0        31.60805 20.800152        NA
#> Nb_Comp_1        17.31431 11.804594 0.4324756
#> Nb_Comp_2        17.01037  6.357437 0.6943562
#> Nb_Comp_3        15.83422  5.699662 0.7259798
#> Nb_Comp_4        13.52676  7.679741 0.6307844
#> Nb_Comp_5        13.60962  6.099077 0.7067773
#> Nb_Comp_6        13.91155  5.205052 0.7497590
#> Nb_Comp_7        14.94390  4.650377 0.7764258
#> Nb_Comp_8        15.25537  4.321314 0.7922461
#> Nb_Comp_9        15.15577  4.307757 0.7928978
#> Nb_Comp_10       15.15490  4.307391 0.7929154
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(x11~.,data=pine,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Component____ 10 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                 AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0  82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1  63.61896  0.38248575  0.0975  0.38248575 12.844390 11.074659
#> Nb_Comp_2  58.47638  0.34836456  0.0975 -0.05525570 11.686597  8.919303
#> Nb_Comp_3  56.55421  0.23688359  0.0975 -0.17107874 10.445206  7.919786
#> Nb_Comp_4  54.35053  0.06999681  0.0975 -0.21869112  9.651773  6.972542
#> Nb_Comp_5  55.99834 -0.07691053  0.0975 -0.15796434  8.073955  6.898523
#> Nb_Comp_6  57.69592 -0.19968885  0.0975 -0.11400977  7.685022  6.835594
#> Nb_Comp_7  59.37953 -0.27722139  0.0975 -0.06462721  7.277359  6.770369
#> Nb_Comp_8  61.21213 -0.30602578  0.0975 -0.02255238  6.923057  6.736112
#> Nb_Comp_9  63.18426 -0.39920228  0.0975 -0.07134354  7.216690  6.730426
#> Nb_Comp_10 65.15982 -0.43743644  0.0975 -0.02732569  6.914340  6.725443
#>                 R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0         NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1  0.4675684 0.4675684   17.03781     19.76046  0.38248575 0.0975
#> Nb_Comp_2  0.5711905 0.5711905   13.72190     17.97925 -0.05525570 0.0975
#> Nb_Comp_3  0.6192438 0.6192438   12.18420     16.06943 -0.17107874 0.0975
#> Nb_Comp_4  0.6647841 0.6647841   10.72691     14.84877 -0.21869112 0.0975
#> Nb_Comp_5  0.6683426 0.6683426   10.61304     12.42138 -0.15796434 0.0975
#> Nb_Comp_6  0.6713681 0.6713681   10.51622     11.82303 -0.11400977 0.0975
#> Nb_Comp_7  0.6745039 0.6745039   10.41588     11.19586 -0.06462721 0.0975
#> Nb_Comp_8  0.6761508 0.6761508   10.36317     10.65078 -0.02255238 0.0975
#> Nb_Comp_9  0.6764242 0.6764242   10.35443     11.10252 -0.07134354 0.0975
#> Nb_Comp_10 0.6766638 0.6766638   10.34676     10.63737 -0.02732569 0.0975
#>            Q2cum_residY  AIC.std   DoF.dof sigmahat.dof   AIC.dof   BIC.dof
#> Nb_Comp_0            NA 96.63448  1.000000    0.8062287 0.6697018 0.6991787
#> Nb_Comp_1    0.38248575 77.83455  3.176360    0.5994089 0.4047616 0.4565153
#> Nb_Comp_2    0.34836456 72.69198  7.133559    0.5761829 0.4138120 0.5212090
#> Nb_Comp_3    0.23688359 70.76981  8.778329    0.5603634 0.4070516 0.5320535
#> Nb_Comp_4    0.06999681 68.56612  8.427874    0.5221703 0.3505594 0.4547689
#> Nb_Comp_5   -0.07691053 70.21393  9.308247    0.5285695 0.3666578 0.4845912
#> Nb_Comp_6   -0.19968885 71.91152  9.291931    0.5259794 0.3629363 0.4795121
#> Nb_Comp_7   -0.27722139 73.59512  9.756305    0.5284535 0.3702885 0.4938445
#> Nb_Comp_8   -0.30602578 75.42772 10.363948    0.5338475 0.3831339 0.5170783
#> Nb_Comp_9   -0.39920228 77.39986 10.732146    0.5378276 0.3920957 0.5328746
#> Nb_Comp_10  -0.43743644 79.37542 11.000000    0.5407500 0.3987417 0.5446065
#>             GMDL.dof DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0  -3.605128         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1  -9.875081         2      0.5977015 0.3788984 0.4112998 -11.451340
#> Nb_Comp_2  -6.985517         3      0.5452615 0.3243383 0.3647862 -12.822703
#> Nb_Comp_3  -6.260610         4      0.5225859 0.3061986 0.3557368 -12.756838
#> Nb_Comp_4  -8.152986         5      0.4990184 0.2867496 0.3432131 -12.811575
#> Nb_Comp_5  -7.111583         6      0.5054709 0.3019556 0.3714754 -11.329638
#> Nb_Comp_6  -7.233043         7      0.5127450 0.3186757 0.4021333  -9.918688
#> Nb_Comp_7  -6.742195         8      0.5203986 0.3364668 0.4347156  -8.592770
#> Nb_Comp_8  -6.038372         9      0.5297842 0.3572181 0.4717708  -7.287834
#> Nb_Comp_9  -5.600237        10      0.5409503 0.3813021 0.5140048  -6.008747
#> Nb_Comp_10 -5.288422        11      0.5529032 0.4076026 0.5600977  -4.799453

data(pineNAX21)
bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
modele="pls-glm-family",family=Gamma(),K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
bbb2 <- cv.plsRglm(x11~.,data=pineNAX21,nt=10,
modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE,verbose=FALSE)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced

#For Jackknife computations
kfolds2coeff(bbb2)
#>            [,1]        [,2]         [,3]       [,4]     [,5]       [,6]
#>  [1,] -12.93267 0.005355558  0.026490638 -0.1361703 1.181354 -0.2711536
#>  [2,] -10.68628 0.005070477  0.045378213 -0.1657417 1.637885 -0.3348348
#>  [3,] -10.23650 0.006104666 -0.005200735  0.1072118 1.399918 -0.3057119
#>  [4,] -15.82110 0.006599564  0.034957995 -0.1981941 1.599511 -0.2941216
#>  [5,] -16.77429 0.006691047  0.025212062 -0.2772409 1.371117 -0.3151450
#>  [6,] -15.11770 0.003988239  0.086199586 -0.2634136 2.770033 -0.5434017
#>  [7,] -17.94653 0.007155095  0.080455708 -0.1815021 1.997504 -0.3825865
#>  [8,] -19.12501 0.007042758  0.030238172 -0.2956674 2.617566 -0.4336634
#>  [9,] -12.60601 0.005492070  0.034701993 -0.1052023 1.717352 -0.3699873
#> [10,] -14.42037 0.005471761  0.056064380 -0.1574416 1.989522 -0.3604814
#>            [,7]        [,8]        [,9]     [,10]      [,11]
#>  [1,]  1.575278  0.19409739 -0.02018638 0.8844664  0.4685616
#>  [2,]  2.173151 -0.36687481 -0.19113342 0.7898988  0.6057814
#>  [3,] -1.054361 -0.01352571 -0.01787996 0.5276048  0.4054358
#>  [4,]  2.676395 -0.30912503 -0.36341991 1.5384597  0.6535078
#>  [5,]  3.127645 -0.02657179 -0.14277232 1.7599397  0.6469898
#>  [6,]  4.342331 -0.79110951 -0.35631988 0.1604043  1.3067018
#>  [7,]  2.400094  0.18002725 -0.34981254 1.2262459  0.4880984
#>  [8,]  4.118963 -1.18427933 -0.51865548 1.1052226  2.0273600
#>  [9,]  1.411636 -0.28570562 -0.07420447 0.5307584  0.5693091
#> [10,]  1.808181  0.44381538 -0.27619117 1.2245908 -0.1806527
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[6]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[7]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[8]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[9]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> [[1]][[10]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Family: Gamma 
#> Link function: inverse 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 56.60919 59.60220           NA     NA        NA                NA
#> Nb_Comp_1 39.08940 43.57892           NA 0.0975        NA                NA
#> Nb_Comp_2 37.36154 43.34757           NA 0.0975        NA                NA
#> Nb_Comp_3 36.81173 44.29427           NA 0.0975        NA                NA
#> Nb_Comp_4 36.53654 45.51559           NA 0.0975        NA                NA
#> Nb_Comp_5 37.24312 47.71867           NA 0.0975        NA                NA
#> Nb_Comp_6 38.18649 50.15855           NA 0.0975        NA                NA
#> Nb_Comp_7 39.35575 52.82432           NA 0.0975        NA                NA
#> Nb_Comp_8 40.86209 55.82716           NA 0.0975        NA                NA
#> Nb_Comp_9 42.80511 59.26669           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y     RSS_Y      R2_Y
#> Nb_Comp_0       31.60805 20.800152        NA
#> Nb_Comp_1       17.30890 12.031518 0.4215659
#> Nb_Comp_2       17.10360  6.183372 0.7027247
#> Nb_Comp_3       15.78579  5.756462 0.7232490
#> Nb_Comp_4       13.49013  7.630460 0.6331536
#> Nb_Comp_5       13.56918  6.303455 0.6969515
#> Nb_Comp_6       14.02295  5.274716 0.7464097
#> Nb_Comp_7       15.05896  4.867806 0.7659726
#> Nb_Comp_8       15.28052  4.317488 0.7924300
#> Nb_Comp_9       15.19429  4.298593 0.7933384
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(x11~.,data=pineNAX21,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> ____There are some NAs in X but not in Y____
#> ____TypeVC____ standard ____
#> ____TypeVC____ standard ____unknown____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Component____ 8 ____
#> ____Component____ 9 ____
#> Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[1,],,drop=FALSE] < 10^{-12}
#> Warning only 9 components could thus be extracted
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#>                AIC     Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0 82.41888          NA      NA          NA        NA 20.800152
#> Nb_Comp_1 63.69250  0.35639805  0.0975  0.35639805 13.387018 11.099368
#> Nb_Comp_2 58.35228  0.28395028  0.0975 -0.11256611 12.348781  8.885823
#> Nb_Comp_3 56.36553  0.07664889  0.0975 -0.28950699 11.458331  7.874634
#> Nb_Comp_4 54.02416 -0.70355579  0.0975 -0.84497074 14.528469  6.903925
#> Nb_Comp_5 55.80450 -0.94905654  0.0975 -0.14411078  7.898855  6.858120
#> Nb_Comp_6 57.45753 -1.27568315  0.0975 -0.16758190  8.007417  6.786392
#> Nb_Comp_7 58.73951 -1.63309014  0.0975 -0.15705481  7.852227  6.640327
#> Nb_Comp_8 60.61227 -1.67907859  0.0975 -0.01746558  6.756304  6.614773
#> Nb_Comp_9 62.25948 -2.15165796  0.0975 -0.17639623  7.781594  6.544432
#>                R2_Y R2_residY RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0        NA        NA   32.00000           NA          NA     NA
#> Nb_Comp_1 0.4663804 0.4663804   17.07583     20.59526  0.35639805 0.0975
#> Nb_Comp_2 0.5728001 0.5728001   13.67040     18.99799 -0.11256611 0.0975
#> Nb_Comp_3 0.6214146 0.6214146   12.11473     17.62807 -0.28950699 0.0975
#> Nb_Comp_4 0.6680830 0.6680830   10.62135     22.35133 -0.84497074 0.0975
#> Nb_Comp_5 0.6702851 0.6702851   10.55088     12.15200 -0.14411078 0.0975
#> Nb_Comp_6 0.6737336 0.6737336   10.44053     12.31901 -0.16758190 0.0975
#> Nb_Comp_7 0.6807558 0.6807558   10.21581     12.08026 -0.15705481 0.0975
#> Nb_Comp_8 0.6819844 0.6819844   10.17650     10.39424 -0.01746558 0.0975
#> Nb_Comp_9 0.6853661 0.6853661   10.06828     11.97160 -0.17639623 0.0975
#>           Q2cum_residY  AIC.std DoF.dof sigmahat.dof AIC.dof BIC.dof GMDL.dof
#> Nb_Comp_0           NA 96.63448      NA           NA      NA      NA       NA
#> Nb_Comp_1   0.35639805 77.90810      NA           NA      NA      NA       NA
#> Nb_Comp_2   0.28395028 72.56787      NA           NA      NA      NA       NA
#> Nb_Comp_3   0.07664889 70.58113      NA           NA      NA      NA       NA
#> Nb_Comp_4  -0.70355579 68.23976      NA           NA      NA      NA       NA
#> Nb_Comp_5  -0.94905654 70.02009      NA           NA      NA      NA       NA
#> Nb_Comp_6  -1.27568315 71.67313      NA           NA      NA      NA       NA
#> Nb_Comp_7  -1.63309014 72.95511      NA           NA      NA      NA       NA
#> Nb_Comp_8  -1.67907859 74.82787      NA           NA      NA      NA       NA
#> Nb_Comp_9  -2.15165796 76.47507      NA           NA      NA      NA       NA
#>           DoF.naive sigmahat.naive AIC.naive BIC.naive GMDL.naive
#> Nb_Comp_0         1      0.8062287 0.6697018 0.6991787  -3.605128
#> Nb_Comp_1         2      0.5983679 0.3797438 0.4122175 -11.413749
#> Nb_Comp_2         3      0.5442372 0.3231208 0.3634169 -12.847656
#> Nb_Comp_3         4      0.5210941 0.3044529 0.3537087 -12.776843
#> Nb_Comp_4         5      0.4965569 0.2839276 0.3398355 -12.891035
#> Nb_Comp_5         6      0.5039885 0.3001871 0.3692997 -11.349498
#> Nb_Comp_6         7      0.5108963 0.3163819 0.3992388  -9.922119
#> Nb_Comp_7         8      0.5153766 0.3300041 0.4263658  -8.696873
#> Nb_Comp_8         9      0.5249910 0.3507834 0.4632727  -7.337679
#> Nb_Comp_9        10      0.5334234 0.3707649 0.4998004  -6.033403
rm(list=c("bbb","bbb2"))



data(Cornell)
summary(cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,modele="pls",verbose=FALSE))
#> ____************************************************____
#> 
#> Model: pls 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y       Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA         NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.8471863  0.0975  0.8471863 71.48572  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8074866  0.0975 -0.2597922 45.02810  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.6262018  0.0975 -0.9416733 21.48773   4.418081 0.9905556
#> Nb_Comp_4 33.76477 -0.3882776  0.0975 -2.7139759 16.40865   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -4.4710117  0.0975 -2.9408628 16.98211   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975         NA       NA   3.496074 0.9925265
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000002    0.7633343  0.9711321  1.1359501 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"

cv.plsRglm(Y~.,data=Cornell,nt=3,
modele="pls-glm-inverse.gaussian",K=12,verbose=FALSE)
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 12
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=inverse.gaussian,K=12,verbose=FALSE)
#> Number of repeated crossvalidations:
#> [1] 1
#> Number of folds for each crossvalidation:
#> [1] 12
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 4    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 2    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[1]][[4]]
#>    [,1] [,2] [,3]
#> 7    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[2]][[2]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[2]][[3]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 5    NA   NA   NA
#> 
#> [[2]][[4]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 10   NA   NA   NA
#> 
#> 
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",
family=inverse.gaussian(),K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[1]][[2]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 3    NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 7   NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 2    NA   NA   NA
#> 
#> [[2]][[3]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 4    NA   NA   NA
#> 
#> [[2]][[4]]
#>    [,1] [,2] [,3]
#> 1    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 8   NA   NA   NA
#> 3   NA   NA   NA
#> 
#> [[2]][[6]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> 
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,
NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 5    NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 3    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[1]][[4]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 9    NA   NA   NA
#> 
#> [[1]][[5]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 2   NA   NA   NA
#> 1   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   [,1] [,2] [,3]
#> 3   NA   NA   NA
#> 4   NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[2]][[3]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[4]]
#>   [,1] [,2] [,3]
#> 6   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 9   NA   NA   NA
#> 5   NA   NA   NA
#> 
#> [[2]][[6]]
#>    [,1] [,2] [,3]
#> 1    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> 
cv.plsRglm(Y~.,data=Cornell,nt=3,modele="pls-glm-family",
family=inverse.gaussian(link = "1/mu^2"),K=6,NK=2,verbose=FALSE)$results_kfolds
#> [[1]]
#> [[1]][[1]]
#>    [,1] [,2] [,3]
#> 10   NA   NA   NA
#> 9    NA   NA   NA
#> 
#> [[1]][[2]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> [[1]][[3]]
#>    [,1] [,2] [,3]
#> 12   NA   NA   NA
#> 4    NA   NA   NA
#> 
#> [[1]][[4]]
#>   [,1] [,2] [,3]
#> 5   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[1]][[5]]
#>    [,1] [,2] [,3]
#> 11   NA   NA   NA
#> 3    NA   NA   NA
#> 
#> [[1]][[6]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 6   NA   NA   NA
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>    [,1] [,2] [,3]
#> 6    NA   NA   NA
#> 12   NA   NA   NA
#> 
#> [[2]][[2]]
#>    [,1] [,2] [,3]
#> 3    NA   NA   NA
#> 11   NA   NA   NA
#> 
#> [[2]][[3]]
#>    [,1] [,2] [,3]
#> 5    NA   NA   NA
#> 10   NA   NA   NA
#> 
#> [[2]][[4]]
#>   [,1] [,2] [,3]
#> 4   NA   NA   NA
#> 8   NA   NA   NA
#> 
#> [[2]][[5]]
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 9   NA   NA   NA
#> 
#> [[2]][[6]]
#>   [,1] [,2] [,3]
#> 7   NA   NA   NA
#> 2   NA   NA   NA
#> 
#> 

bbb2 <- cv.plsRglm(Y~.,data=Cornell,nt=10,
modele="pls-glm-inverse.gaussian",keepcoeffs=TRUE,verbose=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
#>              [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,] 0.0001424507  3.943443e-05 -9.014370e-06  6.703854e-05  2.664711e-05
#> [2,] 0.0013120746 -7.179339e-04 -1.175196e-03 -1.873823e-03 -1.150028e-03
#> [3,] 0.0001388570  4.047152e-04 -5.789572e-06 -5.984148e-04  1.281348e-05
#> [4,] 0.0004184020  8.662012e-04 -2.791667e-04 -2.331776e-03 -3.214869e-04
#> [5,] 0.0001440851  1.776682e-04 -1.477790e-05 -1.835599e-04  5.443790e-06
#>               [,6]          [,7]          [,8]
#> [1,] -1.110429e-05 -4.105864e-05 -1.861825e-04
#> [2,] -1.154620e-03 -1.218259e-03 -1.249177e-03
#> [3,] -2.663498e-05 -3.872940e-05 -3.725056e-05
#> [4,] -3.190328e-04 -3.460991e-04  1.883699e-04
#> [5,]  5.181211e-06 -4.806387e-05 -7.479347e-05
boxplot(kfolds2coeff(bbb2)[,1])


kfolds2Chisqind(bbb2)
#> [[1]]
#> [[1]][[1]]
#> [1] NA NA NA NA NA
#> 
#> [[1]][[2]]
#> [1] NA NA NA NA NA NA
#> 
#> [[1]][[3]]
#> [1] NA NA NA NA NA NA
#> 
#> [[1]][[4]]
#> [1] NA NA NA NA NA NA
#> 
#> [[1]][[5]]
#> [1] NA NA NA NA NA NA
#> 
#> 
kfolds2Chisq(bbb2)
#> [[1]]
#> [1] NA NA NA NA NA
#> 
summary(bbb2)
#> ____************************************************____
#> 
#> Family: inverse.gaussian 
#> Link function: 1/mu^2 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 81.67928 82.64909           NA     NA        NA                NA
#> Nb_Comp_1 49.90521 51.35993           NA 0.0975        NA                NA
#> Nb_Comp_2 31.06918 33.00881           NA 0.0975        NA                NA
#> Nb_Comp_3 28.40632 30.83085           NA 0.0975        NA                NA
#> Nb_Comp_4 27.08522 29.99466           NA 0.0975        NA                NA
#> Nb_Comp_5 28.46056 31.85490           NA 0.0975        NA                NA
#> Nb_Comp_6 29.68366 33.56292           NA 0.0975        NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0   6.729783e-04 467.796667        NA
#> Nb_Comp_1   3.957680e-05  32.478677 0.9305710
#> Nb_Comp_2   7.009452e-06   6.020269 0.9871306
#> Nb_Comp_3   4.727777e-06   3.795855 0.9918857
#> Nb_Comp_4   3.584346e-06   2.699884 0.9942285
#> Nb_Comp_5   3.408069e-06   2.598572 0.9944451
#> Nb_Comp_6   3.195402e-06   2.492371 0.9946721
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
PLS_lm_formula(Y~.,data=Cornell,10,typeVC="standard")$InfCrit
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> Warning : 1 2 3 4 5 6 7 < 10^{-12}
#> Warning only 6 components could thus be extracted
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#>                AIC   Q2cum_Y LimQ2_Y        Q2_Y   PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205        NA      NA          NA        NA 467.796667        NA
#> Nb_Comp_1 53.15173 0.8966556  0.0975  0.89665563 48.344150  35.742486 0.9235940
#> Nb_Comp_2 41.08283 0.9175426  0.0975  0.20210989 28.518576  11.066606 0.9763431
#> Nb_Comp_3 32.06411 0.9399676  0.0975  0.27195907  8.056942   4.418081 0.9905556
#> Nb_Comp_4 33.76477 0.9197009  0.0975 -0.33759604  5.909608   4.309235 0.9907882
#> Nb_Comp_5 33.34373 0.9281373  0.0975  0.10506161  3.856500   3.521924 0.9924713
#> Nb_Comp_6 35.25533 0.9232562  0.0975 -0.06792167  3.761138   3.496074 0.9925265
#>           R2_residY  RSS_residY PRESS_residY   Q2_residY  LimQ2 Q2cum_residY
#> Nb_Comp_0        NA 11.00000000           NA          NA     NA           NA
#> Nb_Comp_1 0.9235940  0.84046633   1.13678803  0.89665563 0.0975    0.8966556
#> Nb_Comp_2 0.9763431  0.26022559   0.67059977  0.20210989 0.0975    0.9175426
#> Nb_Comp_3 0.9905556  0.10388893   0.18945488  0.27195907 0.0975    0.9399676
#> Nb_Comp_4 0.9907882  0.10132947   0.13896142 -0.33759604 0.0975    0.9197009
#> Nb_Comp_5 0.9924713  0.08281624   0.09068364  0.10506161 0.0975    0.9281373
#> Nb_Comp_6 0.9925265  0.08220840   0.08844125 -0.06792167 0.0975    0.9232562
#>              AIC.std  DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof
#> Nb_Comp_0  37.010388 1.000000    6.5212706 46.0708838 47.7893514 27.59461
#> Nb_Comp_1   8.150064 2.740749    1.8665281  4.5699686  4.9558156 21.34020
#> Nb_Comp_2  -3.918831 5.085967    1.1825195  2.1075461  2.3949331 27.40202
#> Nb_Comp_3 -12.937550 5.121086    0.7488308  0.8467795  0.9628191 24.40842
#> Nb_Comp_4 -11.236891 5.103312    0.7387162  0.8232505  0.9357846 24.23105
#> Nb_Comp_5 -11.657929 6.006316    0.7096382  0.7976101  0.9198348 28.21184
#> Nb_Comp_6  -9.746328 7.000001    0.7633342  0.9711319  1.1359499 33.18347
#>           DoF.naive sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0         1      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1         2      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2         3      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3         4      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4         5      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5         6      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6         7      0.8361907  1.1070902  1.3048716   33.63927
rm(list=c("bbb","bbb2"))
#> Warning: object 'bbb' not found


data(bordeaux)
summary(cv.plsRglm(Quality~.,data=bordeaux,10,
modele="pls-glm-polr",K=7))
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> NK: 1 
#> Number of groups : 7 
#> 1 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> 2 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> 3 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> 4 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> 5 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> 6 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> 7 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> Warning :  < 10^{-12}
#> Warning only 4 components could thus be extracted
#> ****________________________________________________****
#> 
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2  Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 78.64736 81.70009           NA     NA         NA                NA
#> Nb_Comp_1 36.50286 41.08194   -0.6194641 0.0975 -0.6194641          100.9466
#>           Chi2_Pearson_Y
#> Nb_Comp_0      62.333333
#> Nb_Comp_1       9.356521
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"

data(bordeauxNA)
summary(cv.plsRglm(Quality~.,data=bordeauxNA,
10,modele="pls-glm-polr",K=10,verbose=FALSE))
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Component____ 1 ____
#> ____Predicting X with NA in X and not in Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2  Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 78.64736 81.70009           NA     NA         NA                NA
#> Nb_Comp_1 36.21263 40.79171   -0.8636295 0.0975 -0.8636295          116.1662
#>           Chi2_Pearson_Y
#> Nb_Comp_0      62.333333
#> Nb_Comp_1       9.454055
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"

summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
modele="pls-glm-polr",method="logistic",verbose=FALSE))
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 78.64736 81.70009           NA     NA        NA                NA
#> Nb_Comp_1 36.50286 41.08194    -1.382892 0.0975 -1.382892          148.5336
#>           Chi2_Pearson_Y
#> Nb_Comp_0      62.333333
#> Nb_Comp_1       9.356521
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
modele="pls-glm-polr",method="probit",verbose=FALSE))
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: probit 
#> 
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 78.64736 81.70009           NA     NA        NA                NA
#> Nb_Comp_1 36.01661 40.59569    -2.488309 0.0975 -2.488309          217.4385
#>           Chi2_Pearson_Y
#> Nb_Comp_0       62.33350
#> Nb_Comp_1        9.71675
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
modele="pls-glm-polr",method="cloglog",verbose=FALSE))
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: cloglog 
#> 
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 78.64736 81.70009           NA     NA        NA                NA
#> Nb_Comp_1 36.92722 41.50630    0.4279518 0.0975 0.4279518          35.65848
#>           Chi2_Pearson_Y
#> Nb_Comp_0       62.33474
#> Nb_Comp_1       10.32213
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
suppressWarnings(summary(cv.plsRglm(Quality~.,data=bordeaux,nt=2,K=7,
modele="pls-glm-polr",method="cauchit",verbose=FALSE)))
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: cauchit 
#> 
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2  Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 79.08163 82.13436           NA     NA         NA                NA
#> Nb_Comp_1 38.11253 42.69161   -0.3558332 0.0975 -0.3558332          84.16807
#>           Chi2_Pearson_Y
#> Nb_Comp_0      62.078483
#> Nb_Comp_1       8.592708
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
# }
```
