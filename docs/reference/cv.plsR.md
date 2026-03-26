# Partial least squares regression models with k-fold cross-validation

This function implements k-fold cross-validation on complete or
incomplete datasets for partial least squares regression models

## Usage

``` r
cv.plsR(object, ...)
# Default S3 method
cv.plsRmodel(object,dataX,nt=2,limQ2set=.0975,modele="pls", 
K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, 
scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12), weights, verbose=TRUE,...)
# S3 method for class 'formula'
cv.plsRmodel(object,data=NULL,nt=2,limQ2set=.0975,modele="pls", 
K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, 
scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL, verbose=TRUE,...)
PLS_lm_kfoldcv(dataY, dataX, nt = 2, limQ2set = 0.0975, modele = "pls", 
K = 5, NK = 1, grouplist = NULL, random = TRUE, scaleX = TRUE, 
scaleY = NULL, keepcoeffs = FALSE, keepfolds = FALSE, keepdataY = TRUE, 
keepMclassed=FALSE, tol_Xi = 10^(-12), weights, verbose=TRUE)
PLS_lm_kfoldcv_formula(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", 
K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, 
scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, 
keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL,verbose=TRUE)
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

  name of the PLS model to be fitted, only (`"pls"` available for this
  fonction.

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

  shall the number of miss classed be returned

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- contrasts:

  an optional list. See the `contrasts.arg` of `model.matrix.default`.

- verbose:

  should info messages be displayed ?

- ...:

  arguments to pass to `cv.plsRmodel.default` or to
  `cv.plsRmodel.formula`

## Details

Predicts 1 group with the `K-1` other groups. Leave one out cross
validation is thus obtained for `K==nrow(dataX)`.

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

An object of class `"cv.plsRmodel"`.  

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

  list of `NK`. Each element of the list sums up the results for a group
  division:

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

- call:

  the call of the function

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frederic Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## Author

Frederic Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Note

Work for complete and incomplete datasets.

## See also

Summary method `summary.cv.plsRmodel`.
[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.md),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.md),
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.md)
and
[`kfolds2CVinfos_lm`](https://fbertran.github.io/plsRglm/reference/kfolds2CVinfos_lm.md)
to extract and transform results from k-fold cross-validation.

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

#Leave one out CV (K=nrow(Cornell)) one time (NK=1)
bbb <- cv.plsR(object=yCornell,dataX=XCornell,nt=6,K=nrow(Cornell),NK=1)
#> NK: 1 
#> Leave One Out
#> Number of groups : 12 
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
#> ____Component____ 6 ____
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
#> 7 
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
#> 8 
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
#> 9 
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
#> 10 
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
#> 11 
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
#> 12 
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
bbb2 <- cv.plsR(Y~.,data=Cornell,nt=6,K=12,NK=1,verbose=FALSE)
(sum1<-summary(bbb2))
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
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y       Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA         NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.8809146  0.0975  0.8809146 55.70774  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8619560  0.0975 -0.1592015 41.43274  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.7471041  0.0975 -0.8319956 20.27397   4.418081 0.9905556
#> Nb_Comp_4 33.76477 -0.2159389  0.0975 -3.8080607 21.24240   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -5.9182568  0.0975 -4.6896417 24.51801   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975         NA       NA   3.496074 0.9925265
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
#> attr(,"computed_nt")
#> [1] 6
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"

#6-fold CV (K=6) two times (NK=2)
#use random=TRUE to randomly create folds for repeated CV
bbb3 <- cv.plsR(object=yCornell,dataX=XCornell,nt=6,K=6,NK=2)
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
#> Warning : 1 2 3 4 5 6 7 < 10^{-12}
#> Warning only 5 components could thus be extracted
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
#> ____Component____ 6 ____
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
#> NK: 2 
#> Number of groups : 6 
#> 1 
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
#> ____Component____ 6 ____
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
bbb4 <- cv.plsR(Y~.,data=Cornell,nt=6,K=6,NK=2,verbose=FALSE)
(sum3<-summary(bbb4))
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
#> NK: 1,  2
#> [[1]]
#>                AIC    Q2cum_Y LimQ2_Y        Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205         NA      NA          NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173  0.8856252  0.0975  0.88562522 53.50414  35.742486 0.9235940
#> Nb_Comp_2 41.08283  0.8838671  0.0975 -0.01537155 36.29190  11.066606 0.9763431
#> Nb_Comp_3 32.06411  0.7763726  0.0975 -0.92561654 21.31004   4.418081 0.9905556
#> Nb_Comp_4 33.76477  0.2710856  0.0975 -2.25950371 14.40075   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -2.4910745  0.0975 -3.78941633 20.63872   3.521924 0.9924713
#> Nb_Comp_6 35.25533         NA  0.0975          NA       NA   3.496074 0.9925265
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
#> attr(,"computed_nt")
#> [1] 6
#> 
#> [[2]]
#>                AIC     Q2cum_Y LimQ2_Y       Q2_Y  PRESS_Y      RSS_Y      R2_Y
#> Nb_Comp_0 82.01205          NA      NA         NA       NA 467.796667        NA
#> Nb_Comp_1 53.15173   0.8557983  0.0975  0.8557983 67.45709  35.742486 0.9235940
#> Nb_Comp_2 41.08283   0.8000384  0.0975 -0.3866793 49.56337  11.066606 0.9763431
#> Nb_Comp_3 32.06411   0.4610525  0.0975 -1.6952557 29.82733   4.418081 0.9905556
#> Nb_Comp_4 33.76477  -3.1852578  0.0975 -6.7656125 34.30911   4.309235 0.9907882
#> Nb_Comp_5 33.34373 -32.0034334  0.0975 -6.8856394 33.98108   3.521924 0.9924713
#> Nb_Comp_6 35.25533          NA  0.0975         NA       NA   3.496074 0.9925265
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
#> attr(,"computed_nt")
#> [1] 6
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
#> 1 2 3 
#> 0 0 1 
cvtable(sum3)
#> 
#> CV Q2 criterion:
#> 0 1 
#> 0 2 
#> 
#> CV Press criterion:
#> 1 2 3 4 
#> 0 0 1 1 
rm(list=c("XCornell","yCornell","bbb","bbb2","bbb3","bbb4"))
```
