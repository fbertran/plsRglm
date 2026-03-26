# Partial least squares Regression models with leave one out cross validation

This function implements Partial least squares Regression models with
leave one out cross validation for complete or incomplete datasets.

## Usage

``` r
plsR(object, ...)
# Default S3 method
plsRmodel(object, dataX, nt = 2, limQ2set = 0.0975, 
dataPredictY = dataX, modele = "pls", family = NULL, typeVC = "none", 
EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12), weights,
sparse = FALSE, sparseStop = TRUE, naive = FALSE,verbose=TRUE,...)
# S3 method for class 'formula'
plsRmodel(object, data, nt = 2, limQ2set = 0.0975,
dataPredictY, modele = "pls", family = NULL, typeVC = "none",
EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12), weights,
subset, contrasts = NULL, sparse = FALSE, sparseStop = TRUE, naive = FALSE,
verbose=TRUE,...)
PLS_lm(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, 
modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, 
scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, 
alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12),
weights,sparse=FALSE,sparseStop=FALSE,naive=FALSE,verbose=TRUE)
PLS_lm_formula(formula,data=NULL,nt=2,limQ2set=.0975,dataPredictY=dataX,
modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,
scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,
tol_Xi=10^(-12),weights,subset,contrasts=NULL,sparse=FALSE,
sparseStop=FALSE,naive=FALSE,verbose=TRUE)
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
  environment from which `plsR` is called.

- nt:

  number of components to be extracted

- limQ2set:

  limit value for the Q2

- dataPredictY:

  predictor(s) (testing) dataset

- modele:

  name of the PLS model to be fitted, only (`"pls"` available for this
  fonction.

- family:

  for the present moment the family argument is ignored and set thanks
  to the value of modele.

- typeVC:

  type of leave one out cross validation. Several procedures are
  available. If cross validation is required, one needs to selects the
  way of predicting the response for left out observations. For complete
  rows, without any missing value, there are two different ways of
  computing these predictions. As a consequence, for mixed datasets,
  with complete and incomplete rows, there are two ways of computing
  prediction : either predicts any row as if there were missing values
  in it (`missingdata`) or selects the prediction method accordingly to
  the completeness of the row (`adaptative`).

  `none`

  :   no cross validation

  `standard`

  :   as in SIMCA for datasets without any missing value. For datasets
      with any missing value, it is the as using `missingdata`

  `missingdata`

  :   all values predicted as those with missing values for datasets
      with any missing values

  `adaptative`

  :   predict a response value for an x with any missing value as those
      with missing values and for an x without any missing value as
      those without missing values.

- EstimXNA:

  only for `modele="pls"`. Set whether the missing X values have to be
  estimated.

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since non always possible for glm
  responses.

- pvals.expli:

  should individual p-values be reported to tune model selection ?

- alpha.pvals.expli:

  level of significance for predictors when pvals.expli=TRUE

- MClassed:

  number of missclassified cases, should only be used for binary
  responses

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

- sparse:

  should the coefficients of non-significant predictors
  (\<`alpha.pvals.expli`) be set to 0

- sparseStop:

  should component extraction stop when no significant predictors
  (\<`alpha.pvals.expli`) are found

- naive:

  Use the naive estimates for the Degrees of Freedom in plsR? Default is
  `FALSE`.

- verbose:

  should info messages be displayed ?

- ...:

  arguments to pass to `plsRmodel.default` or to `plsRmodel.formula`

## Details

There are several ways to deal with missing values that leads to
different computations of leave one out cross validation criteria.

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

The default estimator for Degrees of Freedom is the Kramer and
Sugiyama's one. Information criteria are computed accordingly to these
estimations. Naive Degrees of Freedom and Information Criteria are also
provided for comparison purposes. For more details, see N. Kraemer and
M. Sugiyama. (2011). The Degrees of Freedom of Partial Least Squares
Regression. *Journal of the American Statistical Association*, 106(494),
697-705, 2011.

## Value

- nr:

  Number of observations

- nc:

  Number of predictors

- nt:

  Number of requested components

- ww:

  raw weights (before L2-normalization)

- wwnorm:

  L2 normed weights (to be used with deflated matrices of predictor
  variables)

- wwetoile:

  modified weights (to be used with original matrix of predictor
  variables)

- tt:

  PLS components

- pp:

  loadings of the predictor variables

- CoeffC:

  coefficients of the PLS components

- uscores:

  scores of the response variable

- YChapeau:

  predicted response values for the dataX set

- residYChapeau:

  residuals of the deflated response on the standardized scale

- RepY:

  scaled response vector

- na.miss.Y:

  is there any NA value in the response vector

- YNA:

  indicatrix vector of missing values in RepY

- residY:

  deflated scaled response vector

- ExpliX:

  scaled matrix of predictors

- na.miss.X:

  is there any NA value in the predictor matrix

- XXNA:

  indicator of non-NA values in the predictor matrix

- residXX:

  deflated predictor matrix

- PredictY:

  response values with NA replaced with 0

- press.ind:

  individual PRESS value for each observation (scaled scale)

- press.tot:

  total PRESS value for all observations (scaled scale)

- family:

  glm family used to fit PLSGLR model

- ttPredictY:

  PLS components for the dataset on which prediction was requested

- typeVC:

  type of leave one out cross-validation used

- dataX:

  predictor values

- dataY:

  response values

- computed_nt:

  number of components that were computed

- CoeffCFull:

  matrix of the coefficients of the predictors

- CoeffConstante:

  value of the intercept (scaled scale)

- Std.Coeffs:

  Vector of standardized regression coefficients

- press.ind2:

  individual PRESS value for each observation (original scale)

- RSSresidY:

  residual sum of squares (scaled scale)

- Coeffs:

  Vector of regression coefficients (used with the original data scale)

- Yresidus:

  residuals of the PLS model

- RSS:

  residual sum of squares (original scale)

- residusY:

  residuals of the deflated response on the standardized scale

- AIC.std:

  AIC.std vs number of components (AIC computed for the standardized
  model

- AIC:

  AIC vs number of components

- optional:

  If the response is assumed to be binary:  
  i.e. `MClassed=TRUE`.

  `MissClassed`

  :   Number of miss classed results

  `Probs`

  :   "Probability" predicted by the model. These are not true
      probabilities since they may lay outside of \[0,1\]

  `Probs.trc`

  :   Probability predicted by the model and constrained to belong to
      \[0,1\]

- ttPredictFittedMissingY:

  Description of 'comp2'

- optional:

  If cross validation was requested:  
  i.e. `typeVC="standard"`, `typeVC="missingdata"` or
  `typeVC="adaptative"`.

  `R2residY`

  :   R2 coefficient value on the standardized scale

  `R2`

  :   R2 coefficient value on the original scale

  `press.tot2`

  :   total PRESS value for all observations (original scale)

  `Q2`

  :   Q2 value (standardized scale)

  `limQ2`

  :   limit of the Q2 value

  `Q2_2`

  :   Q2 value (original scale)

  `Q2cum`

  :   cumulated Q2 (standardized scale)

  `Q2cum_2`

  :   cumulated Q2 (original scale)

- InfCrit:

  table of Information Criteria

- Std.ValsPredictY:

  predicted response values for supplementary dataset (standardized
  scale)

- ValsPredictY:

  predicted response values for supplementary dataset (original scale)

- Std.XChapeau:

  estimated values for missing values in the predictor matrix
  (standardized scale)

- XXwotNA:

  predictor matrix with missing values replaced with 0

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

Use [`cv.plsR`](https://fbertran.github.io/plsRglm/reference/cv.plsR.md)
to cross-validate the plsRglm models and
[`bootpls`](https://fbertran.github.io/plsRglm/reference/bootpls.md) to
bootstrap them.

## See also

See also
[`plsRglm`](https://fbertran.github.io/plsRglm/reference/plsRglm.md) to
fit PLSGLR models.

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

#maximum 6 components could be extracted from this dataset
#trying 10 to trigger automatic stopping criterion
modpls10<-plsR(yCornell,XCornell,10)
#> ____************************************************____
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
modpls10
#> Number of required components:
#> [1] 10
#> Number of successfully computed components:
#> [1] 6
#> Coefficients:
#>                  [,1]
#> Intercept  88.7107982
#> X1        -54.3905712
#> X2         -2.7879678
#> X3         52.5411315
#> X4        -11.5306977
#> X5         -0.9605822
#> X6         11.5900307
#> X7         28.2104803
#> Information criteria and Fit statistics:
#>                AIC      RSS_Y      R2_Y R2_residY  RSS_residY    AIC.std
#> Nb_Comp_0 82.01205 467.796667        NA        NA 11.00000000  37.010388
#> Nb_Comp_1 53.15173  35.742486 0.9235940 0.9235940  0.84046633   8.150064
#> Nb_Comp_2 41.08283  11.066606 0.9763431 0.9763431  0.26022559  -3.918831
#> Nb_Comp_3 32.06411   4.418081 0.9905556 0.9905556  0.10388893 -12.937550
#> Nb_Comp_4 33.76477   4.309235 0.9907882 0.9907882  0.10132947 -11.236891
#> Nb_Comp_5 33.34373   3.521924 0.9924713 0.9924713  0.08281624 -11.657929
#> Nb_Comp_6 35.25533   3.496074 0.9925265 0.9925265  0.08220840  -9.746328
#>            DoF.dof sigmahat.dof    AIC.dof    BIC.dof GMDL.dof DoF.naive
#> Nb_Comp_0 1.000000    6.5212706 46.0708838 47.7893514 27.59461         1
#> Nb_Comp_1 2.740749    1.8665281  4.5699686  4.9558156 21.34020         2
#> Nb_Comp_2 5.085967    1.1825195  2.1075461  2.3949331 27.40202         3
#> Nb_Comp_3 5.121086    0.7488308  0.8467795  0.9628191 24.40842         4
#> Nb_Comp_4 5.103312    0.7387162  0.8232505  0.9357846 24.23105         5
#> Nb_Comp_5 6.006316    0.7096382  0.7976101  0.9198348 28.21184         6
#> Nb_Comp_6 7.000001    0.7633342  0.9711319  1.1359499 33.18347         7
#>           sigmahat.naive  AIC.naive  BIC.naive GMDL.naive
#> Nb_Comp_0      6.5212706 46.0708838 47.7893514   27.59461
#> Nb_Comp_1      1.8905683  4.1699567  4.4588195   18.37545
#> Nb_Comp_2      1.1088836  1.5370286  1.6860917   17.71117
#> Nb_Comp_3      0.7431421  0.7363469  0.8256118   19.01033
#> Nb_Comp_4      0.7846050  0.8721072  0.9964867   24.16510
#> Nb_Comp_5      0.7661509  0.8804809  1.0227979   28.64206
#> Nb_Comp_6      0.8361907  1.1070902  1.3048716   33.63927

#With iterated leave one out CV PRESS
modpls6cv<-plsR(Y~.,data=Cornell,6,typeVC="standard")
#> ____************************************************____
#> ____TypeVC____ standard ____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
modpls6cv
#> Number of required components:
#> [1] 6
#> Number of successfully computed components:
#> [1] 6
#> Coefficients:
#>                  [,1]
#> Intercept  88.7107982
#> X1        -54.3905712
#> X2         -2.7879678
#> X3         52.5411315
#> X4        -11.5306977
#> X5         -0.9605822
#> X6         11.5900307
#> X7         28.2104803
#> Leave one out cross validated PRESS, Information criteria and Fit statistics:
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
cv.modpls<-cv.plsR(Y~.,data=Cornell,6,NK=100, verbose=FALSE)
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
#> NK: 1,  2,  3,  4,  5,  6,  7,  8,  9,  10
#> NK: 11,  12,  13,  14,  15,  16,  17,  18,  19,  20
#> NK: 21,  22,  23,  24,  25,  26,  27,  28,  29,  30
#> NK: 31,  32,  33,  34,  35,  36,  37,  38,  39,  40
#> NK: 41,  42,  43,  44,  45,  46,  47,  48,  49,  50
#> NK: 51,  52,  53,  54,  55,  56,  57,  58,  59,  60
#> NK: 61,  62,  63,  64,  65,  66,  67,  68,  69,  70
#> NK: 71,  72,  73,  74,  75,  76,  77,  78,  79,  80
#> NK: 81,  82,  83,  84,  85,  86,  87,  88,  89,  90
#> NK: 91,  92,  93,  94,  95,  96,  97,  98,  99,  100
#> 
#> CV Q2 criterion:
#>  0  1  2 
#>  0 86 14 
#> 
#> CV Press criterion:
#>  1  2  3  4  5 
#>  0  1 31 52 16 
plot(res.cv.modpls)


rm(list=c("XCornell","yCornell","modpls10","modpls6cv"))

# \donttest{
#A binary response example
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
modpls.aze <- plsR(yaze_compl,Xaze_compl,10,MClassed=TRUE,typeVC="standard")
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
modpls.aze
#> Number of required components:
#> [1] 10
#> Number of successfully computed components:
#> [1] 10
#> Coefficients:
#>                   [,1]
#> Intercept  0.308019808
#> D2S138    -0.131218617
#> D18S61     0.450219840
#> D16S422   -0.183848373
#> D17S794    0.269084083
#> D6S264     0.105061098
#> D14S65    -0.052837918
#> D18S53     0.008489326
#> D17S790   -0.213122117
#> D1S225     0.046277290
#> D3S1282   -0.095666162
#> D9S179     0.054547887
#> D5S430    -0.126491043
#> D8S283     0.106373432
#> D11S916    0.111623381
#> D2S159     0.056759714
#> D16S408    0.010288859
#> D5S346     0.233674850
#> D10S191    0.010715856
#> D13S173    0.074148740
#> D6S275    -0.123145693
#> D15S127    0.064566148
#> D1S305     0.190500469
#> D4S394    -0.142585807
#> D20S107   -0.184483600
#> D1S197    -0.284373695
#> D1S207     0.186728597
#> D10S192    0.195516079
#> D3S1283   -0.096309755
#> D4S414     0.017960975
#> D8S264     0.121051206
#> D22S928   -0.049091794
#> TP53      -0.391965015
#> D9S171    -0.012315197
#> Leave one out cross validated PRESS, Information criteria and Fit statistics:
#>                 AIC      Q2cum_Y LimQ2_Y        Q2_Y  PRESS_Y    RSS_Y
#> Nb_Comp_0  154.6179           NA      NA          NA       NA 25.91346
#> Nb_Comp_1  126.4083  -0.09840016  0.0975 -0.09840016 28.46335 19.38086
#> Nb_Comp_2  119.3375  -0.19018163  0.0975 -0.08355923 21.00031 17.76209
#> Nb_Comp_3  114.2313  -0.77332918  0.0975 -0.48996518 26.46489 16.58896
#> Nb_Comp_4  112.3463  -1.64635954  0.0975 -0.49231150 24.75590 15.98071
#> Nb_Comp_5  113.2362  -2.74242209  0.0975 -0.41417749 22.59955 15.81104
#> Nb_Comp_6  114.7620  -4.46009228  0.0975 -0.45897286 23.06788 15.73910
#> Nb_Comp_7  116.5264  -7.36664482  0.0975 -0.53232663 24.11744 15.70350
#> Nb_Comp_8  118.4601 -11.80011367  0.0975 -0.52989806 24.02475 15.69348
#> Nb_Comp_9  120.4452 -17.90787273  0.0975 -0.47716444 23.18185 15.69123
#> Nb_Comp_10 122.4395 -26.50536212  0.0975 -0.45470421 22.82610 15.69037
#>                 R2_Y MissClassed R2_residY RSS_residY PRESS_residY   Q2_residY
#> Nb_Comp_0         NA          49        NA  103.00000           NA          NA
#> Nb_Comp_1  0.2520929          27 0.2520929   77.03443    113.13522 -0.09840016
#> Nb_Comp_2  0.3145613          25 0.3145613   70.60018     83.47137 -0.08355923
#> Nb_Comp_3  0.3598323          27 0.3598323   65.93728    105.19181 -0.48996518
#> Nb_Comp_4  0.3833049          23 0.3833049   63.51960     98.39895 -0.49231150
#> Nb_Comp_5  0.3898523          22 0.3898523   62.84522     89.82798 -0.41417749
#> Nb_Comp_6  0.3926285          21 0.3926285   62.55927     91.68947 -0.45897286
#> Nb_Comp_7  0.3940024          20 0.3940024   62.41775     95.86123 -0.53232663
#> Nb_Comp_8  0.3943888          20 0.3943888   62.37795     95.49280 -0.52989806
#> Nb_Comp_9  0.3944758          19 0.3944758   62.36900     92.14249 -0.47716444
#> Nb_Comp_10 0.3945088          19 0.3945088   62.36560     90.72844 -0.45470421
#>             LimQ2 Q2cum_residY  AIC.std  DoF.dof sigmahat.dof   AIC.dof
#> Nb_Comp_0      NA           NA 298.1344  1.00000    0.5015845 0.2540061
#> Nb_Comp_1  0.0975  -0.09840016 269.9248 22.55372    0.4848429 0.2883114
#> Nb_Comp_2  0.0975  -0.19018163 262.8540 27.31542    0.4781670 0.2908950
#> Nb_Comp_3  0.0975  -0.77332918 257.7478 30.52370    0.4719550 0.2902572
#> Nb_Comp_4  0.0975  -1.64635954 255.8628 34.00000    0.4744263 0.3008285
#> Nb_Comp_5  0.0975  -2.74242209 256.7527 34.00000    0.4719012 0.2976347
#> Nb_Comp_6  0.0975  -4.46009228 258.2785 34.00000    0.4708264 0.2962804
#> Nb_Comp_7  0.0975  -7.36664482 260.0429 33.71066    0.4693382 0.2937976
#> Nb_Comp_8  0.0975 -11.80011367 261.9766 34.00000    0.4701436 0.2954217
#> Nb_Comp_9  0.0975 -17.90787273 263.9617 33.87284    0.4696894 0.2945815
#> Nb_Comp_10 0.0975 -26.50536212 265.9560 34.00000    0.4700970 0.2953632
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

#Direct access to not cross-validated values
modpls.aze$AIC
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]    [,7]     [,8]
#> [1,] 154.6179 126.4083 119.3375 114.2313 112.3463 113.2362 114.762 116.5264
#>          [,9]    [,10]    [,11]
#> [1,] 118.4601 120.4452 122.4395
modpls.aze$AIC.std
#>          [,1]     [,2]    [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 298.1344 269.9248 262.854 257.7478 255.8628 256.7527 258.2785 260.0429
#>          [,9]    [,10]   [,11]
#> [1,] 261.9766 263.9617 265.956
modpls.aze$MissClassed
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,]   49   27   25   27   23   22   21   20   20    19    19

#Raw predicted values (not really probabily since not constrained in [0,1]
modpls.aze$Probs
#>          [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
#> 1   0.4711538  0.46105744  0.63458141  0.67961627  0.69452246  0.64534767
#> 2   0.4711538  0.26911816  0.26581497  0.16989268  0.11760783  0.18096700
#> 3   0.4711538 -0.09080494 -0.05104846 -0.17166916 -0.21455242 -0.21725391
#> 4   0.4711538  0.36370490  0.54112657  0.50724821  0.55508565  0.57773785
#> 5   0.4711538 -0.04408124  0.07399231 -0.07129909 -0.24018962 -0.23445282
#> 6   0.4711538 -0.03776963  0.17275288  0.01806190 -0.02597539 -0.06284454
#> 7   0.4711538 -0.06930728 -0.19928456 -0.09137261  0.01116043  0.06506517
#> 8   0.4711538  0.27158233  0.24933653  0.11611522  0.12804487  0.04118115
#> 9   0.4711538  0.76949497  0.60296556  0.47237794  0.51581382  0.49885092
#> 10  0.4711538  0.22096539  0.34482052  0.34660816  0.38580378  0.43528451
#> 11  0.4711538  0.87147914  0.84865348  0.76372713  0.73582307  0.76725258
#> 12  0.4711538  0.79792975  0.67828859  0.73747065  0.67844373  0.67908585
#> 13  0.4711538  0.09432664 -0.04344681  0.10780023  0.22488457  0.26144110
#> 14  0.4711538  0.28543133  0.29293086  0.37385135  0.37961001  0.30207755
#> 15  0.4711538  0.30637401  0.27816310  0.18074751  0.01510565  0.05074255
#> 16  0.4711538  0.12893721 -0.07276258 -0.05146556 -0.09988241 -0.06790398
#> 17  0.4711538  0.59910292  0.41302582  0.40055026  0.32477692  0.32429673
#> 18  0.4711538  0.60665328  0.51461671  0.70351041  0.63093215  0.60232625
#> 19  0.4711538  0.18381206  0.36596047  0.33591603  0.25289460  0.21859872
#> 20  0.4711538  0.28422822  0.15202852  0.29980632  0.42075827  0.43463142
#> 21  0.4711538  0.35982960  0.40300075  0.63220247  0.58056075  0.55273462
#> 22  0.4711538  0.31574837  0.28422517  0.37116719  0.27156145  0.25529246
#> 23  0.4711538  0.41682757  0.36900849  0.23791176  0.25730930  0.24221472
#> 24  0.4711538  0.30288056  0.15972272  0.19362318  0.07194768  0.07250435
#> 25  0.4711538  0.29650015  0.48867070  0.61025747  0.59737342  0.67704212
#> 26  0.4711538  0.23008536  0.32001822  0.15862645  0.26312675  0.22513847
#> 27  0.4711538  0.67526360  0.68123526  0.58796740  0.51309143  0.44381568
#> 28  0.4711538  0.15222775  0.13544964  0.15605402  0.15868232  0.10574096
#> 29  0.4711538  0.43138914  0.29576924  0.29706087  0.35294305  0.40257625
#> 30  0.4711538  0.13910581  0.26763382  0.10182481  0.12169881  0.13543560
#> 31  0.4711538  0.40295972  0.43810789  0.28684877  0.41632594  0.45388666
#> 32  0.4711538  0.58422149  0.44366239  0.16615851  0.15367980  0.18291151
#> 33  0.4711538  0.69889100  0.72592310  0.57845537  0.50185886  0.51841164
#> 34  0.4711538  0.35960908  0.24234167  0.09364940  0.08428214  0.10528276
#> 35  0.4711538  0.27914959  0.03731133 -0.08896074 -0.06232370 -0.08231459
#> 36  0.4711538  0.38865989  0.39024480  0.44138316  0.47508801  0.42329842
#> 37  0.4711538  0.62200134  0.42145828  0.38142396  0.29675933  0.28947211
#> 38  0.4711538  0.41311694  0.19970983  0.16702613  0.17059545  0.17073272
#> 39  0.4711538  0.31755422  0.28395547  0.17609314  0.23875966  0.25763504
#> 40  0.4711538  0.62628933  0.51627261  0.52025889  0.47789760  0.47304606
#> 41  0.4711538  0.14894845  0.14069540  0.13906223  0.05976750  0.13670893
#> 42  0.4711538  0.64041121  0.49727655  0.49380105  0.53239359  0.51394469
#> 43  0.4711538  0.38696544  0.54930653  0.62650411  0.65244562  0.56755351
#> 44  0.4711538  0.24204195  0.05825611  0.02230584 -0.01790809 -0.03785626
#> 45  0.4711538  0.10349021  0.14957660  0.16304594  0.15564790  0.17065395
#> 46  0.4711538  0.63322787  0.64625855  0.55541948  0.65203351  0.63670168
#> 47  0.4711538  0.20557889  0.23864853  0.24328712  0.13063078  0.09743813
#> 48  0.4711538  0.32352238  0.34894312  0.21162810  0.20487572  0.16461876
#> 49  0.4711538  0.64888519  0.52290405  0.50926772  0.62061797  0.59597941
#> 50  0.4711538  0.44153005  0.49754241  0.32749149  0.24840605  0.32456388
#> 51  0.4711538  0.32562433  0.23887414  0.26764033  0.24950898  0.30432045
#> 52  0.4711538 -0.23250098 -0.28713647 -0.09216174 -0.12709475 -0.18324647
#> 53  0.4711538  0.53388610  0.47710127  0.60836140  0.48273912  0.43334108
#> 54  0.4711538  0.64191356  0.44931093  0.46371798  0.45275305  0.46653696
#> 55  0.4711538  0.05279255  0.06829351  0.15306458  0.25200214  0.21249173
#> 56  0.4711538  0.59808020  0.64333345  0.53741245  0.64108173  0.57876914
#> 57  0.4711538  0.53093147  0.62138656  0.92046148  0.93004391  0.95130430
#> 58  0.4711538  0.64943097  0.57141374  0.66800038  0.64835800  0.65566321
#> 59  0.4711538  0.42541400  0.43027409  0.30117492  0.36183156  0.29992796
#> 60  0.4711538  0.24537249  0.29963849  0.42931558  0.51048830  0.58927966
#> 61  0.4711538  0.64269314  0.62785202  0.75163561  0.68045267  0.67000184
#> 62  0.4711538  0.51277761  0.60877778  0.75493489  0.66735142  0.63862193
#> 63  0.4711538  0.53377378  0.53228159  0.56245626  0.58414332  0.61176055
#> 64  0.4711538  0.79099666  0.90572246  0.92244949  0.93001276  0.93454809
#> 65  0.4711538  0.73768777  0.61339931  0.72362105  0.70536287  0.69970096
#> 66  0.4711538  0.70767466  0.53408924  0.50675818  0.52181506  0.54559559
#> 67  0.4711538  0.96312042  1.17012215  1.08116795  1.22497425  1.21728258
#> 68  0.4711538  0.31575995  0.57179559  0.77297374  0.78532935  0.78484987
#> 69  0.4711538  0.69505872  0.78176548  0.74300700  0.72711033  0.70750770
#> 70  0.4711538  0.72276362  0.90232185  0.89364576  0.84428623  0.92659977
#> 71  0.4711538  0.50950893  0.39503961  0.45591683  0.38297596  0.35086204
#> 72  0.4711538  0.14720074  0.13538571 -0.04473829 -0.05529233  0.02748516
#> 73  0.4711538  0.49275110  0.44937896  0.41856171  0.62470016  0.61654596
#> 74  0.4711538  0.65674324  0.69439259  0.75479685  0.88511667  0.92560996
#> 75  0.4711538  0.68716407  0.57541914  0.59945962  0.54581071  0.55228791
#> 76  0.4711538  0.54839542  0.50508123  0.52627725  0.55765709  0.52543838
#> 77  0.4711538  0.77317727  0.79812663  0.93073165  1.10301473  1.08723742
#> 78  0.4711538  0.85322027  0.76128342  0.81061207  0.85796753  0.87947603
#> 79  0.4711538  0.81659194  0.90228252  0.80744839  0.70383361  0.68468090
#> 80  0.4711538  0.55964651  0.44326524  0.39507689  0.36149039  0.32071350
#> 81  0.4711538  0.87105473  0.86695796  0.89177640  0.74816339  0.69831750
#> 82  0.4711538  0.47715869  0.68930595  0.71280202  0.73606020  0.78321326
#> 83  0.4711538  0.80974821  0.87138779  0.97466313  0.93082943  0.95560886
#> 84  0.4711538  0.67739807  0.85743609  0.98894432  0.96011041  0.90800271
#> 85  0.4711538  0.57131444  0.34250950  0.33855791  0.31118498  0.31383288
#> 86  0.4711538  0.84958765  0.97611051  0.93090902  0.91560248  0.86222031
#> 87  0.4711538  0.57644613  0.41449248  0.48714466  0.54811918  0.57041511
#> 88  0.4711538  0.75932310  0.71214369  0.52234742  0.59011684  0.59023780
#> 89  0.4711538  0.53031516  0.47090892  0.42433053  0.38847912  0.39218094
#> 90  0.4711538  0.76770402  1.07649866  1.00864429  1.06363018  1.09017457
#> 91  0.4711538  0.38643842  0.37696993  0.44452861  0.49450298  0.46628856
#> 92  0.4711538  0.92591633  1.03707888  0.96084369  0.95688931  0.93400393
#> 93  0.4711538  0.66726042  0.89247800  0.87390628  0.87335977  0.95801535
#> 94  0.4711538  0.32634752  0.41373057  0.48066349  0.67273089  0.62115180
#> 95  0.4711538  0.50472276  0.77159222  0.71730564  0.62350221  0.64335334
#> 96  0.4711538  0.34622269  0.33150717  0.49412629  0.44574013  0.46889514
#> 97  0.4711538  0.55805257  0.50280611  0.58541977  0.52239953  0.53556273
#> 98  0.4711538  0.78090964  0.73429355  0.79385683  0.86651416  0.88151677
#> 99  0.4711538  0.21116352  0.10917861  0.02565398  0.18342015  0.15222876
#> 100 0.4711538  0.66672702  0.78264411  0.86306662  0.75733969  0.77632472
#> 101 0.4711538  0.45317545  0.50149615  0.62617428  0.70904267  0.78134354
#> 102 0.4711538  0.74435376  0.66135006  0.72568147  0.70203564  0.77593538
#> 103 0.4711538  0.34690226  0.56605434  0.52782336  0.50951738  0.46795757
#> 104 0.4711538  0.69496014  0.80515138  0.78871059  0.78008789  0.78042831
#>            [,7]         [,8]         [,9]       [,10]       [,11]
#> 1    0.64037279  0.627340571  0.651243676  0.65354280  0.65838797
#> 2    0.21304385  0.202666528  0.206463548  0.21046038  0.20470356
#> 3   -0.22444089 -0.193652144 -0.201652437 -0.20167289 -0.20081532
#> 4    0.58413761  0.600377920  0.608768219  0.60784745  0.60193905
#> 5   -0.26327049 -0.311781941 -0.310765976 -0.31066830 -0.31401382
#> 6   -0.10341096 -0.076840858 -0.080016598 -0.08436785 -0.08777135
#> 7    0.10261786  0.132750517  0.144750243  0.13944940  0.14173177
#> 8    0.04809780  0.063736599  0.063261540  0.07570783  0.07715150
#> 9    0.44543808  0.444943670  0.447021225  0.44582742  0.44748513
#> 10   0.43490588  0.407005593  0.413663760  0.41349216  0.41350336
#> 11   0.78695284  0.777623618  0.783734315  0.78262765  0.77949531
#> 12   0.65654818  0.642154505  0.633436560  0.63197252  0.62875264
#> 13   0.21708341  0.187509114  0.186533541  0.17822088  0.17842513
#> 14   0.28651410  0.260501290  0.279081223  0.27196942  0.26996735
#> 15   0.05784774  0.095949877  0.090894676  0.08865039  0.09080522
#> 16  -0.06239212 -0.039505632 -0.023469898 -0.02045198 -0.02715190
#> 17   0.33482409  0.336193547  0.335468481  0.33501592  0.33670840
#> 18   0.59257402  0.580678556  0.581449676  0.58013547  0.57927640
#> 19   0.24272344  0.246701307  0.240991178  0.23822863  0.23157123
#> 20   0.40555564  0.386074751  0.400419290  0.40843493  0.40860630
#> 21   0.53582205  0.549336181  0.539598759  0.54139679  0.54098845
#> 22   0.27214775  0.265323031  0.262889367  0.27124503  0.27372049
#> 23   0.26310951  0.254264842  0.244111736  0.24044626  0.24169731
#> 24   0.10789143  0.136298989  0.142908568  0.14863967  0.15404679
#> 25   0.62754156  0.650367844  0.644944061  0.64644697  0.64416177
#> 26   0.18477380  0.190402191  0.199245619  0.20152378  0.20582432
#> 27   0.44358238  0.454283166  0.446210009  0.43311399  0.43568840
#> 28   0.13769251  0.118081488  0.117589648  0.11048814  0.11155803
#> 29   0.42759376  0.424018002  0.417256036  0.42362830  0.42160066
#> 30   0.10952676  0.168692282  0.174041713  0.17387211  0.17480429
#> 31   0.45876287  0.432435790  0.425259491  0.43329412  0.43856907
#> 32   0.11742399  0.116981896  0.108778824  0.10636129  0.10689920
#> 33   0.47697376  0.450738829  0.452062165  0.44881973  0.44992624
#> 34   0.09863470  0.095862956  0.092253997  0.09980718  0.10084230
#> 35  -0.08390506 -0.085155814 -0.086690885 -0.08262604 -0.08272575
#> 36   0.45381739  0.433667887  0.449594388  0.44854107  0.44960414
#> 37   0.30022320  0.311627620  0.315020647  0.31958491  0.32067580
#> 38   0.15959719  0.160997289  0.151222415  0.15079738  0.14817778
#> 39   0.27800144  0.253643110  0.255071837  0.25875102  0.25674198
#> 40   0.51555833  0.530347073  0.528402618  0.52790558  0.52538704
#> 41   0.09275154  0.101869980  0.090604733  0.09445259  0.09150467
#> 42   0.50207194  0.483136866  0.487491071  0.48860245  0.49212211
#> 43   0.57336835  0.573091866  0.546609274  0.54001848  0.54366888
#> 44  -0.01451487 -0.004369241 -0.002647953 -0.00844687 -0.01073887
#> 45   0.19657318  0.164988597  0.170522098  0.16781405  0.16715155
#> 46   0.65287559  0.645084587  0.651535848  0.65177845  0.65178893
#> 47   0.12452036  0.119913760  0.117053912  0.11873700  0.11591426
#> 48   0.18275993  0.166998459  0.160605969  0.16869180  0.16855238
#> 49   0.58006784  0.613009102  0.608873759  0.60785525  0.61021103
#> 50   0.33208894  0.314978063  0.311974842  0.30517646  0.31174053
#> 51   0.32816749  0.318321881  0.320797741  0.31281691  0.31540371
#> 52  -0.19584259 -0.184768370 -0.177691131 -0.18225253 -0.18172040
#> 53   0.40794213  0.415136651  0.416044038  0.42408314  0.42075636
#> 54   0.46480524  0.479688395  0.471302624  0.46613648  0.46656474
#> 55   0.23097197  0.213772954  0.189358085  0.18930563  0.19050731
#> 56   0.57331393  0.574558380  0.579775660  0.57922184  0.57800858
#> 57   0.96183872  0.975414458  0.968855613  0.96385549  0.96028958
#> 58   0.64361788  0.631632631  0.641013340  0.63785489  0.64045306
#> 59   0.28643229  0.289148607  0.278791345  0.28784917  0.28782462
#> 60   0.55978204  0.559023840  0.561630138  0.55671889  0.55726966
#> 61   0.65787202  0.659412106  0.650588243  0.64930310  0.65091617
#> 62   0.60705201  0.603977438  0.589826600  0.58812271  0.58793743
#> 63   0.63813827  0.617717657  0.614454772  0.61166047  0.60985403
#> 64   1.00166279  1.008968255  1.008852814  1.01022555  1.00666729
#> 65   0.69297263  0.697694744  0.689810628  0.69214943  0.69038662
#> 66   0.53663407  0.532369158  0.535934708  0.53171176  0.53102146
#> 67   1.22111150  1.202569147  1.194030443  1.19310867  1.19228184
#> 68   0.80071187  0.812462410  0.796883689  0.80360414  0.80689017
#> 69   0.74791867  0.774669023  0.764357440  0.76435231  0.76730276
#> 70   0.93952180  0.926356875  0.916496158  0.91637035  0.92383169
#> 71   0.31179970  0.319584982  0.342996678  0.34334588  0.34349096
#> 72   0.08069763  0.064883087  0.053282640  0.04500807  0.04571864
#> 73   0.63914960  0.656024336  0.647841400  0.64971036  0.64930911
#> 74   0.94540783  0.945531101  0.947433351  0.95548576  0.95086945
#> 75   0.56609663  0.581051527  0.586510546  0.58350518  0.58195688
#> 76   0.49807985  0.507124137  0.503700114  0.50352891  0.50289031
#> 77   1.06463806  1.090636226  1.099163425  1.09692704  1.09684987
#> 78   0.88947890  0.914074686  0.921993535  0.92611779  0.92473772
#> 79   0.70672170  0.689699520  0.693322761  0.69306283  0.69031633
#> 80   0.30181332  0.296560012  0.282907508  0.28672758  0.28334638
#> 81   0.69871492  0.692603330  0.692711040  0.69448149  0.69353087
#> 82   0.73754433  0.755515607  0.753806256  0.76075515  0.76006115
#> 83   0.95554583  0.981114506  0.967899112  0.96789876  0.97167928
#> 84   0.92460814  0.928255751  0.934569121  0.93157477  0.93560641
#> 85   0.31932547  0.326702214  0.334191767  0.33780085  0.33863677
#> 86   0.82950069  0.806588554  0.823792983  0.82202715  0.81981455
#> 87   0.56750119  0.560439488  0.549720853  0.54609905  0.54656224
#> 88   0.57484476  0.588841816  0.582284530  0.57334028  0.57098492
#> 89   0.44028895  0.455790854  0.462215323  0.46149391  0.46635298
#> 90   1.08703587  1.113909884  1.123382027  1.12669615  1.12559691
#> 91   0.50591629  0.528734399  0.542033081  0.53789411  0.54025902
#> 92   0.95078882  0.944122848  0.933909143  0.93429933  0.92818920
#> 93   0.95048831  0.984993815  0.999870929  0.99916201  0.99970721
#> 94   0.63222931  0.617803323  0.604620965  0.60316670  0.59998131
#> 95   0.60173453  0.598838861  0.615665753  0.61062288  0.61003734
#> 96   0.42461082  0.393735123  0.382307431  0.38698306  0.38928062
#> 97   0.55493331  0.545341350  0.548945997  0.55196226  0.55110104
#> 98   0.84875293  0.849088124  0.833445596  0.82896512  0.82882885
#> 99   0.12817884  0.125184677  0.142524174  0.14250387  0.14467219
#> 100  0.76010896  0.740234308  0.744642249  0.75159300  0.75722807
#> 101  0.80270481  0.778624924  0.777761525  0.78102072  0.77766043
#> 102  0.77647730  0.755860583  0.764875857  0.76840548  0.77181270
#> 103  0.50588630  0.488855008  0.495423757  0.50188675  0.50526664
#> 104  0.81071639  0.804180721  0.825464815  0.81861016  0.81635531
#Truncated to [0;1] predicted values (true probabilities)
modpls.aze$Probs.trc
#>          [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
#> 1   0.4711538 0.46105744 0.63458141 0.67961627 0.69452246 0.64534767 0.64037279
#> 2   0.4711538 0.26911816 0.26581497 0.16989268 0.11760783 0.18096700 0.21304385
#> 3   0.4711538 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> 4   0.4711538 0.36370490 0.54112657 0.50724821 0.55508565 0.57773785 0.58413761
#> 5   0.4711538 0.00000000 0.07399231 0.00000000 0.00000000 0.00000000 0.00000000
#> 6   0.4711538 0.00000000 0.17275288 0.01806190 0.00000000 0.00000000 0.00000000
#> 7   0.4711538 0.00000000 0.00000000 0.00000000 0.01116043 0.06506517 0.10261786
#> 8   0.4711538 0.27158233 0.24933653 0.11611522 0.12804487 0.04118115 0.04809780
#> 9   0.4711538 0.76949497 0.60296556 0.47237794 0.51581382 0.49885092 0.44543808
#> 10  0.4711538 0.22096539 0.34482052 0.34660816 0.38580378 0.43528451 0.43490588
#> 11  0.4711538 0.87147914 0.84865348 0.76372713 0.73582307 0.76725258 0.78695284
#> 12  0.4711538 0.79792975 0.67828859 0.73747065 0.67844373 0.67908585 0.65654818
#> 13  0.4711538 0.09432664 0.00000000 0.10780023 0.22488457 0.26144110 0.21708341
#> 14  0.4711538 0.28543133 0.29293086 0.37385135 0.37961001 0.30207755 0.28651410
#> 15  0.4711538 0.30637401 0.27816310 0.18074751 0.01510565 0.05074255 0.05784774
#> 16  0.4711538 0.12893721 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> 17  0.4711538 0.59910292 0.41302582 0.40055026 0.32477692 0.32429673 0.33482409
#> 18  0.4711538 0.60665328 0.51461671 0.70351041 0.63093215 0.60232625 0.59257402
#> 19  0.4711538 0.18381206 0.36596047 0.33591603 0.25289460 0.21859872 0.24272344
#> 20  0.4711538 0.28422822 0.15202852 0.29980632 0.42075827 0.43463142 0.40555564
#> 21  0.4711538 0.35982960 0.40300075 0.63220247 0.58056075 0.55273462 0.53582205
#> 22  0.4711538 0.31574837 0.28422517 0.37116719 0.27156145 0.25529246 0.27214775
#> 23  0.4711538 0.41682757 0.36900849 0.23791176 0.25730930 0.24221472 0.26310951
#> 24  0.4711538 0.30288056 0.15972272 0.19362318 0.07194768 0.07250435 0.10789143
#> 25  0.4711538 0.29650015 0.48867070 0.61025747 0.59737342 0.67704212 0.62754156
#> 26  0.4711538 0.23008536 0.32001822 0.15862645 0.26312675 0.22513847 0.18477380
#> 27  0.4711538 0.67526360 0.68123526 0.58796740 0.51309143 0.44381568 0.44358238
#> 28  0.4711538 0.15222775 0.13544964 0.15605402 0.15868232 0.10574096 0.13769251
#> 29  0.4711538 0.43138914 0.29576924 0.29706087 0.35294305 0.40257625 0.42759376
#> 30  0.4711538 0.13910581 0.26763382 0.10182481 0.12169881 0.13543560 0.10952676
#> 31  0.4711538 0.40295972 0.43810789 0.28684877 0.41632594 0.45388666 0.45876287
#> 32  0.4711538 0.58422149 0.44366239 0.16615851 0.15367980 0.18291151 0.11742399
#> 33  0.4711538 0.69889100 0.72592310 0.57845537 0.50185886 0.51841164 0.47697376
#> 34  0.4711538 0.35960908 0.24234167 0.09364940 0.08428214 0.10528276 0.09863470
#> 35  0.4711538 0.27914959 0.03731133 0.00000000 0.00000000 0.00000000 0.00000000
#> 36  0.4711538 0.38865989 0.39024480 0.44138316 0.47508801 0.42329842 0.45381739
#> 37  0.4711538 0.62200134 0.42145828 0.38142396 0.29675933 0.28947211 0.30022320
#> 38  0.4711538 0.41311694 0.19970983 0.16702613 0.17059545 0.17073272 0.15959719
#> 39  0.4711538 0.31755422 0.28395547 0.17609314 0.23875966 0.25763504 0.27800144
#> 40  0.4711538 0.62628933 0.51627261 0.52025889 0.47789760 0.47304606 0.51555833
#> 41  0.4711538 0.14894845 0.14069540 0.13906223 0.05976750 0.13670893 0.09275154
#> 42  0.4711538 0.64041121 0.49727655 0.49380105 0.53239359 0.51394469 0.50207194
#> 43  0.4711538 0.38696544 0.54930653 0.62650411 0.65244562 0.56755351 0.57336835
#> 44  0.4711538 0.24204195 0.05825611 0.02230584 0.00000000 0.00000000 0.00000000
#> 45  0.4711538 0.10349021 0.14957660 0.16304594 0.15564790 0.17065395 0.19657318
#> 46  0.4711538 0.63322787 0.64625855 0.55541948 0.65203351 0.63670168 0.65287559
#> 47  0.4711538 0.20557889 0.23864853 0.24328712 0.13063078 0.09743813 0.12452036
#> 48  0.4711538 0.32352238 0.34894312 0.21162810 0.20487572 0.16461876 0.18275993
#> 49  0.4711538 0.64888519 0.52290405 0.50926772 0.62061797 0.59597941 0.58006784
#> 50  0.4711538 0.44153005 0.49754241 0.32749149 0.24840605 0.32456388 0.33208894
#> 51  0.4711538 0.32562433 0.23887414 0.26764033 0.24950898 0.30432045 0.32816749
#> 52  0.4711538 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> 53  0.4711538 0.53388610 0.47710127 0.60836140 0.48273912 0.43334108 0.40794213
#> 54  0.4711538 0.64191356 0.44931093 0.46371798 0.45275305 0.46653696 0.46480524
#> 55  0.4711538 0.05279255 0.06829351 0.15306458 0.25200214 0.21249173 0.23097197
#> 56  0.4711538 0.59808020 0.64333345 0.53741245 0.64108173 0.57876914 0.57331393
#> 57  0.4711538 0.53093147 0.62138656 0.92046148 0.93004391 0.95130430 0.96183872
#> 58  0.4711538 0.64943097 0.57141374 0.66800038 0.64835800 0.65566321 0.64361788
#> 59  0.4711538 0.42541400 0.43027409 0.30117492 0.36183156 0.29992796 0.28643229
#> 60  0.4711538 0.24537249 0.29963849 0.42931558 0.51048830 0.58927966 0.55978204
#> 61  0.4711538 0.64269314 0.62785202 0.75163561 0.68045267 0.67000184 0.65787202
#> 62  0.4711538 0.51277761 0.60877778 0.75493489 0.66735142 0.63862193 0.60705201
#> 63  0.4711538 0.53377378 0.53228159 0.56245626 0.58414332 0.61176055 0.63813827
#> 64  0.4711538 0.79099666 0.90572246 0.92244949 0.93001276 0.93454809 1.00000000
#> 65  0.4711538 0.73768777 0.61339931 0.72362105 0.70536287 0.69970096 0.69297263
#> 66  0.4711538 0.70767466 0.53408924 0.50675818 0.52181506 0.54559559 0.53663407
#> 67  0.4711538 0.96312042 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#> 68  0.4711538 0.31575995 0.57179559 0.77297374 0.78532935 0.78484987 0.80071187
#> 69  0.4711538 0.69505872 0.78176548 0.74300700 0.72711033 0.70750770 0.74791867
#> 70  0.4711538 0.72276362 0.90232185 0.89364576 0.84428623 0.92659977 0.93952180
#> 71  0.4711538 0.50950893 0.39503961 0.45591683 0.38297596 0.35086204 0.31179970
#> 72  0.4711538 0.14720074 0.13538571 0.00000000 0.00000000 0.02748516 0.08069763
#> 73  0.4711538 0.49275110 0.44937896 0.41856171 0.62470016 0.61654596 0.63914960
#> 74  0.4711538 0.65674324 0.69439259 0.75479685 0.88511667 0.92560996 0.94540783
#> 75  0.4711538 0.68716407 0.57541914 0.59945962 0.54581071 0.55228791 0.56609663
#> 76  0.4711538 0.54839542 0.50508123 0.52627725 0.55765709 0.52543838 0.49807985
#> 77  0.4711538 0.77317727 0.79812663 0.93073165 1.00000000 1.00000000 1.00000000
#> 78  0.4711538 0.85322027 0.76128342 0.81061207 0.85796753 0.87947603 0.88947890
#> 79  0.4711538 0.81659194 0.90228252 0.80744839 0.70383361 0.68468090 0.70672170
#> 80  0.4711538 0.55964651 0.44326524 0.39507689 0.36149039 0.32071350 0.30181332
#> 81  0.4711538 0.87105473 0.86695796 0.89177640 0.74816339 0.69831750 0.69871492
#> 82  0.4711538 0.47715869 0.68930595 0.71280202 0.73606020 0.78321326 0.73754433
#> 83  0.4711538 0.80974821 0.87138779 0.97466313 0.93082943 0.95560886 0.95554583
#> 84  0.4711538 0.67739807 0.85743609 0.98894432 0.96011041 0.90800271 0.92460814
#> 85  0.4711538 0.57131444 0.34250950 0.33855791 0.31118498 0.31383288 0.31932547
#> 86  0.4711538 0.84958765 0.97611051 0.93090902 0.91560248 0.86222031 0.82950069
#> 87  0.4711538 0.57644613 0.41449248 0.48714466 0.54811918 0.57041511 0.56750119
#> 88  0.4711538 0.75932310 0.71214369 0.52234742 0.59011684 0.59023780 0.57484476
#> 89  0.4711538 0.53031516 0.47090892 0.42433053 0.38847912 0.39218094 0.44028895
#> 90  0.4711538 0.76770402 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#> 91  0.4711538 0.38643842 0.37696993 0.44452861 0.49450298 0.46628856 0.50591629
#> 92  0.4711538 0.92591633 1.00000000 0.96084369 0.95688931 0.93400393 0.95078882
#> 93  0.4711538 0.66726042 0.89247800 0.87390628 0.87335977 0.95801535 0.95048831
#> 94  0.4711538 0.32634752 0.41373057 0.48066349 0.67273089 0.62115180 0.63222931
#> 95  0.4711538 0.50472276 0.77159222 0.71730564 0.62350221 0.64335334 0.60173453
#> 96  0.4711538 0.34622269 0.33150717 0.49412629 0.44574013 0.46889514 0.42461082
#> 97  0.4711538 0.55805257 0.50280611 0.58541977 0.52239953 0.53556273 0.55493331
#> 98  0.4711538 0.78090964 0.73429355 0.79385683 0.86651416 0.88151677 0.84875293
#> 99  0.4711538 0.21116352 0.10917861 0.02565398 0.18342015 0.15222876 0.12817884
#> 100 0.4711538 0.66672702 0.78264411 0.86306662 0.75733969 0.77632472 0.76010896
#> 101 0.4711538 0.45317545 0.50149615 0.62617428 0.70904267 0.78134354 0.80270481
#> 102 0.4711538 0.74435376 0.66135006 0.72568147 0.70203564 0.77593538 0.77647730
#> 103 0.4711538 0.34690226 0.56605434 0.52782336 0.50951738 0.46795757 0.50588630
#> 104 0.4711538 0.69496014 0.80515138 0.78871059 0.78008789 0.78042831 0.81071639
#>           [,8]       [,9]      [,10]      [,11]
#> 1   0.62734057 0.65124368 0.65354280 0.65838797
#> 2   0.20266653 0.20646355 0.21046038 0.20470356
#> 3   0.00000000 0.00000000 0.00000000 0.00000000
#> 4   0.60037792 0.60876822 0.60784745 0.60193905
#> 5   0.00000000 0.00000000 0.00000000 0.00000000
#> 6   0.00000000 0.00000000 0.00000000 0.00000000
#> 7   0.13275052 0.14475024 0.13944940 0.14173177
#> 8   0.06373660 0.06326154 0.07570783 0.07715150
#> 9   0.44494367 0.44702123 0.44582742 0.44748513
#> 10  0.40700559 0.41366376 0.41349216 0.41350336
#> 11  0.77762362 0.78373432 0.78262765 0.77949531
#> 12  0.64215450 0.63343656 0.63197252 0.62875264
#> 13  0.18750911 0.18653354 0.17822088 0.17842513
#> 14  0.26050129 0.27908122 0.27196942 0.26996735
#> 15  0.09594988 0.09089468 0.08865039 0.09080522
#> 16  0.00000000 0.00000000 0.00000000 0.00000000
#> 17  0.33619355 0.33546848 0.33501592 0.33670840
#> 18  0.58067856 0.58144968 0.58013547 0.57927640
#> 19  0.24670131 0.24099118 0.23822863 0.23157123
#> 20  0.38607475 0.40041929 0.40843493 0.40860630
#> 21  0.54933618 0.53959876 0.54139679 0.54098845
#> 22  0.26532303 0.26288937 0.27124503 0.27372049
#> 23  0.25426484 0.24411174 0.24044626 0.24169731
#> 24  0.13629899 0.14290857 0.14863967 0.15404679
#> 25  0.65036784 0.64494406 0.64644697 0.64416177
#> 26  0.19040219 0.19924562 0.20152378 0.20582432
#> 27  0.45428317 0.44621001 0.43311399 0.43568840
#> 28  0.11808149 0.11758965 0.11048814 0.11155803
#> 29  0.42401800 0.41725604 0.42362830 0.42160066
#> 30  0.16869228 0.17404171 0.17387211 0.17480429
#> 31  0.43243579 0.42525949 0.43329412 0.43856907
#> 32  0.11698190 0.10877882 0.10636129 0.10689920
#> 33  0.45073883 0.45206217 0.44881973 0.44992624
#> 34  0.09586296 0.09225400 0.09980718 0.10084230
#> 35  0.00000000 0.00000000 0.00000000 0.00000000
#> 36  0.43366789 0.44959439 0.44854107 0.44960414
#> 37  0.31162762 0.31502065 0.31958491 0.32067580
#> 38  0.16099729 0.15122241 0.15079738 0.14817778
#> 39  0.25364311 0.25507184 0.25875102 0.25674198
#> 40  0.53034707 0.52840262 0.52790558 0.52538704
#> 41  0.10186998 0.09060473 0.09445259 0.09150467
#> 42  0.48313687 0.48749107 0.48860245 0.49212211
#> 43  0.57309187 0.54660927 0.54001848 0.54366888
#> 44  0.00000000 0.00000000 0.00000000 0.00000000
#> 45  0.16498860 0.17052210 0.16781405 0.16715155
#> 46  0.64508459 0.65153585 0.65177845 0.65178893
#> 47  0.11991376 0.11705391 0.11873700 0.11591426
#> 48  0.16699846 0.16060597 0.16869180 0.16855238
#> 49  0.61300910 0.60887376 0.60785525 0.61021103
#> 50  0.31497806 0.31197484 0.30517646 0.31174053
#> 51  0.31832188 0.32079774 0.31281691 0.31540371
#> 52  0.00000000 0.00000000 0.00000000 0.00000000
#> 53  0.41513665 0.41604404 0.42408314 0.42075636
#> 54  0.47968839 0.47130262 0.46613648 0.46656474
#> 55  0.21377295 0.18935809 0.18930563 0.19050731
#> 56  0.57455838 0.57977566 0.57922184 0.57800858
#> 57  0.97541446 0.96885561 0.96385549 0.96028958
#> 58  0.63163263 0.64101334 0.63785489 0.64045306
#> 59  0.28914861 0.27879135 0.28784917 0.28782462
#> 60  0.55902384 0.56163014 0.55671889 0.55726966
#> 61  0.65941211 0.65058824 0.64930310 0.65091617
#> 62  0.60397744 0.58982660 0.58812271 0.58793743
#> 63  0.61771766 0.61445477 0.61166047 0.60985403
#> 64  1.00000000 1.00000000 1.00000000 1.00000000
#> 65  0.69769474 0.68981063 0.69214943 0.69038662
#> 66  0.53236916 0.53593471 0.53171176 0.53102146
#> 67  1.00000000 1.00000000 1.00000000 1.00000000
#> 68  0.81246241 0.79688369 0.80360414 0.80689017
#> 69  0.77466902 0.76435744 0.76435231 0.76730276
#> 70  0.92635688 0.91649616 0.91637035 0.92383169
#> 71  0.31958498 0.34299668 0.34334588 0.34349096
#> 72  0.06488309 0.05328264 0.04500807 0.04571864
#> 73  0.65602434 0.64784140 0.64971036 0.64930911
#> 74  0.94553110 0.94743335 0.95548576 0.95086945
#> 75  0.58105153 0.58651055 0.58350518 0.58195688
#> 76  0.50712414 0.50370011 0.50352891 0.50289031
#> 77  1.00000000 1.00000000 1.00000000 1.00000000
#> 78  0.91407469 0.92199353 0.92611779 0.92473772
#> 79  0.68969952 0.69332276 0.69306283 0.69031633
#> 80  0.29656001 0.28290751 0.28672758 0.28334638
#> 81  0.69260333 0.69271104 0.69448149 0.69353087
#> 82  0.75551561 0.75380626 0.76075515 0.76006115
#> 83  0.98111451 0.96789911 0.96789876 0.97167928
#> 84  0.92825575 0.93456912 0.93157477 0.93560641
#> 85  0.32670221 0.33419177 0.33780085 0.33863677
#> 86  0.80658855 0.82379298 0.82202715 0.81981455
#> 87  0.56043949 0.54972085 0.54609905 0.54656224
#> 88  0.58884182 0.58228453 0.57334028 0.57098492
#> 89  0.45579085 0.46221532 0.46149391 0.46635298
#> 90  1.00000000 1.00000000 1.00000000 1.00000000
#> 91  0.52873440 0.54203308 0.53789411 0.54025902
#> 92  0.94412285 0.93390914 0.93429933 0.92818920
#> 93  0.98499381 0.99987093 0.99916201 0.99970721
#> 94  0.61780332 0.60462097 0.60316670 0.59998131
#> 95  0.59883886 0.61566575 0.61062288 0.61003734
#> 96  0.39373512 0.38230743 0.38698306 0.38928062
#> 97  0.54534135 0.54894600 0.55196226 0.55110104
#> 98  0.84908812 0.83344560 0.82896512 0.82882885
#> 99  0.12518468 0.14252417 0.14250387 0.14467219
#> 100 0.74023431 0.74464225 0.75159300 0.75722807
#> 101 0.77862492 0.77776152 0.78102072 0.77766043
#> 102 0.75586058 0.76487586 0.76840548 0.77181270
#> 103 0.48885501 0.49542376 0.50188675 0.50526664
#> 104 0.80418072 0.82546481 0.81861016 0.81635531
modpls.aze$Probs-modpls.aze$Probs.trc
#>     [,1]        [,2]        [,3]         [,4]        [,5]        [,6]
#> 1      0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 2      0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 3      0 -0.09080494 -0.05104846 -0.171669164 -0.21455242 -0.21725391
#> 4      0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 5      0 -0.04408124  0.00000000 -0.071299085 -0.24018962 -0.23445282
#> 6      0 -0.03776963  0.00000000  0.000000000 -0.02597539 -0.06284454
#> 7      0 -0.06930728 -0.19928456 -0.091372612  0.00000000  0.00000000
#> 8      0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 9      0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 10     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 11     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 12     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 13     0  0.00000000 -0.04344681  0.000000000  0.00000000  0.00000000
#> 14     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 15     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 16     0  0.00000000 -0.07276258 -0.051465557 -0.09988241 -0.06790398
#> 17     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 18     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 19     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 20     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 21     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 22     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 23     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 24     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 25     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 26     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 27     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 28     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 29     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 30     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 31     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 32     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 33     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 34     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 35     0  0.00000000  0.00000000 -0.088960742 -0.06232370 -0.08231459
#> 36     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 37     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 38     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 39     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 40     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 41     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 42     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 43     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 44     0  0.00000000  0.00000000  0.000000000 -0.01790809 -0.03785626
#> 45     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 46     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 47     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 48     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 49     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 50     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 51     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 52     0 -0.23250098 -0.28713647 -0.092161735 -0.12709475 -0.18324647
#> 53     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 54     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 55     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 56     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 57     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 58     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 59     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 60     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 61     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 62     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 63     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 64     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 65     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 66     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 67     0  0.00000000  0.17012215  0.081167947  0.22497425  0.21728258
#> 68     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 69     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 70     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 71     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 72     0  0.00000000  0.00000000 -0.044738287 -0.05529233  0.00000000
#> 73     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 74     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 75     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 76     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 77     0  0.00000000  0.00000000  0.000000000  0.10301473  0.08723742
#> 78     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 79     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 80     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 81     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 82     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 83     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 84     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 85     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 86     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 87     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 88     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 89     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 90     0  0.00000000  0.07649866  0.008644293  0.06363018  0.09017457
#> 91     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 92     0  0.00000000  0.03707888  0.000000000  0.00000000  0.00000000
#> 93     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 94     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 95     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 96     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 97     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 98     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 99     0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 100    0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 101    0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 102    0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 103    0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#> 104    0  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000
#>             [,7]         [,8]         [,9]       [,10]        [,11]
#> 1    0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 2    0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 3   -0.224440890 -0.193652144 -0.201652437 -0.20167289 -0.200815325
#> 4    0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 5   -0.263270486 -0.311781941 -0.310765976 -0.31066830 -0.314013823
#> 6   -0.103410961 -0.076840858 -0.080016598 -0.08436785 -0.087771351
#> 7    0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 8    0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 9    0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 10   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 11   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 12   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 13   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 14   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 15   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 16  -0.062392124 -0.039505632 -0.023469898 -0.02045198 -0.027151896
#> 17   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 18   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 19   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 20   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 21   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 22   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 23   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 24   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 25   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 26   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 27   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 28   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 29   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 30   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 31   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 32   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 33   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 34   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 35  -0.083905063 -0.085155814 -0.086690885 -0.08262604 -0.082725745
#> 36   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 37   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 38   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 39   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 40   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 41   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 42   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 43   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 44  -0.014514870 -0.004369241 -0.002647953 -0.00844687 -0.010738872
#> 45   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 46   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 47   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 48   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 49   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 50   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 51   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 52  -0.195842587 -0.184768370 -0.177691131 -0.18225253 -0.181720401
#> 53   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 54   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 55   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 56   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 57   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 58   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 59   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 60   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 61   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 62   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 63   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 64   0.001662791  0.008968255  0.008852814  0.01022555  0.006667285
#> 65   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 66   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 67   0.221111495  0.202569147  0.194030443  0.19310867  0.192281835
#> 68   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 69   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 70   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 71   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 72   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 73   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 74   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 75   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 76   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 77   0.064638063  0.090636226  0.099163425  0.09692704  0.096849873
#> 78   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 79   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 80   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 81   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 82   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 83   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 84   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 85   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 86   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 87   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 88   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 89   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 90   0.087035865  0.113909884  0.123382027  0.12669615  0.125596908
#> 91   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 92   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 93   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 94   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 95   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 96   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 97   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 98   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 99   0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 100  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 101  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 102  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 103  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> 104  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000

#Repeated cross validation of the model (NK=100 times)
cv.modpls.aze<-cv.plsR(y~.,data=aze_compl,10,NK=100, verbose=FALSE)
res.cv.modpls.aze<-cvtable(summary(cv.modpls.aze,MClassed=TRUE))
#> ____************************************************____
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
#> 
#> NK: 1,  2,  3,  4,  5,  6,  7,  8,  9,  10
#> NK: 11,  12,  13,  14,  15,  16,  17,  18,  19,  20
#> NK: 21,  22,  23,  24,  25,  26,  27,  28,  29,  30
#> NK: 31,  32,  33,  34,  35,  36,  37,  38,  39,  40
#> NK: 41,  42,  43,  44,  45,  46,  47,  48,  49,  50
#> NK: 51,  52,  53,  54,  55,  56,  57,  58,  59,  60
#> NK: 61,  62,  63,  64,  65,  66,  67,  68,  69,  70
#> NK: 71,  72,  73,  74,  75,  76,  77,  78,  79,  80
#> NK: 81,  82,  83,  84,  85,  86,  87,  88,  89,  90
#> NK: 91,  92,  93,  94,  95,  96,  97,  98,  99,  100
#> 
#> CV MissClassed criterion:
#>  1  2  3  4  5  6  7  8  9 10 
#> 24 10 25 11  5  5  7  8  2  3 
#> 
#> CV Q2 criterion:
#>   0 
#> 100 
#> 
#> CV Press criterion:
#>  1  2 
#> 81 19 
#High discrepancy in the number of component choice using repeated cross validation
#and missclassed criterion
plot(res.cv.modpls.aze)


rm(list=c("Xaze_compl","yaze_compl","modpls.aze","cv.modpls.aze","res.cv.modpls.aze"))

#24 predictors
dimX <- 24
#2 components
Astar <- 2
simul_data_UniYX(dimX,Astar)
#>          Y         X1         X2         X3         X4         X5         X6 
#> 11.6445768  0.6175008  0.6056228  5.1301635  0.6107449  0.6138664  5.1302654 
#>         X7         X8         X9        X10        X11        X12        X13 
#>  0.6402793  0.6103912  5.1466522  0.6048168  0.6205430  5.1257361  0.6101874 
#>        X14        X15        X16        X17        X18        X19        X20 
#>  0.6103582  5.1152508  0.6286553  0.6139872  5.1205242  0.6145876  0.6165494 
#>        X21        X22        X23        X24 
#>  5.1176202  0.6228756  0.6110989  5.1274519 
dataAstar2 <- data.frame(t(replicate(250,simul_data_UniYX(dimX,Astar))))
modpls.A2<- plsR(Y~.,data=dataAstar2,10,typeVC="standard")
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
modpls.A2
#> Number of required components:
#> [1] 10
#> Number of successfully computed components:
#> [1] 10
#> Coefficients:
#>                   [,1]
#> Intercept -0.009308745
#> X1        -1.252039537
#> X2         0.251769360
#> X3         1.069303144
#> X4        -1.523206574
#> X5         0.429339251
#> X6        -0.103727035
#> X7         0.514838746
#> X8        -0.008132813
#> X9        -1.109769998
#> X10        0.697707208
#> X11        0.110179634
#> X12        0.263432702
#> X13        0.048981023
#> X14        1.102957195
#> X15        2.998639072
#> X16       -1.030175931
#> X17        0.005090035
#> X18       -0.108388491
#> X19        1.071595090
#> X20        0.619985470
#> X21       -0.996526239
#> X22       -0.639107098
#> X23        0.152756799
#> X24        0.266087715
#> Leave one out cross validated PRESS, Information criteria and Fit statistics:
#>                  AIC   Q2cum_Y LimQ2_Y        Q2_Y     PRESS_Y        RSS_Y
#> Nb_Comp_0  1707.6611        NA      NA          NA          NA 13336.080111
#> Nb_Comp_1  1268.8401 0.8246912  0.0975  0.82469120 2337.932230  2286.881298
#> Nb_Comp_2  -202.5201 0.9995127  0.0975  0.99722014    6.357207     6.306012
#> Nb_Comp_3  -219.4997 0.9994751  0.0975 -0.07700506    6.791607     5.844991
#> Nb_Comp_4  -219.7889 0.9993831  0.0975 -0.17536483    6.869997     5.791713
#> Nb_Comp_5  -218.0111 0.9992724  0.0975 -0.17939626    6.830725     5.786567
#> Nb_Comp_6  -216.0264 0.9991613  0.0975 -0.15268975    6.670117     5.786214
#> Nb_Comp_7  -214.0278 0.9990390  0.0975 -0.14581108    6.629908     5.786181
#> Nb_Comp_8  -212.0280 0.9989079  0.0975 -0.13650388    6.576018     5.786178
#> Nb_Comp_9  -210.0280 0.9987662  0.0975 -0.12969385    6.536609     5.786177
#> Nb_Comp_10 -208.0280 0.9986200  0.0975 -0.11849061    6.471785     5.786177
#>                 R2_Y R2_residY  RSS_residY PRESS_residY   Q2_residY  LimQ2
#> Nb_Comp_0         NA        NA 249.0000000           NA          NA     NA
#> Nb_Comp_1  0.8285192 0.8285192  42.6987119   43.6518917  0.82469120 0.0975
#> Nb_Comp_2  0.9995271 0.9995271   0.1177405    0.1186964  0.99722014 0.0975
#> Nb_Comp_3  0.9995617 0.9995617   0.1091327    0.1268071 -0.07700506 0.0975
#> Nb_Comp_4  0.9995657 0.9995657   0.1081380    0.1282708 -0.17536483 0.0975
#> Nb_Comp_5  0.9995661 0.9995661   0.1080419    0.1275375 -0.17939626 0.0975
#> Nb_Comp_6  0.9995661 0.9995661   0.1080353    0.1245388 -0.15268975 0.0975
#> Nb_Comp_7  0.9995661 0.9995661   0.1080347    0.1237880 -0.14581108 0.0975
#> Nb_Comp_8  0.9995661 0.9995661   0.1080346    0.1227818 -0.13650388 0.0975
#> Nb_Comp_9  0.9995661 0.9995661   0.1080346    0.1220460 -0.12969385 0.0975
#> Nb_Comp_10 0.9995661 0.9995661   0.1080346    0.1208357 -0.11849061 0.0975
#>            Q2cum_residY    AIC.std   DoF.dof sigmahat.dof     AIC.dof
#> Nb_Comp_0            NA   712.4673  1.000000    7.3183710 53.77278888
#> Nb_Comp_1     0.8246912   273.6462  2.555426    3.0339405  9.33570257
#> Nb_Comp_2     0.9995127 -1197.7140  3.000068    0.1594599  0.02583432
#> Nb_Comp_3     0.9994751 -1214.6936 22.574476    0.1599630  0.02800108
#> Nb_Comp_4     0.9993831 -1214.9828 25.000000    0.1600845  0.02829226
#> Nb_Comp_5     0.9992724 -1213.2050 25.000000    0.1600134  0.02826713
#> Nb_Comp_6     0.9991613 -1211.2203 25.000000    0.1600085  0.02826540
#> Nb_Comp_7     0.9990390 -1209.2217 25.000000    0.1600080  0.02826524
#> Nb_Comp_8     0.9989079 -1207.2219 25.000000    0.1600080  0.02826522
#> Nb_Comp_9     0.9987662 -1205.2219 25.000000    0.1600080  0.02826522
#> Nb_Comp_10    0.9986200 -1203.2219 25.000000    0.1600080  0.02826522
#>                BIC.dof  GMDL.dof DoF.naive sigmahat.naive   AIC.naive
#> Nb_Comp_0  54.52720631  499.7901         1      7.3183710 53.77278888
#> Nb_Comp_1   9.66703221  288.0890         2      3.0366586  9.29506592
#> Nb_Comp_2   0.02690885 -438.1210         3      0.1597824  0.02583678
#> Nb_Comp_3   0.03613762 -342.0295         4      0.1541432  0.02414029
#> Nb_Comp_4   0.03731673 -330.9487         5      0.1537519  0.02411244
#> Nb_Comp_5   0.03728357 -331.0487         6      0.1539982  0.02428461
#> Nb_Comp_6   0.03728129 -331.0556         7      0.1543100  0.02447830
#> Nb_Comp_7   0.03728109 -331.0562         8      0.1546281  0.02467496
#> Nb_Comp_8   0.03728106 -331.0563         9      0.1549485  0.02487336
#> Nb_Comp_9   0.03728106 -331.0563        10      0.1552710  0.02507344
#> Nb_Comp_10  0.03728106 -331.0563        11      0.1555955  0.02527518
#>              BIC.naive GMDL.naive
#> Nb_Comp_0  54.52720631   499.7901
#> Nb_Comp_1   9.55484538   286.8472
#> Nb_Comp_2   0.02691563  -437.6224
#> Nb_Comp_3   0.02547901  -441.0025
#> Nb_Comp_4   0.02577736  -436.2568
#> Nb_Comp_5   0.02628892  -430.5960
#> Nb_Comp_6   0.02682615  -424.9195
#> Nb_Comp_7   0.02736928  -419.3100
#> Nb_Comp_8   0.02791705  -413.7646
#> Nb_Comp_9   0.02846940  -408.2768
#> Nb_Comp_10  0.02902638  -402.8412
cv.modpls.A2<-cv.plsR(Y~.,data=dataAstar2,10,NK=100, verbose=FALSE)
res.cv.modpls.A2<-cvtable(summary(cv.modpls.A2,verbose=FALSE))
#> Error in eval(mf, parent.frame()): object 'dataAstar2' not found
#Perfect choice for the Q2 criterion in PLSR
plot(res.cv.modpls.A2)
#> Error: object 'res.cv.modpls.A2' not found

#Binarized data.frame
simbin1 <- data.frame(dicho(dataAstar2))
modpls.B2 <- plsR(Y~.,data=simbin1,10,typeVC="standard",MClassed=TRUE, verbose=FALSE)
modpls.B2
#> Number of required components:
#> [1] 10
#> Number of successfully computed components:
#> [1] 8
#> Coefficients:
#>                    [,1]
#> Intercept -2.232996e-02
#> X1        -2.564600e-04
#> X2         1.827582e-02
#> X3         3.200076e-01
#> X4         1.827582e-02
#> X5         8.548667e-05
#> X6        -7.441941e-02
#> X7         1.827582e-02
#> X8         8.548667e-05
#> X9         1.600038e-01
#> X10       -6.653264e-03
#> X11        1.794648e-02
#> X12        5.540587e-01
#> X13        1.794648e-02
#> X14       -6.307629e-03
#> X15       -7.441941e-02
#> X16        1.827582e-02
#> X17        1.827582e-02
#> X18       -7.441941e-02
#> X19       -1.961416e-02
#> X20       -6.653264e-03
#> X21       -7.441941e-02
#> X22        8.548667e-05
#> X23        1.827582e-02
#> X24        1.600038e-01
#> Leave one out cross validated PRESS, Information criteria and Fit statistics:
#>                 AIC   Q2cum_Y LimQ2_Y          Q2_Y   PRESS_Y     RSS_Y
#> Nb_Comp_0 366.87968        NA      NA            NA        NA 62.496000
#> Nb_Comp_1  26.83779 0.7381855  0.0975  0.7381855034 16.362359 15.909796
#> Nb_Comp_2 -76.52049 0.8252118  0.0975  0.3323966613 10.621433 10.438510
#> Nb_Comp_3 -91.51975 0.8283713  0.0975  0.0180764406 10.249819  9.752317
#> Nb_Comp_4 -89.75294 0.8284216  0.0975  0.0002930630  9.749459  9.743224
#> Nb_Comp_5 -87.81866 0.8284417  0.0975  0.0001171978  9.742083  9.740663
#> Nb_Comp_6 -85.82121 0.8284216  0.0975 -0.0001170090  9.741803  9.740564
#> Nb_Comp_7 -83.82127 0.8283970  0.0975 -0.0001436374  9.741963  9.740562
#> Nb_Comp_8 -81.82127 0.8283740  0.0975 -0.0001343106  9.741870  9.740562
#>                R2_Y MissClassed R2_residY RSS_residY PRESS_residY     Q2_residY
#> Nb_Comp_0        NA         124        NA  249.00000           NA            NA
#> Nb_Comp_1 0.7454270          11 0.7454270   63.38868     65.19181  0.7381855034
#> Nb_Comp_2 0.8329731          13 0.8329731   41.58969     42.31850  0.3323966613
#> Nb_Comp_3 0.8439529          11 0.8439529   38.85572     40.83789  0.0180764406
#> Nb_Comp_4 0.8440984          11 0.8440984   38.81949     38.84433  0.0002930630
#> Nb_Comp_5 0.8441394          11 0.8441394   38.80929     38.81494  0.0001171978
#> Nb_Comp_6 0.8441410          11 0.8441410   38.80889     38.81383 -0.0001170090
#> Nb_Comp_7 0.8441410          11 0.8441410   38.80888     38.81447 -0.0001436374
#> Nb_Comp_8 0.8441410          11 0.8441410   38.80888     38.81409 -0.0001343106
#>            LimQ2 Q2cum_residY  AIC.std   DoF.dof sigmahat.dof    AIC.dof
#> Nb_Comp_0     NA           NA 712.4673  1.000000    0.5009870 0.25199190
#> Nb_Comp_1 0.0975    0.7381855 372.4254  3.112471    0.2533407 0.06523729
#> Nb_Comp_2 0.0975    0.8252118 269.0671  3.041696    0.2051776 0.04277843
#> Nb_Comp_3 0.0975    0.8283713 254.0678 10.053091    0.2011839 0.04226445
#> Nb_Comp_4 0.0975    0.8284216 255.8346  8.294626    0.2003603 0.04163675
#> Nb_Comp_5 0.0975    0.8284417 257.7689  8.533333    0.2004325 0.04170514
#> Nb_Comp_6 0.0975    0.8284216 259.7664  8.819577    0.2005499 0.04180006
#> Nb_Comp_7 0.0975    0.8283970 261.7663  8.983400    0.2006178 0.04185472
#> Nb_Comp_8 0.0975    0.8283740 263.7663  8.999573    0.2006245 0.04186012
#>              BIC.dof  GMDL.dof DoF.naive sigmahat.naive  AIC.naive  BIC.naive
#> Nb_Comp_0 0.25552728 -167.2823         1      0.5009870 0.25199190 0.25552728
#> Nb_Comp_1 0.06805112 -330.7000         2      0.2532832 0.06466562 0.06647290
#> Nb_Comp_2 0.04458211 -382.8861         3      0.2055752 0.04276831 0.04455416
#> Nb_Comp_3 0.04799596 -369.7858         4      0.1991069 0.04027786 0.04251151
#> Nb_Comp_4 0.04632708 -374.9363         5      0.1994198 0.04056363 0.04336448
#> Nb_Comp_5 0.04653393 -374.2734         6      0.1998018 0.04087885 0.04425275
#> Nb_Comp_6 0.04679668 -373.4463         7      0.2002115 0.04120700 0.04515938
#> Nb_Comp_7 0.04694759 -372.9744         8      0.2006247 0.04153826 0.04607393
#> Nb_Comp_8 0.04696250 -372.9279         9      0.2010405 0.04187229 0.04699609
#>           GMDL.naive
#> Nb_Comp_0  -167.2823
#> Nb_Comp_1  -333.8147
#> Nb_Comp_2  -382.5286
#> Nb_Comp_3  -387.5578
#> Nb_Comp_4  -384.4408
#> Nb_Comp_5  -381.3439
#> Nb_Comp_6  -378.3019
#> Nb_Comp_7  -375.3325
#> Nb_Comp_8  -372.4277
modpls.B2$Probs
#>      [,1]        [,2]         [,3]          [,4]        [,5]         [,6]
#> 1   0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 2   0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 3   0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 4   0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 5   0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 6   0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 7   0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 8   0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 9   0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 10  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 11  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 12  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 13  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 14  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 15  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 16  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 17  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 18  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 19  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 20  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 21  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 22  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 23  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 24  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 25  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 26  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 27  0.496  0.03606045  0.002955833  0.0001432602 -0.01599899  0.007925162
#> 28  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 29  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 30  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 31  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 32  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 33  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 34  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 35  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 36  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 37  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 38  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 39  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 40  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 41  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 42  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 43  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 44  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 45  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 46  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 47  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 48  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 49  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 50  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 51  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 52  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 53  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 54  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 55  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 56  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 57  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 58  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 59  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 60  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 61  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 62  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 63  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 64  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 65  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 66  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 67  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 68  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 69  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 70  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 71  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 72  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 73  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 74  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 75  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 76  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 77  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 78  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 79  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 80  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 81  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 82  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 83  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 84  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 85  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 86  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 87  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 88  0.496  1.00500492  0.973254333  1.0029362526  0.99021061  0.997293320
#> 89  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 90  0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 91  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 92  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 93  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 94  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 95  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 96  0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 97  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 98  0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 99  0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 100 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 101 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 102 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 103 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 104 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 105 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 106 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 107 0.496  0.90947840  0.945577585  1.0369797504  0.97619695  0.995322751
#> 108 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 109 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 110 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 111 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 112 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 113 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 114 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 115 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 116 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 117 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 118 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 119 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 120 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 121 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 122 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 123 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 124 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 125 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 126 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 127 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 128 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 129 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 130 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 131 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 132 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 133 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 134 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 135 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 136 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 137 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 138 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 139 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 140 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 141 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 142 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 143 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 144 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 145 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 146 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 147 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 148 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 149 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 150 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 151 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 152 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 153 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 154 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 155 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 156 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 157 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 158 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 159 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 160 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 161 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 162 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 163 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 164 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 165 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 166 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 167 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 168 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 169 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 170 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 171 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 172 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 173 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 174 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 175 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 176 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 177 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 178 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 179 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 180 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 181 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 182 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 183 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 184 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 185 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 186 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 187 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 188 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 189 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 190 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 191 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 192 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 193 0.496  0.45907245  0.636656724  0.0755600483  0.02597481  0.001514519
#> 194 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 195 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 196 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 197 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 198 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 199 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 200 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 201 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 202 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 203 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 204 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 205 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 206 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 207 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 208 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 209 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 210 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 211 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 212 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 213 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 214 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 215 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 216 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 217 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 218 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 219 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 220 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 221 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 222 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 223 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 224 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 225 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 226 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 227 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 228 0.496  0.95759591  0.959727077  0.9862193637  1.01721602  0.999177807
#> 229 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 230 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 231 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 232 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 233 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 234 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 235 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 236 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 237 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 238 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 239 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 240 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 241 0.496  0.36839434  0.524751547 -0.0668244029 -0.02771783 -0.002376501
#> 242 0.496  0.64117170  0.862240206  0.8736969788  0.87418616  0.874112455
#> 243 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 244 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#> 245 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 246 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 247 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 248 0.496  0.30502486  0.086944139  0.0840816281  0.08415168  0.083947110
#> 249 0.496  1.02826346  0.979507152  0.9800896768  0.98049332  0.980469640
#> 250 0.496 -0.08206690 -0.030322806 -0.0223110698 -0.02215549 -0.022410075
#>              [,7]          [,8]          [,9]
#> 1    8.740665e-01  8.740690e-01  8.740663e-01
#> 2    9.803926e-01  9.803899e-01  9.803858e-01
#> 3    9.803926e-01  9.803899e-01  9.803858e-01
#> 4    8.740665e-01  8.740690e-01  8.740663e-01
#> 5    9.803926e-01  9.803899e-01  9.803858e-01
#> 6    9.803926e-01  9.803899e-01  9.803858e-01
#> 7    8.740665e-01  8.740690e-01  8.740663e-01
#> 8   -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 9   -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 10   9.803926e-01  9.803899e-01  9.803858e-01
#> 11  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 12   9.803926e-01  9.803899e-01  9.803858e-01
#> 13  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 14   8.740665e-01  8.740690e-01  8.740663e-01
#> 15  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 16   8.399153e-02  8.399119e-02  8.398958e-02
#> 17   9.803926e-01  9.803899e-01  9.803858e-01
#> 18  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 19   8.399153e-02  8.399119e-02  8.398958e-02
#> 20   8.399153e-02  8.399119e-02  8.398958e-02
#> 21   9.803926e-01  9.803899e-01  9.803858e-01
#> 22   8.399153e-02  8.399119e-02  8.398958e-02
#> 23   8.399153e-02  8.399119e-02  8.398958e-02
#> 24  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 25   9.803926e-01  9.803899e-01  9.803858e-01
#> 26   9.803926e-01  9.803899e-01  9.803858e-01
#> 27   1.814101e-04 -9.424424e-05 -1.054712e-15
#> 28   8.740665e-01  8.740690e-01  8.740663e-01
#> 29   9.803926e-01  9.803899e-01  9.803858e-01
#> 30  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 31   9.803926e-01  9.803899e-01  9.803858e-01
#> 32   9.803926e-01  9.803899e-01  9.803858e-01
#> 33   8.740665e-01  8.740690e-01  8.740663e-01
#> 34  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 35   8.740665e-01  8.740690e-01  8.740663e-01
#> 36   9.803926e-01  9.803899e-01  9.803858e-01
#> 37   9.803926e-01  9.803899e-01  9.803858e-01
#> 38  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 39  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 40  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 41   8.740665e-01  8.740690e-01  8.740663e-01
#> 42   8.399153e-02  8.399119e-02  8.398958e-02
#> 43   8.740665e-01  8.740690e-01  8.740663e-01
#> 44   8.399153e-02  8.399119e-02  8.398958e-02
#> 45   8.740665e-01  8.740690e-01  8.740663e-01
#> 46  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 47   9.803926e-01  9.803899e-01  9.803858e-01
#> 48   8.740665e-01  8.740690e-01  8.740663e-01
#> 49   8.740665e-01  8.740690e-01  8.740663e-01
#> 50   8.399153e-02  8.399119e-02  8.398958e-02
#> 51   8.399153e-02  8.399119e-02  8.398958e-02
#> 52  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 53   9.803926e-01  9.803899e-01  9.803858e-01
#> 54   8.740665e-01  8.740690e-01  8.740663e-01
#> 55   9.803926e-01  9.803899e-01  9.803858e-01
#> 56   9.803926e-01  9.803899e-01  9.803858e-01
#> 57   9.803926e-01  9.803899e-01  9.803858e-01
#> 58   8.399153e-02  8.399119e-02  8.398958e-02
#> 59  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 60   9.803926e-01  9.803899e-01  9.803858e-01
#> 61   8.740665e-01  8.740690e-01  8.740663e-01
#> 62   8.740665e-01  8.740690e-01  8.740663e-01
#> 63  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 64  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 65   9.803926e-01  9.803899e-01  9.803858e-01
#> 66   9.803926e-01  9.803899e-01  9.803858e-01
#> 67   9.803926e-01  9.803899e-01  9.803858e-01
#> 68  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 69   8.399153e-02  8.399119e-02  8.398958e-02
#> 70  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 71   8.399153e-02  8.399119e-02  8.398958e-02
#> 72   8.399153e-02  8.399119e-02  8.398958e-02
#> 73  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 74   8.399153e-02  8.399119e-02  8.398958e-02
#> 75  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 76   8.740665e-01  8.740690e-01  8.740663e-01
#> 77   8.399153e-02  8.399119e-02  8.398958e-02
#> 78   8.399153e-02  8.399119e-02  8.398958e-02
#> 79   9.803926e-01  9.803899e-01  9.803858e-01
#> 80   8.740665e-01  8.740690e-01  8.740663e-01
#> 81   8.740665e-01  8.740690e-01  8.740663e-01
#> 82   8.740665e-01  8.740690e-01  8.740663e-01
#> 83   8.399153e-02  8.399119e-02  8.398958e-02
#> 84   8.740665e-01  8.740690e-01  8.740663e-01
#> 85   8.399153e-02  8.399119e-02  8.398958e-02
#> 86   8.399153e-02  8.399119e-02  8.398958e-02
#> 87   9.803926e-01  9.803899e-01  9.803858e-01
#> 88   9.986012e-01  9.998504e-01  1.000000e+00
#> 89   9.803926e-01  9.803899e-01  9.803858e-01
#> 90   8.740665e-01  8.740690e-01  8.740663e-01
#> 91   8.399153e-02  8.399119e-02  8.398958e-02
#> 92  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 93  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 94  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 95   9.803926e-01  9.803899e-01  9.803858e-01
#> 96   8.399153e-02  8.399119e-02  8.398958e-02
#> 97  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 98  -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 99   9.803926e-01  9.803899e-01  9.803858e-01
#> 100  8.740665e-01  8.740690e-01  8.740663e-01
#> 101  9.803926e-01  9.803899e-01  9.803858e-01
#> 102  8.399153e-02  8.399119e-02  8.398958e-02
#> 103  9.803926e-01  9.803899e-01  9.803858e-01
#> 104 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 105  9.803926e-01  9.803899e-01  9.803858e-01
#> 106 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 107  1.000509e+00  9.999525e-01  1.000000e+00
#> 108  8.740665e-01  8.740690e-01  8.740663e-01
#> 109  8.399153e-02  8.399119e-02  8.398958e-02
#> 110 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 111  9.803926e-01  9.803899e-01  9.803858e-01
#> 112  9.803926e-01  9.803899e-01  9.803858e-01
#> 113  8.740665e-01  8.740690e-01  8.740663e-01
#> 114  8.399153e-02  8.399119e-02  8.398958e-02
#> 115  9.803926e-01  9.803899e-01  9.803858e-01
#> 116  8.740665e-01  8.740690e-01  8.740663e-01
#> 117  8.740665e-01  8.740690e-01  8.740663e-01
#> 118 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 119  8.399153e-02  8.399119e-02  8.398958e-02
#> 120  9.803926e-01  9.803899e-01  9.803858e-01
#> 121 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 122  8.740665e-01  8.740690e-01  8.740663e-01
#> 123 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 124  8.399153e-02  8.399119e-02  8.398958e-02
#> 125  8.399153e-02  8.399119e-02  8.398958e-02
#> 126  9.803926e-01  9.803899e-01  9.803858e-01
#> 127  9.803926e-01  9.803899e-01  9.803858e-01
#> 128  9.803926e-01  9.803899e-01  9.803858e-01
#> 129  9.803926e-01  9.803899e-01  9.803858e-01
#> 130  8.399153e-02  8.399119e-02  8.398958e-02
#> 131  9.803926e-01  9.803899e-01  9.803858e-01
#> 132  8.399153e-02  8.399119e-02  8.398958e-02
#> 133  9.803926e-01  9.803899e-01  9.803858e-01
#> 134  8.740665e-01  8.740690e-01  8.740663e-01
#> 135  9.803926e-01  9.803899e-01  9.803858e-01
#> 136  8.740665e-01  8.740690e-01  8.740663e-01
#> 137  9.803926e-01  9.803899e-01  9.803858e-01
#> 138  8.740665e-01  8.740690e-01  8.740663e-01
#> 139  8.399153e-02  8.399119e-02  8.398958e-02
#> 140 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 141 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 142 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 143  8.740665e-01  8.740690e-01  8.740663e-01
#> 144  8.740665e-01  8.740690e-01  8.740663e-01
#> 145  8.740665e-01  8.740690e-01  8.740663e-01
#> 146  8.399153e-02  8.399119e-02  8.398958e-02
#> 147  8.740665e-01  8.740690e-01  8.740663e-01
#> 148  9.803926e-01  9.803899e-01  9.803858e-01
#> 149 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 150  8.740665e-01  8.740690e-01  8.740663e-01
#> 151 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 152 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 153 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 154  9.803926e-01  9.803899e-01  9.803858e-01
#> 155  9.803926e-01  9.803899e-01  9.803858e-01
#> 156  9.803926e-01  9.803899e-01  9.803858e-01
#> 157 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 158  8.399153e-02  8.399119e-02  8.398958e-02
#> 159 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 160 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 161  8.399153e-02  8.399119e-02  8.398958e-02
#> 162 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 163 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 164  9.803926e-01  9.803899e-01  9.803858e-01
#> 165  8.399153e-02  8.399119e-02  8.398958e-02
#> 166  8.740665e-01  8.740690e-01  8.740663e-01
#> 167 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 168  8.740665e-01  8.740690e-01  8.740663e-01
#> 169  8.399153e-02  8.399119e-02  8.398958e-02
#> 170  9.803926e-01  9.803899e-01  9.803858e-01
#> 171  8.399153e-02  8.399119e-02  8.398958e-02
#> 172  8.740665e-01  8.740690e-01  8.740663e-01
#> 173  8.399153e-02  8.399119e-02  8.398958e-02
#> 174  9.803926e-01  9.803899e-01  9.803858e-01
#> 175  9.803926e-01  9.803899e-01  9.803858e-01
#> 176 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 177  9.803926e-01  9.803899e-01  9.803858e-01
#> 178  8.740665e-01  8.740690e-01  8.740663e-01
#> 179 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 180 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 181  8.399153e-02  8.399119e-02  8.398958e-02
#> 182  9.803926e-01  9.803899e-01  9.803858e-01
#> 183  8.399153e-02  8.399119e-02  8.398958e-02
#> 184  9.803926e-01  9.803899e-01  9.803858e-01
#> 185 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 186  9.803926e-01  9.803899e-01  9.803858e-01
#> 187  8.740665e-01  8.740690e-01  8.740663e-01
#> 188  8.740665e-01  8.740690e-01  8.740663e-01
#> 189  9.803926e-01  9.803899e-01  9.803858e-01
#> 190 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 191  9.803926e-01  9.803899e-01  9.803858e-01
#> 192 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 193 -4.099320e-05 -1.534115e-05  9.436896e-16
#> 194  8.740665e-01  8.740690e-01  8.740663e-01
#> 195  8.399153e-02  8.399119e-02  8.398958e-02
#> 196  9.803926e-01  9.803899e-01  9.803858e-01
#> 197  8.740665e-01  8.740690e-01  8.740663e-01
#> 198  9.803926e-01  9.803899e-01  9.803858e-01
#> 199  8.740665e-01  8.740690e-01  8.740663e-01
#> 200 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 201  8.740665e-01  8.740690e-01  8.740663e-01
#> 202  8.399153e-02  8.399119e-02  8.398958e-02
#> 203  8.740665e-01  8.740690e-01  8.740663e-01
#> 204  8.740665e-01  8.740690e-01  8.740663e-01
#> 205  8.399153e-02  8.399119e-02  8.398958e-02
#> 206 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 207  9.803926e-01  9.803899e-01  9.803858e-01
#> 208  8.399153e-02  8.399119e-02  8.398958e-02
#> 209 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 210  9.803926e-01  9.803899e-01  9.803858e-01
#> 211  8.399153e-02  8.399119e-02  8.398958e-02
#> 212 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 213  9.803926e-01  9.803899e-01  9.803858e-01
#> 214  9.803926e-01  9.803899e-01  9.803858e-01
#> 215  9.803926e-01  9.803899e-01  9.803858e-01
#> 216  8.399153e-02  8.399119e-02  8.398958e-02
#> 217  8.399153e-02  8.399119e-02  8.398958e-02
#> 218  9.803926e-01  9.803899e-01  9.803858e-01
#> 219 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 220  8.399153e-02  8.399119e-02  8.398958e-02
#> 221 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 222  8.399153e-02  8.399119e-02  8.398958e-02
#> 223  8.740665e-01  8.740690e-01  8.740663e-01
#> 224 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 225  8.740665e-01  8.740690e-01  8.740663e-01
#> 226  8.740665e-01  8.740690e-01  8.740663e-01
#> 227 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 228  1.000372e+00  9.997744e-01  1.000000e+00
#> 229  9.803926e-01  9.803899e-01  9.803858e-01
#> 230  9.803926e-01  9.803899e-01  9.803858e-01
#> 231 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 232  9.803926e-01  9.803899e-01  9.803858e-01
#> 233 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 234  9.803926e-01  9.803899e-01  9.803858e-01
#> 235 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 236 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 237  8.399153e-02  8.399119e-02  8.398958e-02
#> 238  9.803926e-01  9.803899e-01  9.803858e-01
#> 239  9.803926e-01  9.803899e-01  9.803858e-01
#> 240  8.740665e-01  8.740690e-01  8.740663e-01
#> 241  5.931236e-05 -7.182277e-06  7.216450e-16
#> 242  8.740665e-01  8.740690e-01  8.740663e-01
#> 243  9.803926e-01  9.803899e-01  9.803858e-01
#> 244 -2.233453e-02 -2.232974e-02 -2.232996e-02
#> 245  8.399153e-02  8.399119e-02  8.398958e-02
#> 246  8.399153e-02  8.399119e-02  8.398958e-02
#> 247  8.399153e-02  8.399119e-02  8.398958e-02
#> 248  8.399153e-02  8.399119e-02  8.398958e-02
#> 249  9.803926e-01  9.803899e-01  9.803858e-01
#> 250 -2.233453e-02 -2.232974e-02 -2.232996e-02
modpls.B2$Probs.trc
#>      [,1]       [,2]        [,3]         [,4]       [,5]        [,6]
#> 1   0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 2   0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 3   0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 4   0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 5   0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 6   0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 7   0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 8   0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 9   0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 10  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 11  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 12  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 13  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 14  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 15  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 16  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 17  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 18  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 19  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 20  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 21  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 22  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 23  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 24  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 25  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 26  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 27  0.496 0.03606045 0.002955833 0.0001432602 0.00000000 0.007925162
#> 28  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 29  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 30  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 31  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 32  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 33  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 34  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 35  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 36  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 37  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 38  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 39  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 40  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 41  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 42  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 43  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 44  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 45  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 46  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 47  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 48  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 49  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 50  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 51  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 52  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 53  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 54  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 55  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 56  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 57  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 58  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 59  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 60  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 61  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 62  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 63  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 64  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 65  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 66  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 67  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 68  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 69  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 70  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 71  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 72  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 73  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 74  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 75  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 76  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 77  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 78  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 79  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 80  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 81  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 82  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 83  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 84  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 85  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 86  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 87  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 88  0.496 1.00000000 0.973254333 1.0000000000 0.99021061 0.997293320
#> 89  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 90  0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 91  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 92  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 93  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 94  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 95  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 96  0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 97  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 98  0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 99  0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 100 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 101 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 102 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 103 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 104 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 105 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 106 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 107 0.496 0.90947840 0.945577585 1.0000000000 0.97619695 0.995322751
#> 108 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 109 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 110 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 111 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 112 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 113 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 114 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 115 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 116 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 117 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 118 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 119 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 120 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 121 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 122 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 123 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 124 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 125 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 126 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 127 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 128 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 129 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 130 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 131 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 132 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 133 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 134 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 135 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 136 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 137 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 138 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 139 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 140 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 141 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 142 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 143 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 144 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 145 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 146 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 147 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 148 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 149 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 150 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 151 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 152 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 153 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 154 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 155 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 156 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 157 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 158 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 159 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 160 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 161 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 162 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 163 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 164 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 165 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 166 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 167 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 168 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 169 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 170 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 171 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 172 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 173 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 174 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 175 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 176 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 177 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 178 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 179 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 180 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 181 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 182 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 183 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 184 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 185 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 186 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 187 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 188 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 189 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 190 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 191 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 192 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 193 0.496 0.45907245 0.636656724 0.0755600483 0.02597481 0.001514519
#> 194 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 195 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 196 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 197 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 198 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 199 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 200 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 201 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 202 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 203 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 204 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 205 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 206 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 207 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 208 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 209 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 210 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 211 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 212 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 213 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 214 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 215 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 216 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 217 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 218 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 219 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 220 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 221 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 222 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 223 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 224 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 225 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 226 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 227 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 228 0.496 0.95759591 0.959727077 0.9862193637 1.00000000 0.999177807
#> 229 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 230 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 231 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 232 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 233 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 234 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 235 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 236 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 237 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 238 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 239 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 240 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 241 0.496 0.36839434 0.524751547 0.0000000000 0.00000000 0.000000000
#> 242 0.496 0.64117170 0.862240206 0.8736969788 0.87418616 0.874112455
#> 243 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 244 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#> 245 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 246 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 247 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 248 0.496 0.30502486 0.086944139 0.0840816281 0.08415168 0.083947110
#> 249 0.496 1.00000000 0.979507152 0.9800896768 0.98049332 0.980469640
#> 250 0.496 0.00000000 0.000000000 0.0000000000 0.00000000 0.000000000
#>             [,7]       [,8]         [,9]
#> 1   8.740665e-01 0.87406897 8.740663e-01
#> 2   9.803926e-01 0.98038990 9.803858e-01
#> 3   9.803926e-01 0.98038990 9.803858e-01
#> 4   8.740665e-01 0.87406897 8.740663e-01
#> 5   9.803926e-01 0.98038990 9.803858e-01
#> 6   9.803926e-01 0.98038990 9.803858e-01
#> 7   8.740665e-01 0.87406897 8.740663e-01
#> 8   0.000000e+00 0.00000000 0.000000e+00
#> 9   0.000000e+00 0.00000000 0.000000e+00
#> 10  9.803926e-01 0.98038990 9.803858e-01
#> 11  0.000000e+00 0.00000000 0.000000e+00
#> 12  9.803926e-01 0.98038990 9.803858e-01
#> 13  0.000000e+00 0.00000000 0.000000e+00
#> 14  8.740665e-01 0.87406897 8.740663e-01
#> 15  0.000000e+00 0.00000000 0.000000e+00
#> 16  8.399153e-02 0.08399119 8.398958e-02
#> 17  9.803926e-01 0.98038990 9.803858e-01
#> 18  0.000000e+00 0.00000000 0.000000e+00
#> 19  8.399153e-02 0.08399119 8.398958e-02
#> 20  8.399153e-02 0.08399119 8.398958e-02
#> 21  9.803926e-01 0.98038990 9.803858e-01
#> 22  8.399153e-02 0.08399119 8.398958e-02
#> 23  8.399153e-02 0.08399119 8.398958e-02
#> 24  0.000000e+00 0.00000000 0.000000e+00
#> 25  9.803926e-01 0.98038990 9.803858e-01
#> 26  9.803926e-01 0.98038990 9.803858e-01
#> 27  1.814101e-04 0.00000000 0.000000e+00
#> 28  8.740665e-01 0.87406897 8.740663e-01
#> 29  9.803926e-01 0.98038990 9.803858e-01
#> 30  0.000000e+00 0.00000000 0.000000e+00
#> 31  9.803926e-01 0.98038990 9.803858e-01
#> 32  9.803926e-01 0.98038990 9.803858e-01
#> 33  8.740665e-01 0.87406897 8.740663e-01
#> 34  0.000000e+00 0.00000000 0.000000e+00
#> 35  8.740665e-01 0.87406897 8.740663e-01
#> 36  9.803926e-01 0.98038990 9.803858e-01
#> 37  9.803926e-01 0.98038990 9.803858e-01
#> 38  0.000000e+00 0.00000000 0.000000e+00
#> 39  0.000000e+00 0.00000000 0.000000e+00
#> 40  0.000000e+00 0.00000000 0.000000e+00
#> 41  8.740665e-01 0.87406897 8.740663e-01
#> 42  8.399153e-02 0.08399119 8.398958e-02
#> 43  8.740665e-01 0.87406897 8.740663e-01
#> 44  8.399153e-02 0.08399119 8.398958e-02
#> 45  8.740665e-01 0.87406897 8.740663e-01
#> 46  0.000000e+00 0.00000000 0.000000e+00
#> 47  9.803926e-01 0.98038990 9.803858e-01
#> 48  8.740665e-01 0.87406897 8.740663e-01
#> 49  8.740665e-01 0.87406897 8.740663e-01
#> 50  8.399153e-02 0.08399119 8.398958e-02
#> 51  8.399153e-02 0.08399119 8.398958e-02
#> 52  0.000000e+00 0.00000000 0.000000e+00
#> 53  9.803926e-01 0.98038990 9.803858e-01
#> 54  8.740665e-01 0.87406897 8.740663e-01
#> 55  9.803926e-01 0.98038990 9.803858e-01
#> 56  9.803926e-01 0.98038990 9.803858e-01
#> 57  9.803926e-01 0.98038990 9.803858e-01
#> 58  8.399153e-02 0.08399119 8.398958e-02
#> 59  0.000000e+00 0.00000000 0.000000e+00
#> 60  9.803926e-01 0.98038990 9.803858e-01
#> 61  8.740665e-01 0.87406897 8.740663e-01
#> 62  8.740665e-01 0.87406897 8.740663e-01
#> 63  0.000000e+00 0.00000000 0.000000e+00
#> 64  0.000000e+00 0.00000000 0.000000e+00
#> 65  9.803926e-01 0.98038990 9.803858e-01
#> 66  9.803926e-01 0.98038990 9.803858e-01
#> 67  9.803926e-01 0.98038990 9.803858e-01
#> 68  0.000000e+00 0.00000000 0.000000e+00
#> 69  8.399153e-02 0.08399119 8.398958e-02
#> 70  0.000000e+00 0.00000000 0.000000e+00
#> 71  8.399153e-02 0.08399119 8.398958e-02
#> 72  8.399153e-02 0.08399119 8.398958e-02
#> 73  0.000000e+00 0.00000000 0.000000e+00
#> 74  8.399153e-02 0.08399119 8.398958e-02
#> 75  0.000000e+00 0.00000000 0.000000e+00
#> 76  8.740665e-01 0.87406897 8.740663e-01
#> 77  8.399153e-02 0.08399119 8.398958e-02
#> 78  8.399153e-02 0.08399119 8.398958e-02
#> 79  9.803926e-01 0.98038990 9.803858e-01
#> 80  8.740665e-01 0.87406897 8.740663e-01
#> 81  8.740665e-01 0.87406897 8.740663e-01
#> 82  8.740665e-01 0.87406897 8.740663e-01
#> 83  8.399153e-02 0.08399119 8.398958e-02
#> 84  8.740665e-01 0.87406897 8.740663e-01
#> 85  8.399153e-02 0.08399119 8.398958e-02
#> 86  8.399153e-02 0.08399119 8.398958e-02
#> 87  9.803926e-01 0.98038990 9.803858e-01
#> 88  9.986012e-01 0.99985040 1.000000e+00
#> 89  9.803926e-01 0.98038990 9.803858e-01
#> 90  8.740665e-01 0.87406897 8.740663e-01
#> 91  8.399153e-02 0.08399119 8.398958e-02
#> 92  0.000000e+00 0.00000000 0.000000e+00
#> 93  0.000000e+00 0.00000000 0.000000e+00
#> 94  0.000000e+00 0.00000000 0.000000e+00
#> 95  9.803926e-01 0.98038990 9.803858e-01
#> 96  8.399153e-02 0.08399119 8.398958e-02
#> 97  0.000000e+00 0.00000000 0.000000e+00
#> 98  0.000000e+00 0.00000000 0.000000e+00
#> 99  9.803926e-01 0.98038990 9.803858e-01
#> 100 8.740665e-01 0.87406897 8.740663e-01
#> 101 9.803926e-01 0.98038990 9.803858e-01
#> 102 8.399153e-02 0.08399119 8.398958e-02
#> 103 9.803926e-01 0.98038990 9.803858e-01
#> 104 0.000000e+00 0.00000000 0.000000e+00
#> 105 9.803926e-01 0.98038990 9.803858e-01
#> 106 0.000000e+00 0.00000000 0.000000e+00
#> 107 1.000000e+00 0.99995246 1.000000e+00
#> 108 8.740665e-01 0.87406897 8.740663e-01
#> 109 8.399153e-02 0.08399119 8.398958e-02
#> 110 0.000000e+00 0.00000000 0.000000e+00
#> 111 9.803926e-01 0.98038990 9.803858e-01
#> 112 9.803926e-01 0.98038990 9.803858e-01
#> 113 8.740665e-01 0.87406897 8.740663e-01
#> 114 8.399153e-02 0.08399119 8.398958e-02
#> 115 9.803926e-01 0.98038990 9.803858e-01
#> 116 8.740665e-01 0.87406897 8.740663e-01
#> 117 8.740665e-01 0.87406897 8.740663e-01
#> 118 0.000000e+00 0.00000000 0.000000e+00
#> 119 8.399153e-02 0.08399119 8.398958e-02
#> 120 9.803926e-01 0.98038990 9.803858e-01
#> 121 0.000000e+00 0.00000000 0.000000e+00
#> 122 8.740665e-01 0.87406897 8.740663e-01
#> 123 0.000000e+00 0.00000000 0.000000e+00
#> 124 8.399153e-02 0.08399119 8.398958e-02
#> 125 8.399153e-02 0.08399119 8.398958e-02
#> 126 9.803926e-01 0.98038990 9.803858e-01
#> 127 9.803926e-01 0.98038990 9.803858e-01
#> 128 9.803926e-01 0.98038990 9.803858e-01
#> 129 9.803926e-01 0.98038990 9.803858e-01
#> 130 8.399153e-02 0.08399119 8.398958e-02
#> 131 9.803926e-01 0.98038990 9.803858e-01
#> 132 8.399153e-02 0.08399119 8.398958e-02
#> 133 9.803926e-01 0.98038990 9.803858e-01
#> 134 8.740665e-01 0.87406897 8.740663e-01
#> 135 9.803926e-01 0.98038990 9.803858e-01
#> 136 8.740665e-01 0.87406897 8.740663e-01
#> 137 9.803926e-01 0.98038990 9.803858e-01
#> 138 8.740665e-01 0.87406897 8.740663e-01
#> 139 8.399153e-02 0.08399119 8.398958e-02
#> 140 0.000000e+00 0.00000000 0.000000e+00
#> 141 0.000000e+00 0.00000000 0.000000e+00
#> 142 0.000000e+00 0.00000000 0.000000e+00
#> 143 8.740665e-01 0.87406897 8.740663e-01
#> 144 8.740665e-01 0.87406897 8.740663e-01
#> 145 8.740665e-01 0.87406897 8.740663e-01
#> 146 8.399153e-02 0.08399119 8.398958e-02
#> 147 8.740665e-01 0.87406897 8.740663e-01
#> 148 9.803926e-01 0.98038990 9.803858e-01
#> 149 0.000000e+00 0.00000000 0.000000e+00
#> 150 8.740665e-01 0.87406897 8.740663e-01
#> 151 0.000000e+00 0.00000000 0.000000e+00
#> 152 0.000000e+00 0.00000000 0.000000e+00
#> 153 0.000000e+00 0.00000000 0.000000e+00
#> 154 9.803926e-01 0.98038990 9.803858e-01
#> 155 9.803926e-01 0.98038990 9.803858e-01
#> 156 9.803926e-01 0.98038990 9.803858e-01
#> 157 0.000000e+00 0.00000000 0.000000e+00
#> 158 8.399153e-02 0.08399119 8.398958e-02
#> 159 0.000000e+00 0.00000000 0.000000e+00
#> 160 0.000000e+00 0.00000000 0.000000e+00
#> 161 8.399153e-02 0.08399119 8.398958e-02
#> 162 0.000000e+00 0.00000000 0.000000e+00
#> 163 0.000000e+00 0.00000000 0.000000e+00
#> 164 9.803926e-01 0.98038990 9.803858e-01
#> 165 8.399153e-02 0.08399119 8.398958e-02
#> 166 8.740665e-01 0.87406897 8.740663e-01
#> 167 0.000000e+00 0.00000000 0.000000e+00
#> 168 8.740665e-01 0.87406897 8.740663e-01
#> 169 8.399153e-02 0.08399119 8.398958e-02
#> 170 9.803926e-01 0.98038990 9.803858e-01
#> 171 8.399153e-02 0.08399119 8.398958e-02
#> 172 8.740665e-01 0.87406897 8.740663e-01
#> 173 8.399153e-02 0.08399119 8.398958e-02
#> 174 9.803926e-01 0.98038990 9.803858e-01
#> 175 9.803926e-01 0.98038990 9.803858e-01
#> 176 0.000000e+00 0.00000000 0.000000e+00
#> 177 9.803926e-01 0.98038990 9.803858e-01
#> 178 8.740665e-01 0.87406897 8.740663e-01
#> 179 0.000000e+00 0.00000000 0.000000e+00
#> 180 0.000000e+00 0.00000000 0.000000e+00
#> 181 8.399153e-02 0.08399119 8.398958e-02
#> 182 9.803926e-01 0.98038990 9.803858e-01
#> 183 8.399153e-02 0.08399119 8.398958e-02
#> 184 9.803926e-01 0.98038990 9.803858e-01
#> 185 0.000000e+00 0.00000000 0.000000e+00
#> 186 9.803926e-01 0.98038990 9.803858e-01
#> 187 8.740665e-01 0.87406897 8.740663e-01
#> 188 8.740665e-01 0.87406897 8.740663e-01
#> 189 9.803926e-01 0.98038990 9.803858e-01
#> 190 0.000000e+00 0.00000000 0.000000e+00
#> 191 9.803926e-01 0.98038990 9.803858e-01
#> 192 0.000000e+00 0.00000000 0.000000e+00
#> 193 0.000000e+00 0.00000000 9.436896e-16
#> 194 8.740665e-01 0.87406897 8.740663e-01
#> 195 8.399153e-02 0.08399119 8.398958e-02
#> 196 9.803926e-01 0.98038990 9.803858e-01
#> 197 8.740665e-01 0.87406897 8.740663e-01
#> 198 9.803926e-01 0.98038990 9.803858e-01
#> 199 8.740665e-01 0.87406897 8.740663e-01
#> 200 0.000000e+00 0.00000000 0.000000e+00
#> 201 8.740665e-01 0.87406897 8.740663e-01
#> 202 8.399153e-02 0.08399119 8.398958e-02
#> 203 8.740665e-01 0.87406897 8.740663e-01
#> 204 8.740665e-01 0.87406897 8.740663e-01
#> 205 8.399153e-02 0.08399119 8.398958e-02
#> 206 0.000000e+00 0.00000000 0.000000e+00
#> 207 9.803926e-01 0.98038990 9.803858e-01
#> 208 8.399153e-02 0.08399119 8.398958e-02
#> 209 0.000000e+00 0.00000000 0.000000e+00
#> 210 9.803926e-01 0.98038990 9.803858e-01
#> 211 8.399153e-02 0.08399119 8.398958e-02
#> 212 0.000000e+00 0.00000000 0.000000e+00
#> 213 9.803926e-01 0.98038990 9.803858e-01
#> 214 9.803926e-01 0.98038990 9.803858e-01
#> 215 9.803926e-01 0.98038990 9.803858e-01
#> 216 8.399153e-02 0.08399119 8.398958e-02
#> 217 8.399153e-02 0.08399119 8.398958e-02
#> 218 9.803926e-01 0.98038990 9.803858e-01
#> 219 0.000000e+00 0.00000000 0.000000e+00
#> 220 8.399153e-02 0.08399119 8.398958e-02
#> 221 0.000000e+00 0.00000000 0.000000e+00
#> 222 8.399153e-02 0.08399119 8.398958e-02
#> 223 8.740665e-01 0.87406897 8.740663e-01
#> 224 0.000000e+00 0.00000000 0.000000e+00
#> 225 8.740665e-01 0.87406897 8.740663e-01
#> 226 8.740665e-01 0.87406897 8.740663e-01
#> 227 0.000000e+00 0.00000000 0.000000e+00
#> 228 1.000000e+00 0.99977443 1.000000e+00
#> 229 9.803926e-01 0.98038990 9.803858e-01
#> 230 9.803926e-01 0.98038990 9.803858e-01
#> 231 0.000000e+00 0.00000000 0.000000e+00
#> 232 9.803926e-01 0.98038990 9.803858e-01
#> 233 0.000000e+00 0.00000000 0.000000e+00
#> 234 9.803926e-01 0.98038990 9.803858e-01
#> 235 0.000000e+00 0.00000000 0.000000e+00
#> 236 0.000000e+00 0.00000000 0.000000e+00
#> 237 8.399153e-02 0.08399119 8.398958e-02
#> 238 9.803926e-01 0.98038990 9.803858e-01
#> 239 9.803926e-01 0.98038990 9.803858e-01
#> 240 8.740665e-01 0.87406897 8.740663e-01
#> 241 5.931236e-05 0.00000000 7.216450e-16
#> 242 8.740665e-01 0.87406897 8.740663e-01
#> 243 9.803926e-01 0.98038990 9.803858e-01
#> 244 0.000000e+00 0.00000000 0.000000e+00
#> 245 8.399153e-02 0.08399119 8.398958e-02
#> 246 8.399153e-02 0.08399119 8.398958e-02
#> 247 8.399153e-02 0.08399119 8.398958e-02
#> 248 8.399153e-02 0.08399119 8.398958e-02
#> 249 9.803926e-01 0.98038990 9.803858e-01
#> 250 0.000000e+00 0.00000000 0.000000e+00
modpls.B2$MissClassed
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#> [1,]  124   11   13   11   11   11   11   11   11
plsR(simbin1$Y,dataAstar2[,-1],10,typeVC="standard",MClassed=TRUE,verbose=FALSE)$InfCrit
#>                  AIC     Q2cum_Y LimQ2_Y        Q2_Y  PRESS_Y    RSS_Y
#> Nb_Comp_0  366.87968          NA      NA          NA       NA 62.49600
#> Nb_Comp_1  179.27607  0.52469907  0.0975  0.52469907 29.70441 29.27366
#> Nb_Comp_2  111.50920  0.63784869  0.0975  0.23805890 22.30480 22.14520
#> Nb_Comp_3   92.19471  0.61461931  0.0975 -0.06414274 23.56565 20.33539
#> Nb_Comp_4   92.80958  0.54850120  0.0975 -0.17156571 23.82425 20.22303
#> Nb_Comp_5   94.61150  0.46505988  0.0975 -0.18480962 23.96044 20.20702
#> Nb_Comp_6   96.58935  0.37811037  0.0975 -0.16254065 23.49148 20.20523
#> Nb_Comp_7   98.58715  0.28547861  0.0975 -0.14895208 23.21484 20.20505
#> Nb_Comp_8  100.58707  0.18505584  0.0975 -0.14054551 23.04478 20.20504
#> Nb_Comp_9  102.58707  0.07664463  0.0975 -0.13302900 22.89290 20.20504
#> Nb_Comp_10 104.58707 -0.03839400  0.0975 -0.12458760 22.72234 20.20504
#>                 R2_Y MissClassed R2_residY RSS_residY PRESS_residY   Q2_residY
#> Nb_Comp_0         NA         124        NA  249.00000           NA          NA
#> Nb_Comp_1  0.5315915          30 0.5315915  116.63372    118.34993  0.52469907
#> Nb_Comp_2  0.6456542           6 0.6456542   88.23211     88.86803  0.23805890
#> Nb_Comp_3  0.6746129           8 0.6746129   81.02138     93.89156 -0.06414274
#> Nb_Comp_4  0.6764108          11 0.6764108   80.57372     94.92187 -0.17156571
#> Nb_Comp_5  0.6766671          10 0.6766671   80.50990     95.46452 -0.18480962
#> Nb_Comp_6  0.6766957          10 0.6766957   80.50277     93.59604 -0.16254065
#> Nb_Comp_7  0.6766985          10 0.6766985   80.50206     92.49383 -0.14895208
#> Nb_Comp_8  0.6766986          10 0.6766986   80.50204     91.81627 -0.14054551
#> Nb_Comp_9  0.6766987          10 0.6766987   80.50203     91.21114 -0.13302900
#> Nb_Comp_10 0.6766987          10 0.6766987   80.50203     90.53159 -0.12458760
#>             LimQ2 Q2cum_residY  AIC.std   DoF.dof sigmahat.dof    AIC.dof
#> Nb_Comp_0      NA           NA 712.4673  1.000000    0.5009870 0.25199190
#> Nb_Comp_1  0.0975   0.52469907 524.8637  2.621153    0.3433059 0.11956605
#> Nb_Comp_2  0.0975   0.63784869 457.0968  3.000068    0.2988230 0.09072392
#> Nb_Comp_3  0.0975   0.61461931 437.7823 21.497047    0.2976680 0.09657974
#> Nb_Comp_4  0.0975   0.54850120 438.3972 22.449662    0.2974625 0.09678361
#> Nb_Comp_5  0.0975   0.46505988 440.1991 22.773576    0.2975556 0.09695892
#> Nb_Comp_6  0.0975   0.37811037 442.1769 22.958265    0.2976629 0.09709431
#> Nb_Comp_7  0.0975   0.28547861 444.1747 23.044744    0.2977180 0.09716095
#> Nb_Comp_8  0.0975   0.18505584 446.1747 23.145754    0.2977840 0.09723982
#> Nb_Comp_9  0.0975   0.07664463 448.1747 23.186981    0.2978109 0.09727204
#> Nb_Comp_10 0.0975  -0.03839400 450.1747 23.388199    0.2979425 0.09742948
#>              BIC.dof  GMDL.dof DoF.naive sigmahat.naive  AIC.naive  BIC.naive
#> Nb_Comp_0  0.2555273 -167.2823         1      0.5009870 0.25199190 0.25552728
#> Nb_Comp_1  0.1239175 -257.0188         2      0.3435680 0.11898326 0.12230862
#> Nb_Comp_2  0.0944974 -290.3040         3      0.2994272 0.09073255 0.09452122
#> Nb_Comp_3  0.1234101 -257.2237         4      0.2875138 0.08398681 0.08864439
#> Nb_Comp_4  0.1247642 -255.9517         5      0.2873030 0.08419385 0.09000729
#> Nb_Comp_5  0.1253610 -255.4013         6      0.2877771 0.08480321 0.09180238
#> Nb_Comp_6  0.1257474 -255.0480         7      0.2883558 0.08547725 0.09367583
#> Nb_Comp_7  0.1259326 -254.8792         8      0.2889497 0.08616367 0.09557211
#> Nb_Comp_8  0.1261504 -254.6812         9      0.2895485 0.08685653 0.09748493
#> Nb_Comp_9  0.1262394 -254.6005        10      0.2901511 0.08755518 0.09941372
#> Nb_Comp_10 0.1266740 -254.2069        11      0.2907575 0.08825968 0.10135865
#>            GMDL.naive
#> Nb_Comp_0   -167.2823
#> Nb_Comp_1   -258.3373
#> Nb_Comp_2   -289.8052
#> Nb_Comp_3   -297.3647
#> Nb_Comp_4   -295.2257
#> Nb_Comp_5   -292.6062
#> Nb_Comp_6   -289.9866
#> Nb_Comp_7   -287.4310
#> Nb_Comp_8   -284.9391
#> Nb_Comp_9   -282.5049
#> Nb_Comp_10  -280.1229
cv.modpls.B2<-cv.plsR(Y~.,data=simbin1,2,NK=100,verbose=FALSE)
res.cv.modpls.B2<-cvtable(summary(cv.modpls.B2,MClassed=TRUE))
#> ____************************************************____
#> Error in eval(mf, parent.frame()): object 'simbin1' not found
#Only one component found by repeated CV missclassed criterion
plot(res.cv.modpls.B2)
#> Error: object 'res.cv.modpls.B2' not found

rm(list=c("dimX","Astar","dataAstar2","modpls.A2","cv.modpls.A2",
"res.cv.modpls.A2","simbin1","modpls.B2","cv.modpls.B2","res.cv.modpls.B2"))
#> Warning: object 'res.cv.modpls.A2' not found
#> Warning: object 'res.cv.modpls.B2' not found
# }
```
