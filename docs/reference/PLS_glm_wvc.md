# Light version of PLS_glm for cross validation purposes

Light version of `PLS_glm` for cross validation purposes either on
complete or incomplete datasets.

## Usage

``` r
PLS_glm_wvc(
  dataY,
  dataX,
  nt = 2,
  dataPredictY = dataX,
  modele = "pls",
  family = NULL,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepstd.coeffs = FALSE,
  tol_Xi = 10^(-12),
  weights,
  method = "logistic",
  fit_backend = "stats",
  verbose = TRUE
)
```

## Arguments

- dataY:

  response (training) dataset

- dataX:

  predictor(s) (training) dataset

- nt:

  number of components to be extracted

- dataPredictY:

  predictor(s) (testing) dataset

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

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since non always possible for glm
  responses.

- keepcoeffs:

  whether the coefficients of the linear fit on link scale of
  unstandardized eXplanatory variables should be returned or not.

- keepstd.coeffs:

  whether the coefficients of the linear fit on link scale of
  standardized eXplanatory variables should be returned or not.

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- method:

  logistic, probit, complementary log-log or cauchit (corresponding to a
  Cauchy latent variable).

- fit_backend:

  backend used for repeated non-ordinal score-space GLM fits. Use
  `"stats"` for the compatibility path or `"fastglm"` to opt into the
  accelerated complete-data backend. Unsupported cases fall back to
  `"stats"` with a warning.

- verbose:

  should info messages be displayed ?

## Value

- valsPredict:

  `nrow(dataPredictY) * nt` matrix of the predicted values

- list("coeffs"):

  If the coefficients of the eXplanatory variables were requested:  
  i.e. `keepcoeffs=TRUE`.  
  `ncol(dataX) * 1` matrix of the coefficients of the the eXplanatory
  variables

## Details

This function is called by
[`PLS_glm_kfoldcv_formula`](https://fbertran.github.io/plsRglm/reference/cv.plsRglm.md)
in order to perform cross-validation either on complete or incomplete
datasets.

There are seven different predefined models with predefined link
functions available :

- list("\\pls\\"):

  ordinary pls models

- list("\\pls-glm-Gamma\\"):

  glm gaussian with inverse link pls models

- list("\\pls-glm-gaussian\\"):

  glm gaussian with identity link pls models

- list("\\pls-glm-inverse-gamma\\"):

  glm binomial with square inverse link pls models

- list("\\pls-glm-logistic\\"):

  glm binomial with logit link pls models

- list("\\pls-glm-poisson\\"):

  glm poisson with log link pls models

- list("\\pls-glm-polr\\"):

  glm polr with logit link pls models

Using the `"family="` option and setting `"modele=pls-glm-family"`
allows changing the family and link function the same way as for the
[`glm`](https://rdrr.io/r/stats/glm.html) function. As a consequence
user-specified families can also be used.

- The :

  accepts the links (as names) `identity`, `log` and `inverse`.

- list("gaussian"):

  accepts the links (as names) `identity`, `log` and `inverse`.

- family:

  accepts the links (as names) `identity`, `log` and `inverse`.

- The :

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- list("binomial"):

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- family:

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- The :

  accepts the links `inverse`, `identity` and `log`.

- list("Gamma"):

  accepts the links `inverse`, `identity` and `log`.

- family:

  accepts the links `inverse`, `identity` and `log`.

- The :

  accepts the links `log`, `identity`, and `sqrt`.

- list("poisson"):

  accepts the links `log`, `identity`, and `sqrt`.

- family:

  accepts the links `log`, `identity`, and `sqrt`.

- The :

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- list("inverse.gaussian"):

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- family:

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- The :

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- list("quasi"):

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- family:

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- The function :

  can be used to create a power link function.

- list("power"):

  can be used to create a power link function.

Non-NULL weights can be used to indicate that different observations
have different dispersions (with the values in weights being inversely
proportional to the dispersions); or equivalently, when the elements of
weights are positive integers w_i, that each response y_i is the mean of
w_i unit-weight observations.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`PLS_glm`](https://fbertran.github.io/plsRglm/reference/plsRglm.md) for
more detailed results,
[`PLS_glm_kfoldcv`](https://fbertran.github.io/plsRglm/reference/cv.plsRglm.md)
for cross-validating models and
[`PLS_lm_wvc`](https://fbertran.github.io/plsRglm/reference/PLS_lm_wvc.md)
for the same function dedicated to plsR models

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",
dataPredictY=XCornell[1,])
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
#> $valsPredict
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",
family=gaussian(),dataPredictY=XCornell[1,], verbose=FALSE)
#> $valsPredict
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-gaussian",
dataPredictY=XCornell[1,], verbose=FALSE)
#> $valsPredict
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-family",
family=gaussian(),dataPredictY=XCornell[1,], verbose=FALSE)
#> $valsPredict
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
rm("XCornell","yCornell")

# \donttest{
## With an incomplete dataset (X[1,2] is NA)
data(pine)
ypine <- pine[,11]
data(XpineNAX21)
PLS_glm_wvc(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian")
#> ____************************************************____
#> Only naive DoF can be used with missing data
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> ____There are some NAs in X but not in Y____
#> ____Predicting X with NA in X and not in Y____
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
#> ****________________________________________________****
#> 
#> $valsPredict
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [2,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [3,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [4,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [5,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [6,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [7,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [8,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>  [9,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [10,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [11,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [12,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [13,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [14,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [15,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [16,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [17,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [18,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [19,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [20,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [21,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [22,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [23,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [24,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [25,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [26,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [27,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [28,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [29,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [30,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [31,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [32,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> [33,]   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> 
rm("XpineNAX21","ypine")
#> Warning: object 'XpineNAX21' not found

data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_glm_wvc(ypine,Xpine,10,modele="pls", verbose=FALSE)
#> $valsPredict
#>           [,1]        [,2]        [,3]        [,4]         [,5]        [,6]
#> 1   1.57673729  2.00433040  2.00216849  2.01939633  1.937245401  1.88550891
#> 2   1.01203030  1.15199511  1.11571861  1.24570213  1.235802532  1.20179188
#> 3   1.48993594  1.24149458  1.36818432  1.53585542  1.486704711  1.51671911
#> 4   0.82287026  0.63250520  0.93956109  1.04161870  1.007232341  1.04125064
#> 5   1.03843648  0.73560935  0.56129947  0.58902449  0.528441708  0.57134580
#> 6   1.35675180  1.40730081  1.25451432  1.37731103  1.366862539  1.39997555
#> 7   0.32014006 -0.32532009 -0.14347163 -0.17938500 -0.165426534 -0.11535438
#> 8   0.71839021  0.73099520  0.30866945  0.39516648  0.412751496  0.32418945
#> 9   1.48896024  1.52799376  1.67587131  1.60139990  1.652717167  1.65007715
#> 10  0.99804354  0.88959103  0.97569244  0.86352520  0.947912813  0.94847153
#> 11  0.43405362  0.86608054  0.97012047  0.64182858  0.616856901  0.63882421
#> 12 -0.07021385  0.26999887  0.59312200  0.92157715  1.008253831  0.98380367
#> 13  1.76424139  1.90090137  2.11336407  2.19169300  2.166416346  2.19486369
#> 14  1.37638458  1.70611149  1.58504677  1.84688715  1.857927378  1.90455582
#> 15  1.29114690  1.38753932  1.35365960  1.43462512  1.476219895  1.51027552
#> 16 -0.31402446 -0.47052602 -0.50369169 -0.48552205 -0.489472457 -0.51351160
#> 17  0.07868727  0.09889902 -0.01564596  0.15763382  0.083278591  0.10574381
#> 18  1.03444868  1.14007872  1.09595610  0.88458167  0.852892800  0.78377367
#> 19  1.08137124  0.67811651  0.73915351  0.82968208  0.879143668  0.88040097
#> 20  0.27838475  0.03080907 -0.03681897 -0.06160100 -0.054332827 -0.07763891
#> 21  0.04732997  0.24640318  0.51964685  0.71563687  0.653313657  0.61594736
#> 22  0.64993047  0.65252960  0.76864397  0.76705111  0.732874483  0.71165796
#> 23  0.91576469  0.95152885  0.70877390  0.68206508  0.726537887  0.71111241
#> 24  0.28714240  0.25747771  0.32732647  0.36607751  0.340816493  0.27266379
#> 25  0.55621083  0.20897499 -0.13621852 -0.01161858 -0.029168653 -0.04249364
#> 26  0.92489518  0.94649310  0.77909889  0.64060202  0.670826485  0.67364464
#> 27  0.01702572  0.38181268  0.20327363  0.07976920  0.023726459  0.15754257
#> 28  1.54352449  1.13003525  1.14773845  0.98446655  0.935861794  0.90956750
#> 29  1.45181367  1.42783969  1.54487387  1.46033916  1.495705868  1.51049933
#> 30  0.73123225  0.66412607  0.84569291  0.49408354  0.482259156  0.42326984
#> 31  0.98197420  1.24443724  1.14351551  1.31718138  1.406786445  1.39942166
#> 32  0.66850881  0.94811355  0.92220097  0.50456716  0.531183973  0.55048739
#> 33  0.21787107  0.10572384  0.04295933 -0.08122119 -0.008152348  0.04161271
#>           [,7]         [,8]        [,9]       [,10]
#> 1   1.95676227  1.948450849  1.93896093  1.94144962
#> 2   1.21592652  1.235442425  1.24794773  1.24303207
#> 3   1.54977163  1.524520187  1.53266664  1.53900851
#> 4   1.05581321  1.021862826  1.03178760  1.03548430
#> 5   0.58624206  0.562391104  0.56933856  0.58993331
#> 6   1.38175170  1.360508864  1.35664763  1.36575593
#> 7  -0.07128741 -0.047428997 -0.05482975 -0.06507565
#> 8   0.35423193  0.306999026  0.28300642  0.26870841
#> 9   1.63319069  1.631509959  1.64240206  1.64555642
#> 10  0.90453538  0.837015730  0.82600711  0.82493397
#> 11  0.66440603  0.619480777  0.61806414  0.61605314
#> 12  1.06477362  1.031255436  1.02700955  1.02592619
#> 13  2.16399005  2.163911302  2.15049944  2.14429918
#> 14  1.84909183  1.845501860  1.84534215  1.84361941
#> 15  1.53220512  1.567310225  1.56848600  1.56902383
#> 16 -0.51165239 -0.500336508 -0.52095379 -0.50662867
#> 17  0.01828375 -0.002607768  0.01667444 -0.02206645
#> 18  0.88281941  0.918670514  0.94195909  0.93284893
#> 19  0.88632519  0.909007733  0.91379654  0.90631555
#> 20 -0.03748213 -0.025963454 -0.03697219 -0.04597628
#> 21  0.56329746  0.606867697  0.60725951  0.63204655
#> 22  0.71534817  0.734795585  0.73328612  0.71985709
#> 23  0.66101344  0.690882417  0.69856666  0.69360847
#> 24  0.17883038  0.221072111  0.21360194  0.22005958
#> 25 -0.02445300 -0.070467998 -0.05897223 -0.03210573
#> 26  0.69268993  0.763387338  0.78043784  0.79028544
#> 27  0.19700010  0.215177839  0.20405047  0.19800511
#> 28  0.87932761  0.893834308  0.86843397  0.85997782
#> 29  1.52036274  1.552065698  1.53445170  1.53023569
#> 30  0.38164488  0.326448839  0.35187384  0.34979385
#> 31  1.37872533  1.385981563  1.39371317  1.39564352
#> 32  0.49904917  0.493233383  0.48196396  0.49873722
#> 33  0.04746533  0.049219132  0.06349274  0.06165365
#> 
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-Gamma", verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-family",family=Gamma(), verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-gaussian", verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-family",family=gaussian(log), verbose=FALSE)
#> Warning: glm.fit: algorithm did not converge
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
PLS_glm_wvc(round(ypine),Xpine,10,modele="pls-glm-poisson", verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
PLS_glm_wvc(round(ypine),Xpine,10,modele="pls-glm-family",family=poisson(log), verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
rm(list=c("pine","ypine","Xpine"))
#> Warning: object 'pine' not found


data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm_wvc(yCornell,XCornell,10,modele="pls-glm-inverse.gaussian", verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6]
#> 1    NA   NA   NA   NA   NA   NA
#> 2    NA   NA   NA   NA   NA   NA
#> 3    NA   NA   NA   NA   NA   NA
#> 4    NA   NA   NA   NA   NA   NA
#> 5    NA   NA   NA   NA   NA   NA
#> 6    NA   NA   NA   NA   NA   NA
#> 7    NA   NA   NA   NA   NA   NA
#> 8    NA   NA   NA   NA   NA   NA
#> 9    NA   NA   NA   NA   NA   NA
#> 10   NA   NA   NA   NA   NA   NA
#> 11   NA   NA   NA   NA   NA   NA
#> 12   NA   NA   NA   NA   NA   NA
#> 
PLS_glm_wvc(yCornell,XCornell,10,modele="pls-glm-family",
family=inverse.gaussian(), verbose=FALSE)
#> $valsPredict
#>    [,1] [,2] [,3] [,4] [,5] [,6]
#> 1    NA   NA   NA   NA   NA   NA
#> 2    NA   NA   NA   NA   NA   NA
#> 3    NA   NA   NA   NA   NA   NA
#> 4    NA   NA   NA   NA   NA   NA
#> 5    NA   NA   NA   NA   NA   NA
#> 6    NA   NA   NA   NA   NA   NA
#> 7    NA   NA   NA   NA   NA   NA
#> 8    NA   NA   NA   NA   NA   NA
#> 9    NA   NA   NA   NA   NA   NA
#> 10   NA   NA   NA   NA   NA   NA
#> 11   NA   NA   NA   NA   NA   NA
#> 12   NA   NA   NA   NA   NA   NA
#> 
rm(list=c("XCornell","yCornell"))


data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",
dataPredictY=XCornell[1,], verbose=FALSE)
#> $valsPredict
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-gaussian",
dataPredictY=XCornell[1,], verbose=FALSE)
#> $valsPredict
#>   [,1] [,2] [,3]
#> 1   NA   NA   NA
#> 
rm("XCornell","yCornell")

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
PLS_glm(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic",typeVC="none", verbose=FALSE)$InfCrit
#>                 AIC      BIC Missclassed Chi2_Pearson_Y    RSS_Y      R2_Y
#> Nb_Comp_0  145.8283 148.4727          49      104.00000 25.91346        NA
#> Nb_Comp_1  118.1398 123.4285          28      100.53823 19.32272 0.2543365
#> Nb_Comp_2  109.9553 117.8885          26       99.17955 17.33735 0.3309519
#> Nb_Comp_3  105.1591 115.7366          22      123.37836 15.58198 0.3986915
#> Nb_Comp_4  103.8382 117.0601          21      114.77551 15.14046 0.4157299
#> Nb_Comp_5  104.7338 120.6001          21      105.35382 15.08411 0.4179043
#> Nb_Comp_6  105.6770 124.1878          21       98.87767 14.93200 0.4237744
#> Nb_Comp_7  107.2828 128.4380          20       97.04072 14.87506 0.4259715
#> Nb_Comp_8  109.0172 132.8167          22       98.90110 14.84925 0.4269676
#> Nb_Comp_9  110.9354 137.3793          21      100.35563 14.84317 0.4272022
#> Nb_Comp_10 112.9021 141.9904          20      102.85214 14.79133 0.4292027
#>             R2_residY RSS_residY
#> Nb_Comp_0          NA   25.91346
#> Nb_Comp_1   -6.004879  181.52066
#> Nb_Comp_2   -9.617595  275.13865
#> Nb_Comp_3  -12.332217  345.48389
#> Nb_Comp_4  -15.496383  427.47839
#> Nb_Comp_5  -15.937183  438.90105
#> Nb_Comp_6  -16.700929  458.69233
#> Nb_Comp_7  -16.908851  464.08033
#> Nb_Comp_8  -17.555867  480.84675
#> Nb_Comp_9  -17.834439  488.06552
#> Nb_Comp_10 -17.999267  492.33678
PLS_glm_wvc(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic", keepcoeffs=TRUE, verbose=FALSE)
#> $valsPredict
#>     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> 1     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 2     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 3     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 4     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 5     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 6     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 7     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 8     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 9     NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 10    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 11    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 12    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 13    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 14    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 15    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 16    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 17    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 18    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 19    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 20    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 21    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 22    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 23    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 24    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 25    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 26    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 27    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 28    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 29    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 30    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 31    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 32    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 33    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 34    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 35    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 36    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 37    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 38    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 39    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 40    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 41    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 42    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 43    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 44    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 45    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 46    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 47    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 48    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 49    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 50    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 51    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 52    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 53    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 54    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 55    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 56    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 57    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 58    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 59    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 60    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 61    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 62    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 63    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 64    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 65    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 66    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 67    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 68    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 69    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 70    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 71    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 72    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 73    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 74    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 75    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 76    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 77    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 78    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 79    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 80    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 81    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 82    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 83    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 84    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 85    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 86    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 87    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 88    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 89    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 90    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 91    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 92    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 93    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 94    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 95    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 96    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 97    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 98    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 99    NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 100   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 101   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 102   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 103   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 104   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA
#> 
#> $coeffs
#>                       [,1]
#> tempConstante -2.276982302
#>               -1.068275295
#>                3.509231595
#>               -1.651869135
#>                2.207538418
#>                0.568523938
#>               -0.059691869
#>               -0.214529856
#>               -1.405223273
#>                0.396973880
#>               -0.782167532
#>                0.677591817
#>               -0.972259676
#>                0.650745841
#>                0.723667343
#>                0.477540145
#>                0.638755948
#>                1.666070158
#>               -0.005938234
#>                0.482766293
#>               -0.904425334
#>                0.300460249
#>                1.367992779
#>               -1.201977825
#>               -1.536120691
#>               -1.983144986
#>                1.544435411
#>                1.410302156
#>               -0.495400138
#>                0.454129717
#>                1.240250301
#>               -0.222933455
#>               -2.822712745
#>                0.026369914
#> 
rm("Xaze_compl","yaze_compl")
# }
```
