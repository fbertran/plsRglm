# Light version of PLS_lm for cross validation purposes

Light version of `PLS_lm` for cross validation purposes either on
complete or incomplete datasets.

## Usage

``` r
PLS_lm_wvc(
  dataY,
  dataX,
  nt = 2,
  dataPredictY = dataX,
  modele = "pls",
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepstd.coeffs = FALSE,
  tol_Xi = 10^(-12),
  weights,
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

  name of the PLS model to be fitted, only (`"pls"` available for this
  fonction.

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since non always possible for glm
  responses.

- keepcoeffs:

  whether the coefficients of unstandardized eXplanatory variables
  should be returned or not.

- keepstd.coeffs:

  whether the coefficients of standardized eXplanatory variables should
  be returned or not.

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

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
[`PLS_lm_kfoldcv`](https://fbertran.github.io/plsRglm/reference/cv.plsR.md)
in order to perform cross-validation either on complete or incomplete
datasets.

Non-NULL weights can be used to indicate that different observations
have different dispersions (with the values in weights being inversely
proportional to the dispersions); or equivalently, when the elements of
weights are positive integers w_i, that each response y_i is the mean of
w_i unit-weight observations.

## Note

Use
[`PLS_lm_kfoldcv`](https://fbertran.github.io/plsRglm/reference/cv.plsR.md)
for a wrapper in view of cross-validation.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`PLS_lm`](https://fbertran.github.io/plsRglm/reference/plsR.md) for
more detailed results,
[`PLS_lm_kfoldcv`](https://fbertran.github.io/plsRglm/reference/cv.plsR.md)
for cross-validating models and
[`PLS_glm_wvc`](https://fbertran.github.io/plsRglm/reference/PLS_glm_wvc.md)
for the same function dedicated to plsRglm models

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm_wvc(dataY=yCornell,dataX=XCornell,nt=3,dataPredictY=XCornell[1,])
#> ____************************************************____
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
#> $valsPredict
#>       [,1]     [,2]     [,3]
#> 1 95.03164 96.49529 97.55864
#> 
PLS_lm_wvc(dataY=yCornell[-c(1,2)],dataX=XCornell[-c(1,2),],nt=3,dataPredictY=XCornell[c(1,2),],
verbose=FALSE)
#> $valsPredict
#>       [,1]     [,2]     [,3]
#> 1 93.21260 93.91414 95.73889
#> 2 94.71773 95.68684 96.79566
#> 
PLS_lm_wvc(dataY=yCornell[-c(1,2)],dataX=XCornell[-c(1,2),],nt=3,dataPredictY=XCornell[c(1,2),],
keepcoeffs=TRUE, verbose=FALSE)
#> $valsPredict
#>       [,1]     [,2]     [,3]
#> 1 93.21260 93.91414 95.73889
#> 2 94.71773 95.68684 96.79566
#> 
#> $coeffs
#>                     [,1]
#> tempConstante  90.036430
#> X1             -8.146409
#> X2             -4.000226
#> X3            -13.728783
#> X4             -5.883018
#> X5              6.937067
#> X6             10.148166
#> X7            -29.571048
#> 
rm("XCornell","yCornell")

## With an incomplete dataset (X[1,2] is NA)
data(pine)
ypine <- pine[,11]
data(XpineNAX21)
PLS_lm_wvc(dataY=ypine[-1],dataX=XpineNAX21[-1,],nt=3, verbose=FALSE)
#> $valsPredict
#>           [,1]        [,2]          [,3]
#> 2   0.94076404  1.05174832  1.0374616113
#> 3   1.39047885  1.13453980  1.3242356673
#> 4   0.79354603  0.61363359  0.9634719738
#> 5   0.97860870  0.65429731  0.5150093803
#> 6   1.27787513  1.34466033  1.2141683374
#> 7   0.38074833 -0.19960869 -0.0297875729
#> 8   0.63365198  0.54155986  0.2091351671
#> 9   1.43295859  1.55738575  1.6867033620
#> 10  0.99298825  0.98045144  1.0732360517
#> 11  0.39646989  0.77799528  0.8878866176
#> 12 -0.06668771  0.22196126  0.6489885307
#> 13  1.66741434  1.86115065  2.0625005103
#> 14  1.30880185  1.68986944  1.5677761570
#> 15  1.23730000  1.38593444  1.3409666162
#> 16 -0.25157291 -0.43528542 -0.4743325478
#> 17  0.11040475  0.12039192  0.0007238673
#> 18  0.92497736  0.93090871  0.9139450905
#> 19  1.07369519  0.75362958  0.8101882969
#> 20  0.28649448  0.01372797 -0.0381486634
#> 21  0.07611695  0.25760359  0.4947474762
#> 22  0.62440377  0.61179280  0.7213780672
#> 23  0.90767675  1.00727900  0.7150222566
#> 24  0.32767392  0.33311540  0.3319260445
#> 25  0.52077680  0.11553149 -0.1540043757
#> 26  0.90533669  0.96712444  0.7408564424
#> 27  0.03199415  0.35664666  0.1497469324
#> 28  1.47696322  1.10783727  1.0880798820
#> 29  1.39868360  1.45100576  1.5321709218
#> 30  0.71071212  0.65722183  0.8296817343
#> 31  0.95366454  1.27340221  1.1747095867
#> 32  0.67802770  1.02264538  0.9118337341
#> 33  0.27905267  0.23984263  0.1497228439
#> 
PLS_lm_wvc(dataY=ypine[-1],dataX=XpineNAX21[-1,],nt=3,dataPredictY=XpineNAX21[1,], verbose=FALSE)
#> list()
PLS_lm_wvc(dataY=ypine[-2],dataX=XpineNAX21[-2,],nt=3,dataPredictY=XpineNAX21[2,], verbose=FALSE)
#> $valsPredict
#>        [,1]     [,2]     [,3]
#> 2 0.9914237 1.127694 1.076871
#> 
PLS_lm_wvc(dataY=ypine,dataX=XpineNAX21,nt=3, verbose=FALSE)
#> $valsPredict
#>                [,1]        [,2]        [,3]
#>  [1,]  1.4589428372  2.11593441  2.26945463
#>  [2,]  0.9829728848  1.18076301  1.03620773
#>  [3,]  1.5465793485  1.18650144  1.27588027
#>  [4,]  0.8655200475  0.46560150  0.87695905
#>  [5,]  1.1072893206  0.85251283  0.57894070
#>  [6,]  1.3475354913  1.50304781  1.20343678
#>  [7,]  0.4642031730 -0.42135911 -0.09955758
#>  [8,]  0.7227037635  0.98797643  0.21698630
#>  [9,]  1.4817373396  1.44514851  1.73358055
#> [10,]  1.0254800667  0.83934144  1.03617380
#> [11,]  0.3387272475  0.80689483  1.16252267
#> [12,] -0.1404993319  0.07995556  0.28848700
#> [13,]  1.7330197324  1.78515568  2.09996156
#> [14,]  1.3042461768  1.77459300  1.45332343
#> [15,]  1.2715747058  1.40954157  1.32075175
#> [16,] -0.2777168878 -0.44870259 -0.51434295
#> [17,]  0.0732828651  0.16650220 -0.07993936
#> [18,]  1.0136353777  1.18186572  1.21367764
#> [19,]  1.1731943384  0.65036199  0.69946245
#> [20,]  0.3361806444  0.07862672 -0.02520022
#> [21,]  0.0005399781  0.08838604  0.42008495
#> [22,]  0.6493952920  0.59285401  0.78398392
#> [23,]  0.9096512309  1.09049865  0.76180402
#> [24,]  0.2922013348  0.21706236  0.34621654
#> [25,]  0.6382767422  0.42340484 -0.21649161
#> [26,]  0.9211966658  1.04536567  0.90202594
#> [27,] -0.0650662222  0.48221488  0.31920563
#> [28,]  1.6349387353  1.13681799  1.30914490
#> [29,]  1.4578617140  1.36487443  1.62299471
#> [30,]  0.7469063021  0.56602059  1.06652525
#> [31,]  0.9273242679  1.29767879  1.02671404
#> [32,]  0.6059700763  0.95665914  1.21783413
#> [33,]  0.2450768211  0.13732755  0.12831479
#> 
rm("ypine")
```
