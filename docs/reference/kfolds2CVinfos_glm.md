# Extracts and computes information criteria and fits statistics for k-fold cross validated partial least squares glm models

This function extracts and computes information criteria and fits
statistics for k-fold cross validated partial least squares glm models
for both formula or classic specifications of the model.

## Usage

``` r
kfolds2CVinfos_glm(pls_kfolds, MClassed = FALSE, verbose = TRUE)
```

## Arguments

- pls_kfolds:

  an object computed using
  [`cv.plsRglm`](https://fbertran.github.io/plsRglm/reference/cv.plsRglm.md)

- MClassed:

  should number of miss classed be computed ?

- verbose:

  should infos be displayed ?

## Value

- list:

  table of fit statistics for first group partition

- list():

  ...

- list:

  table of fit statistics for last group partition

## Details

The Mclassed option should only set to `TRUE` if the response is binary.

## Note

Use [`summary`](https://rdrr.io/r/base/summary.html) and
[`cv.plsRglm`](https://fbertran.github.io/plsRglm/reference/cv.plsRglm.md)
instead.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/47/>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.md),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.md),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.md),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.md)
and
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.md)
to extract and transforms results from k-fold cross-validation.

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data(Cornell)
summary(cv.plsRglm(Y~.,data=Cornell,
nt=6,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls",verbose=FALSE),MClassed=TRUE)
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
#>                AIC MissClassed CV_MissClassed    Q2cum_Y LimQ2_Y       Q2_Y
#> Nb_Comp_0 82.01205          12             NA         NA      NA         NA
#> Nb_Comp_1 53.15173          12             12  0.8809146  0.0975  0.8809146
#> Nb_Comp_2 41.08283          12             12  0.8619560  0.0975 -0.1592015
#> Nb_Comp_3 32.06411          12             12  0.7471041  0.0975 -0.8319956
#> Nb_Comp_4 33.76477          12             12 -0.2159389  0.0975 -3.8080607
#> Nb_Comp_5 33.34373          12             12 -5.9182568  0.0975 -4.6896417
#> Nb_Comp_6 35.25533          12             NA         NA  0.0975         NA
#>            PRESS_Y      RSS_Y      R2_Y    AIC.std  DoF.dof sigmahat.dof
#> Nb_Comp_0       NA 467.796667        NA  37.010388 1.000000    6.5212706
#> Nb_Comp_1 55.70774  35.742486 0.9235940   8.150064 2.740749    1.8665281
#> Nb_Comp_2 41.43274  11.066606 0.9763431  -3.918831 5.085967    1.1825195
#> Nb_Comp_3 20.27397   4.418081 0.9905556 -12.937550 5.121086    0.7488308
#> Nb_Comp_4 21.24240   4.309235 0.9907882 -11.236891 5.103312    0.7387162
#> Nb_Comp_5 24.51801   3.521924 0.9924713 -11.657929 6.006316    0.7096382
#> Nb_Comp_6       NA   3.496074 0.9925265  -9.746328 7.000002    0.7633343
#>              AIC.dof    BIC.dof GMDL.dof DoF.naive sigmahat.naive  AIC.naive
#> Nb_Comp_0 46.0708838 47.7893514 27.59461         1      6.5212706 46.0708838
#> Nb_Comp_1  4.5699686  4.9558156 21.34020         2      1.8905683  4.1699567
#> Nb_Comp_2  2.1075461  2.3949331 27.40202         3      1.1088836  1.5370286
#> Nb_Comp_3  0.8467795  0.9628191 24.40842         4      0.7431421  0.7363469
#> Nb_Comp_4  0.8232505  0.9357846 24.23105         5      0.7846050  0.8721072
#> Nb_Comp_5  0.7976101  0.9198348 28.21184         6      0.7661509  0.8804809
#> Nb_Comp_6  0.9711321  1.1359501 33.18347         7      0.8361907  1.1070902
#>            BIC.naive GMDL.naive
#> Nb_Comp_0 47.7893514   27.59461
#> Nb_Comp_1  4.4588195   18.37545
#> Nb_Comp_2  1.6860917   17.71117
#> Nb_Comp_3  0.8256118   19.01033
#> Nb_Comp_4  0.9964867   24.16510
#> Nb_Comp_5  1.0227979   28.64206
#> Nb_Comp_6  1.3048716   33.63927
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"


data(aze_compl)
summary(cv.plsR(y~.,data=aze_compl,nt=10,K=8,modele="pls",verbose=FALSE),
MClassed=TRUE,verbose=FALSE)
#> [[1]]
#>                 AIC MissClassed CV_MissClassed       Q2cum_Y LimQ2_Y       Q2_Y
#> Nb_Comp_0  154.6179          49             NA            NA      NA         NA
#> Nb_Comp_1  126.4083          27             45    -0.1267597  0.0975 -0.1267597
#> Nb_Comp_2  119.3375          25             50    -0.7883936  0.0975 -0.5872005
#> Nb_Comp_3  114.2313          27             44    -2.4260707  0.0975 -0.9157252
#> Nb_Comp_4  112.3463          23             48    -6.6030021  0.0975 -1.2191609
#> Nb_Comp_5  113.2362          22             52   -16.9740350  0.0975 -1.3640708
#> Nb_Comp_6  114.7620          21             53   -42.5282283  0.0975 -1.4217283
#> Nb_Comp_7  116.5264          20             52  -104.3561948  0.0975 -1.4204108
#> Nb_Comp_8  118.4601          20             52  -253.9485578  0.0975 -1.4198725
#> Nb_Comp_9  120.4452          19             52  -617.9161650  0.0975 -1.4276119
#> Nb_Comp_10 122.4395          19             53 -1507.0570241  0.0975 -1.4366095
#>             PRESS_Y    RSS_Y      R2_Y  AIC.std  DoF.dof sigmahat.dof   AIC.dof
#> Nb_Comp_0        NA 25.91346        NA 298.1344  1.00000    0.5015845 0.2540061
#> Nb_Comp_1  29.19824 19.38086 0.2520929 269.9248 22.55372    0.4848429 0.2883114
#> Nb_Comp_2  30.76131 17.76209 0.3145613 262.8540 27.31542    0.4781670 0.2908950
#> Nb_Comp_3  34.02728 16.58896 0.3598323 257.7478 30.52370    0.4719550 0.2902572
#> Nb_Comp_4  36.81358 15.98071 0.3833049 255.8628 34.00000    0.4744263 0.3008285
#> Nb_Comp_5  37.77952 15.81104 0.3898523 256.7527 34.00000    0.4719012 0.2976347
#> Nb_Comp_6  38.29004 15.73910 0.3926285 258.2785 34.00000    0.4708264 0.2962804
#> Nb_Comp_7  38.09509 15.70350 0.3940024 260.0429 33.71066    0.4693382 0.2937976
#> Nb_Comp_8  38.00046 15.69348 0.3943888 261.9766 34.00000    0.4701436 0.2954217
#> Nb_Comp_9  38.09768 15.69123 0.3944758 263.9617 33.87284    0.4696894 0.2945815
#> Nb_Comp_10 38.23340 15.69037 0.3945088 265.9560 34.00000    0.4700970 0.2953632
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
#> attr(,"computed_nt")
#> [1] 10
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRmodel"
summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,modele="pls",verbose=FALSE),
MClassed=TRUE,verbose=FALSE)
#> [[1]]
#>                 AIC MissClassed CV_MissClassed       Q2cum_Y LimQ2_Y       Q2_Y
#> Nb_Comp_0  154.6179          49             NA            NA      NA         NA
#> Nb_Comp_1  126.4083          27             52    -0.2605261  0.0975 -0.2605261
#> Nb_Comp_2  119.3375          25             49    -1.2082024  0.0975 -0.7518101
#> Nb_Comp_3  114.2313          27             48    -3.4138555  0.0975 -0.9988456
#> Nb_Comp_4  112.3463          23             52    -9.5433577  0.0975 -1.3886957
#> Nb_Comp_5  113.2362          22             48   -25.6668906  0.0975 -1.5292598
#> Nb_Comp_6  114.7620          21             48   -66.4979067  0.0975 -1.5311502
#> Nb_Comp_7  116.5264          20             49  -168.4398626  0.0975 -1.5102980
#> Nb_Comp_8  118.4601          20             48  -424.6365094  0.0975 -1.5120211
#> Nb_Comp_9  120.4452          19             49 -1069.4054719  0.0975 -1.5148347
#> Nb_Comp_10 122.4395          19             49 -2691.6606109  0.0975 -1.5155520
#>             PRESS_Y    RSS_Y      R2_Y  AIC.std  DoF.dof sigmahat.dof   AIC.dof
#> Nb_Comp_0        NA 25.91346        NA 298.1344  1.00000    0.5015845 0.2540061
#> Nb_Comp_1  32.66459 19.38086 0.2520929 269.9248 22.55372    0.4848429 0.2883114
#> Nb_Comp_2  33.95159 17.76209 0.3145613 262.8540 27.31542    0.4781670 0.2908950
#> Nb_Comp_3  35.50367 16.58896 0.3598323 257.7478 30.52370    0.4719550 0.2902572
#> Nb_Comp_4  39.62598 15.98071 0.3833049 255.8628 34.00000    0.4744263 0.3008285
#> Nb_Comp_5  40.41935 15.81104 0.3898523 256.7527 34.00000    0.4719012 0.2976347
#> Nb_Comp_6  40.02012 15.73910 0.3926285 258.2785 34.00000    0.4708264 0.2962804
#> Nb_Comp_7  39.50983 15.70350 0.3940024 260.0429 33.71066    0.4693382 0.2937976
#> Nb_Comp_8  39.44751 15.69348 0.3943888 261.9766 34.00000    0.4701436 0.2954217
#> Nb_Comp_9  39.46651 15.69123 0.3944758 263.9617 33.87284    0.4696894 0.2945815
#> Nb_Comp_10 39.47210 15.69037 0.3945088 265.9560 34.00000    0.4700970 0.2953632
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
summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,
modele="pls-glm-family",
family=gaussian(),verbose=FALSE),
MClassed=TRUE,verbose=FALSE)
#> [[1]]
#>                 AIC      BIC MissClassed CV_MissClassed Q2Chisqcum_Y  limQ2
#> Nb_Comp_0  154.6179 159.9067          49             NA           NA     NA
#> Nb_Comp_1  126.4083 134.3415          27             NA           NA 0.0975
#> Nb_Comp_2  119.2021 129.7796          28             NA           NA 0.0975
#> Nb_Comp_3  113.9553 127.1773          26             NA           NA 0.0975
#> Nb_Comp_4  112.4466 128.3130          25             NA           NA 0.0975
#> Nb_Comp_5  113.2280 131.7387          23             NA           NA 0.0975
#> Nb_Comp_6  114.7095 135.8646          21             NA           NA 0.0975
#> Nb_Comp_7  116.5144 140.3139          20             NA           NA 0.0975
#> Nb_Comp_8  118.4615 144.9054          20             NA           NA 0.0975
#> Nb_Comp_9  120.4453 149.5336          19             NA           NA 0.0975
#> Nb_Comp_10 122.4403 154.1729          19             NA           NA 0.0975
#>            Q2Chisq_Y PREChi2_Pearson_Y Chi2_Pearson_Y    RSS_Y      R2_Y
#> Nb_Comp_0         NA                NA       25.91346 25.91346        NA
#> Nb_Comp_1         NA                NA       19.38086 19.38086 0.2520929
#> Nb_Comp_2         NA                NA       17.73898 17.73898 0.3154532
#> Nb_Comp_3         NA                NA       16.54501 16.54501 0.3615285
#> Nb_Comp_4         NA                NA       15.99613 15.99613 0.3827095
#> Nb_Comp_5         NA                NA       15.80978 15.80978 0.3899009
#> Nb_Comp_6         NA                NA       15.73116 15.73116 0.3929346
#> Nb_Comp_7         NA                NA       15.70168 15.70168 0.3940726
#> Nb_Comp_8         NA                NA       15.69369 15.69369 0.3943807
#> Nb_Comp_9         NA                NA       15.69125 15.69125 0.3944749
#> Nb_Comp_10        NA                NA       15.69049 15.69049 0.3945043
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,
modele="pls-glm-logistic",
verbose=FALSE),MClassed=TRUE,verbose=FALSE)
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
summary(cv.plsRglm(y~.,data=aze_compl,nt=10,K=8,
modele="pls-glm-family",
family=binomial(),verbose=FALSE),
MClassed=TRUE,verbose=FALSE)
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


if(require(chemometrics)){
data(hyptis)
hyptis
yhyptis <- factor(hyptis$Group,ordered=TRUE)
Xhyptis <- as.data.frame(hyptis[,c(1:6)])
options(contrasts = c("contr.treatment", "contr.poly"))
modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
modpls2$Coeffsmodel_vals
modpls2$InfCrit
modpls2$Coeffs
modpls2$std.coeffs

table(yhyptis,predict(modpls2$FinalModel,type="class"))

modpls3 <- PLS_glm(yhyptis[-c(1,2,3)],Xhyptis[-c(1,2,3),],3,modele="pls-glm-polr",
dataPredictY=Xhyptis[c(1,2,3),],verbose=FALSE)

summary(cv.plsRglm(factor(Group,ordered=TRUE)~.,data=hyptis[,-c(7,8)],nt=4,K=10,
random=TRUE,modele="pls-glm-polr",keepcoeffs=TRUE,verbose=FALSE),
MClassed=TRUE,verbose=FALSE)
}
#> ____************************************************____
#> 
#> Model: pls-glm-polr 
#> Method: logistic 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> [[1]]
#>                AIC     BIC MissClassed CV_MissClassed Q2Chisqcum_Y  limQ2
#> Nb_Comp_0 86.87461 91.0782          20             NA           NA     NA
#> Nb_Comp_1 72.73191 78.3367          13             15    -5.092085 0.0975
#>           Q2Chisq_Y PREChi2_Pearson_Y Chi2_Pearson_Y
#> Nb_Comp_0        NA                NA       60.00011
#> Nb_Comp_1 -5.092085          365.5258       30.46894
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
# }
```
