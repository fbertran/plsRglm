# Summary method for plsRglm models

This function provides a summary method for the class
`"cv.plsRglmmodel"`

## Usage

``` r
# S3 method for class 'cv.plsRglmmodel'
summary(object, ...)
```

## Arguments

- object:

  an object of the class `"cv.plsRglmmodel"`

- ...:

  further arguments to be passed to or from methods.

## Value

An object of class `"summary.cv.plsRmodel"` if `model` is missing or
`model="pls"`. Otherwise an object of class `"summary.cv.plsRglmmodel"`.

## References

Nicolas Meyer, Myriam Maumy-Bertrand et Frédéric Bertrand (2010).
Comparing the linear and the logistic PLS regression with qualitative
predictors: application to allelotyping data. *Journal de la Societe
Francaise de Statistique*, 151(2), pages 1-18.
<https://www.numdam.org/item/JSFS_2010__151_2_1_0/>

## See also

[`summary`](https://rdrr.io/r/base/summary.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
summary(cv.plsRglm(Y~.,data=Cornell,nt=10,NK=1,
modele="pls-glm-family",family=gaussian(), verbose=FALSE))
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
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
#> 
#> NK: 1
#> [[1]]
#>                AIC      BIC Q2Chisqcum_Y  limQ2   Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0 82.01205 82.98186           NA     NA          NA                NA
#> Nb_Comp_1 53.15173 54.60645   0.89977127 0.0975   0.8997713          46.88666
#> Nb_Comp_2 31.46903 33.40866   0.95110461 0.0975   0.5121619          17.43654
#> Nb_Comp_3 31.54404 33.96857   0.85210738 0.0975  -2.0246742          15.02305
#> Nb_Comp_4 33.20141 36.11085  -0.03500905 0.0975  -5.9983822          29.60801
#> Nb_Comp_5 33.25554 36.64989 -11.79853168 0.0975 -11.3656230          50.84259
#> Nb_Comp_6 35.25533 39.13459           NA 0.0975          NA                NA
#>           Chi2_Pearson_Y      RSS_Y      R2_Y
#> Nb_Comp_0     467.796667 467.796667        NA
#> Nb_Comp_1      35.742486  35.742486 0.9235940
#> Nb_Comp_2       4.966831   4.966831 0.9893825
#> Nb_Comp_3       4.230693   4.230693 0.9909561
#> Nb_Comp_4       4.111608   4.111608 0.9912107
#> Nb_Comp_5       3.496135   3.496135 0.9925264
#> Nb_Comp_6       3.496074   3.496074 0.9925265
#> 
#> attr(,"class")
#> [1] "summary.cv.plsRglmmodel"
rm(list=c("XCornell","yCornell","bbb"))
#> Warning: object 'bbb' not found
```
