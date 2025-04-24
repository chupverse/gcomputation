survivalNET: an R Package for Flexible Relative and Net Survival 
================

## Description

The R package ‘survivalNET’ contains a variety of functions to estimate relative or net
survival models. S3 methods are included to evaluate the predictive capacities, as
well as predictiions from new observations.

## Basic Usage

``` r
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The parametric estimation of the relative survival model (weibull distribution)

sNET <- survivalNET(Surv(time, event) ~ stade + delay + sex, data = dataK, ratetable=fr.ratetable,
                     age="age", sex="sexchara", year="year", dist="weibull",
                     strata=NULL, weights=NULL)
 
sNET
#>              coef exp(coef) se(coef)       z      p
#> log sigma  6.8098  906.7017   0.1854 36.7295 0.0000
#> log nu    -0.0551    0.9464   0.0423 -1.3020 0.1929
#> stade      0.4996    1.6480   0.1108  4.5069 0.0000
#> delay      0.3983    1.4893   0.1089  3.6587 0.0003
#> sex       -0.6864    0.5034   0.1110 -6.1840 0.0000

#> Likelihood ratio test=70.5213 on 3 df, p=0
#> n=1000, number of events=369
```


## Calibration plot of the model at 2 years

``` r
plot(sNET, n.groups=3, pro.time=2*365.24, 
  ratetable=fr.ratetable, age="age", sex="sexchara", year="diagnum")
```

## Prediction for new individuals

``` r
# For a patient with 
predictions <- predict(sNET, newdata=data.frame( stade=0, delay=0, sex=2 ) )
plot(predictions$times/365.24, predictions$predictions, type="l",
  ylab="Predicted survival", xlab="Time in years")
```

## Installation

To install the latest release from CRAN:

``` r
install.packages("survivalNET")
```

To install the development version from GitHub:

``` r
remotes::install_github("chupverse/survivalNET")
```

## Reporting bugs

You can report any issues at this
[link](https://github.com/chupverse/survivalNET/issues).
