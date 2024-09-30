
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kappaGold

<!-- badges: start -->
<!-- badges: end -->

The R-package `kappaGold` estimates agreement of a group of raters with
a gold standard. It builds on the idea of Conger that the multi-rater
kappa due to Light (1971) is actually a mean of *all* pairwise Cohen’s
kappas. In the situation of a gold standard, we only consider the
pairwise Cohen’s kappas of each rater with that gold standard.

The implementation of this measure of agreement with a gold standard is
found in the function `kappam_gold`. This function expects a matrix of
ratings with observations in the row and raters in the columns. The
rater in the 1<sup>st</sup> column is taken to be the gold standard. The
delete-1 jackknife method is used to get an estimate of bias and
standard error.

## Example

In medicine, staging is the process of assessing the extent to which a
tumour has grown. Staging affets treatment choice, for instance, if
radiation is used or not. Pathological assessment is typically the
gold-standard while non-invasive imaging allows for easier and earlier
tumour staging by radiologists. Inspired by the
[OCUM-trial](https://link.springer.com/article/10.1245/s10434-019-07696-y)
on colorectal tumour staging the data set `stagingData` carries the
fictitious staging of 21 colorectal tumour patients by a pathologist
based on a histological sample (gold standard) and 5 different
radiologists. The agreement of the radiologists (columns 2 to 6) with
the pathological staging as gold standard can be estimated by
`kappam_gold`:

``` r
library("kappaGold")

# 1st column corresponds to gold-standard
kappam_gold(kappaGold::stagingData)
#> $method
#> [1] "Averaged Cohen's Kappa with gold standard"
#> 
#> $subjects
#> [1] 21
#> 
#> $raters
#> [1] 5
#> 
#> $categories
#> [1] 3
#> 
#> $agreem
#> [1] 0.60952
#> 
#> $value0
#> [1] 0.41429
#> 
#> $value
#> [1] 0.42552
#> 
#> $se_j
#> [1] 0.074303
#> 
#> $conf.level
#> [1] 0.95
#> 
#> $ci.lo
#> [1] 0.27989
#> 
#> $ci.hi
#> [1] 0.57115
#> 
#> $ci.width
#> [1] 0.29126
```

Entry `agreem` is the mean pairwise agreement between the raters (to be
evaluated) and the gold standard rating. The entry `value0` shows the
mean of all pairwise Cohen’s kappa between the raters and the gold
standard. Delete-1 jackknife gives an estimate for bias and standard
error. These quantities are used to get the bias-corrected estimate
`value` which can be used as point estimate and a 95% confidence
interval.

## Installation

Package `kappaGold` was recently released on CRAN. For installation,
simply issue `install.packages("kappaGold")` in your R-session.

The development of the R-package `kappaGold` is going on at
[Gitlab](https://gitlab.com/imb-dev/kappa_gold).

With the help of the `remotes`-package you can install the development
version of package `kappaGold` via:

``` r
remotes::install_gitlab("imb-dev/kappa_gold@develop")
```
