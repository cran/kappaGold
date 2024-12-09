
# kappaGold 0.4.0
* new function `kappa_test_corr` add test for difference in kappa in correlated groups (same subjects) through bootstrap
* modify API:
    * rename entries to `se` (instead of `SE`)
    * `kappam_gold` has parameter `refIdx=` to specify the column of gold standard
    * `kappam_vanbelle` works with single matrix of ratings (`ratings=`) and `refIdx=` allows to specify one of the groups of raters

# kappaGold 0.3.2
* enhance description of package (upon request from CRAN)
* little code optimizations

# kappaGold 0.3.1
* add datasets
* drop dependency on obsolete `bootstrap`-package: last bits

# kappaGold 0.3.0
* `kappam.gold`
    * use z-quantile (std normal) for confidence interval (not t-quantile)
    * nRaters refers to number of evaluated raters (excluding gold standard)
* `kappam.fleiss` implemented
    * based on `irr::kappam.fleiss`
    * allow to set `ratingScale=`
    * jackknife method for bias-correction and SE
* `kappaTest` implements a homogeneity test for kappa's from independent groups of raters
* rename: `kappaInference` becomes `simulKappa` (more clear)
* drop dependency on obsolete `bootstrap`-package: use own implementation of jackknife


# kappaGold 0.2.3
* `kappa2`: acknowledge that we have SE only for classical kappa2 (not for robust variant)

# kappaGold 0.2.2
* `kappa2`: add std. error for unweighted kappa

# kappaGold 0.2.1
* Cohen's kappa `kappa2`
    * simplified calculation of random agreement
    * return more information as list
    * allow to set `ratingScale=`
    * allow for `robust=` variant

