# mkuhn, 2024-12-09
# testing inference methods


test_that("kappa homogeneity significance test for independent samples", {
  # Cohen's kappa on three studies
  k2_studies <- lapply(agreem_binary, kappa2)
  expect_named(k2_studies, c("study1", "study2", "study3"))
  
  # cf Solutions in Appendix C to problem 18.3
  # in Fleiss, Statistical Methods.. 3rd ed (2003)
  expect_identical(purrr::map_dbl(k2_studies, "subjects"),
                   expected = c(20, 20, 30),
                   ignore_attr = "names")
  expect_equal(purrr::map_dbl(k2_studies, "value"),
               expected = c(.39, .48, .35),
               tolerance = .025, ignore_attr = "names")
  expect_equal(purrr::map_dbl(k2_studies, "se"),
               expected = c(.21, .25, .17),
               tolerance = .025, ignore_attr = "names")
  
  # homogeneity test for kappa
  k2hom <- kappa_test(kappas = k2_studies,
                      val = "value", se = "se")
  
  expect_s3_class(k2hom, class = "htest")
  expect_named(k2hom, expected = c("method", "data.name", "estimate",
                                   "statistic", "p.value", 
                                   "alternative", "parameter", "conf.int"))
  expect_equal(k2hom$estimate, expected = c(`overall kappa`=0.39),
               tolerance = .01)
  expect_equal(k2hom$statistic, expected = c(`X^2`=0.18),
               tolerance = .1)
  expect_equal(k2hom$conf.int, expected = c(.16, .62),
               tolerance = .1, ignore_attr = "conf.level")
})


test_that("kappa homogeneity test for correlated samples", {
  
  set.seed(2024-12-09)
  k2hom <- kappa_test_corr(ratings = depression,
                           grpIdx = list(c(1, 2), c(1, 3)),
                           kappaF = kappa2,
                           kappaF_args = list(ratingScale = c("neg", "pos")),
                           B = 2500)
  
  # see Vanbelle (2008), section 5.2
  expect_equal(as.numeric(k2hom$statistic), 2.19, tolerance = .05)
  expect_equal(k2hom$p.value, 0.14, tolerance = .05)
})