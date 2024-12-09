# test agreement estimates between groups of raters due to Vanbelle

test_that("kappa vanbelle", {

  SCt_idx_gr1 <- which(startsWith(colnames(SC_test), "S"))
  # taken from Vanbelle & Albert paper, see Table 5, first row "proposed method"
  #+they seem to have used linear weights
  kvb_sct <- kappam_vanbelle(ratings = SC_test,
                             refIdx = SCt_idx_gr1,
                             ratingScale = -2:2,
                             weights = "linear")
  
  expect_named(kvb_sct,
               expected = c("method", "subjects", "raters", "categories",
                            "weights", "value0", "value", 
                            "se", "conf.level", "ci.lo", "ci.hi", "ci.width"))
  expect_equal(kvb_sct$value, expected = 0.72, tolerance = 1e-2)
  expect_equal(kvb_sct$se, expected = 0.049, tolerance = 1e-2)

  # same results when inverting group1 <=> group2
  expect_equal(kappam_vanbelle(ratings = SC_test,
                                   refIdx = -SCt_idx_gr1,
                                   ratingScale = -2:2,
                                   weights = "linear"),
               expected = kvb_sct)
})
