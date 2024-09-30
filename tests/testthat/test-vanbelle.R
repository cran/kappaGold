# test agreement estimates between groups of raters due to Vanbelle

test_that("kappa vanbelle", {

  # taken from vanbelle, albert paper, see Table 5, first row "proposed method"
  #+they seem to have used linear weights
  kvb_sct <- kappam_vanbelle(ratingsGr1 = SC_test[, startsWith(colnames(SC_test), "S")],
                             ratingsGr2 = SC_test[, startsWith(colnames(SC_test), "E")],
                             ratingScale = -2:2,
                             weights = "linear")
  
  expect_named(kvb_sct,
               expected = c("method", "subjects", "raters", "categories",
                            "weights", "value0", "value", 
                            "se", "conf.level", "ci.lo", "ci.hi", "ci.width"))
  expect_equal(kvb_sct$value, expected = 0.72, tolerance = 1e-2)
  expect_equal(kvb_sct$se, expected = 0.049, tolerance = 1e-2)
})
