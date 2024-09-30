

# data -----------------------------------------------------------------

# cf dataset anxiety from irr-package
anxiety <- data.frame(rater1 = c(3L, 3L, 3L, 4L, 5L, 5L, 2L, 3L, 5L, 2L, 2L, 6L, 1L, 5L, 2L, 2L, 1L, 2L, 4L, 3L),
                      rater2 = c(3L, 6L, 4L, 6L, 2L, 4L, 2L, 4L, 3L, 3L, 2L, 3L, 3L, 3L, 2L, 2L, 1L, 3L, 3L, 4L),
                      rater3 = c(2L, 1L, 4L, 4L, 3L, 2L, 1L, 6L, 1L, 1L, 1L, 2L, 3L, 3L, 1L, 1L, 3L, 3L, 2L, 2L))

# example data for 3-level categorical UICC rating
# gold standard in first column with no NA
# ordinary raters have many NAs
uicc <- structure(list(
  GS = c("I", "II", "II", "III", "I", "III", "I", 
         "I", "III", "II", "III", "II", "I", "II", "II", "III", "III", 
         "I", "I", "I", "III", "I", "III", "III", "III", "II", "III", 
         "III", "I", "I", "I", "I", "III", "I", "I", "I", "I", "II", "II", 
         "I", "I", "III", "III", "II", "I", "I", "III", "II", "III", "I", 
         "I"),
  rater1 = c("III", "II", "II", "II", "III", "III", NA, "I", 
             NA, "III", NA, "I", "I", "III", NA, "III", "III", "I", "I", "I", 
             NA, "I", NA, NA, NA, NA, NA, "I", "I", "II", NA, "I", "II", NA, 
             "II", "III", "I", "II", "III", "II", NA, "II", NA, NA, "I", "I", 
             "III", NA, "III", NA, "I"),
  rater2 = c("II", "II", "II", "II", 
             NA, "III", "II", NA, "II", "II", "III", "I", "II", NA, NA, "III", 
             "I", "I", "I", "I", NA, "III", "III", "II", "II", NA, "I", "II", 
             "I", "II", "I", "I", "I", NA, NA, "II", "I", "II", NA, "II", 
             "I", NA, "II", "II", "II", "I", NA, "III", NA, "II", "II"),
  rater3 = c(NA, 
             "II", "II", "I", "III", "III", "I", "III", "I", "III", "III", 
             "I", "II", "III", "III", "III", "I", "I", "II", "II", "III", 
             "III", "III", "II", "III", "III", "I", NA, "I", "I", NA, "I", 
             "I", NA, "III", "III", "II", NA, "III", "III", "III", "III", 
             "III", "III", NA, NA, "II", "III", "III", "III", "II"),
  rater4 = c(NA, 
             NA, NA, NA, NA, NA, "I", NA, NA, "II", NA, "I", NA, NA, NA, "I", 
             NA, NA, "I", "I", "II", "III", "III", "II", "II", NA, "I", "II", 
             "I", "I", "I", "I", "II", NA, NA, "II", "I", "II", "III", "I", 
             "II", "III", "II", "III", "II", "I", "III", "II", "III", NA, 
             "I"),
  rater5 = c("I", "II", "I", "I", NA, NA, "I", NA, NA, "III", 
             NA, "I", "II", "III", NA, "I", "I", "I", "I", "II", NA, "III", 
             NA, "III", NA, "II", "III", "II", "I", "I", NA, "I", "I", "I", 
             "I", NA, "II", "II", "II", "II", "II", "III", "II", "II", "I", 
             "II", "II", "II", "III", "II", "II"),
  rater6 = c(NA, "II", "II", 
             "II", NA, "III", "II", "I", "II", "III", "III", "I", "I", "III", 
             "III", NA, "I", "I", "II", "I", NA, "III", "III", NA, NA, "II", 
             "I", "II", "I", "II", "III", "II", "II", "I", NA, "II", "II", 
             "II", "III", "II", "II", "III", NA, "II", "I", "I", "III", "II", 
             "III", NA, "II"),
  rater7 = c(NA, "II", "II", "II", "III", "III", 
             "III", "I", "I", "II", "II", "II", "II", "III", "III", "III", 
             NA, "II", "I", "II", NA, "III", "III", "II", "II", "III", "III", 
             NA, "I", "II", "III", "I", "II", NA, "III", "III", "I", "II", 
             "III", "III", "III", "III", "III", "III", "I", "I", "II", "III", 
             "III", NA, "I"),
  rater8 = c("I", NA, "I", NA, "II", "III", "III", 
             "III", "I", "III", "III", "I", NA, "III", "III", NA, "III", NA, 
             NA, NA, "III", "III", "III", NA, NA, "III", "III", NA, NA, NA, 
             NA, NA, "I", "III", "III", "III", NA, NA, "III", "III", NA, "III", 
             "III", "III", NA, "III", "III", "III", "III", "III", NA),
  rater9 = c(NA, 
             NA, "I", NA, NA, NA, "II", NA, NA, "III", "III", NA, NA, "III", 
             "II", "III", "I", "II", "I", "III", "III", "III", "III", "II", 
             "II", "III", "III", "III", "I", "I", NA, NA, "I", "I", "III", 
             NA, "II", "I", "III", "II", "III", "III", NA, NA, "II", NA, "III", 
             NA, "III", "II", "II"),
  rater10 = c("II", "I", "I", "II", "III", 
              "III", "III", "II", "I", NA, "I", NA, "I", NA, NA, "II", NA, 
              "I", "I", NA, NA, "III", NA, "II", "I", NA, "I", "I", "II", "I", 
              NA, "I", NA, "III", NA, "III", "I", NA, "III", "II", NA, "III", 
              NA, NA, "I", NA, "III", "I", NA, "III", "I"),
  rater11 = c("II", 
              "I", "II", "I", "I", "III", NA, NA, "III", NA, NA, "I", "I", 
              "III", "II", "III", NA, "III", NA, "I", "III", NA, "III", "II", 
              "I", "II", "I", "I", "I", "I", NA, "I", NA, NA, NA, "III", "I", 
              "II", "III", "III", "II", "III", NA, NA, "II", "I", NA, NA, NA, 
              NA, NA),
  rater12 = c(NA, NA, NA, "III", "III", "III", "III", 
              NA, "III", "III", "III", NA, "I", "III", NA, "III", NA, "I", 
              NA, "III", "III", "III", "III", "III", "III", "II", "III", "II", 
              "I", "I", NA, NA, NA, "III", NA, "III", "II", "III", "III", NA, 
              NA, "III", "III", NA, NA, "II", "III", NA, "III", NA, "II"), 
  rater13 = c("I", "I", "I", "I", "III", "III", NA, NA, "II", 
              "III", "I", "I", NA, NA, "III", "II", "I", "III", "I", "I", 
              "III", NA, NA, NA, NA, NA, NA, NA, NA, "I", NA, "I", "II", 
              "II", NA, "III", "I", NA, "III", NA, "III", "II", "III", 
              NA, "I", "I", "III", "III", NA, "II", "I"),
  rater14 = c("II", "II", "II", "II", "II", "II", "II", "II", "II", NA, "II", 
              "II", "II", "II", NA, "II", "I", "II", "II", "II", "II", 
              NA, NA, "II", "II", "II", "I", "I", "I", "II", "I", "I", 
              "II", "I", "I", "II", "II", "II", "II", "I", "II", "III", 
              "II", "I", "II", NA, "II", "II", NA, NA, "II"), 
  rater15 = c("II", "II", "II", "II", NA, NA, "I", "II", "II", 
              "II", "III", "I", "II", "III", "III", "III", "I", "II", "I", 
              "III", NA, "III", "III", "III", "II", "II", "I", "I", "I", 
              "I", "I", "I", "II", "III", "I", NA, NA, "II", "III", "II", 
              "II", "III", "II", "II", "I", "I", "II", "II", NA, "II", 
              "I")),
  row.names = c("sj1", "sj2", "sj3", "sj4", 
                "sj5", "sj6", "sj7", "sj8", "sj9", "sj10", "sj11", 
                "sj12", "sj13", "sj14", "sj15", "sj16", "sj17", "sj18", 
                "sj19", "sj20", "sj21", "sj22", "sj23", "sj24", "sj25", 
                "sj26", "sj27", "sj28", "sj29", "sj30", "sj31", "sj32", 
                "sj33", "sj34", "sj35", "sj36", "sj37", "sj38", "sj39", 
                "sj40", "sj41", "sj42", "sj43", "sj44", "sj45", "sj46", 
                "sj47", "sj48", "sj49", "sj50", "sj51"),
  class = "data.frame")


# cf Fleiss, Cohen, Everitt "large sample std err of kappa", 1969
# Table1 (and also Table2)
FCEdat <- tibble(n = as.integer(200 * c(.53, .11, .01,
                                        .05, .14, .06,
                                        .02, .05, .03)),
                  R1 = rep.int(c("I", "II", "III"), times = 3),
                  R2 = rep(c("I", "II", "III"), each = 3)) |> 
  tidyr::uncount(weights = n)


test_that("ratings data sets", {
  expect_true(is.matrix(diagnoses))
  expect_type(diagnoses, type = "character")
  expect_s3_class(anxiety, class = "data.frame")
  expect_s3_class(uicc, class = "data.frame")
  expect_s3_class(FCEdat, class = "data.frame")
  expect_identical(NROW(FCEdat), expected = 200L)
  expect_identical(NCOL(FCEdat), expected = 2L)
  
  # 1st column is gold standard
  expect_identical(names(uicc)[[1L]], expected = "GS")
  expect_identical(sum(is.na(uicc$GS)), expected = 0L)
})



# kappa -------------------------------------------------------------------


test_that("kappa2 implementation", {
  # select two columns manually
  k2_diag23 <- kappa2(diagnoses[, 2:3])
  expect_type(k2_diag23, type = "list")
  expect_named(k2_diag23,
               expected = c("method", "subjects", "raters", "categories",
                            "robust", "agreem", "value", "SE"))
  expect_equal(k2_diag23$value,
               expected = irr::kappa2(diagnoses[,2:3])$value)
  
  # all sets of two columns
  purrr::walk(.x = utils::combn(NCOL(diagnoses), 2L, simplify = FALSE),
              .f = ~ expect_equal(kappa2(diagnoses[, .x])$value,
                                 expected = irr::kappa2(diagnoses[, .x])$value))
  
  # cf Fleiss, Cohen, Everitt "large sample std err of kappa", 1969
  # unweighted kappa (after Table2)
  k2_FCE <- kappa2(FCEdat)
  expect_equal(k2_FCE$value, expected = .429, tolerance = 1e-3)
  expect_equal(k2_FCE$SE^2, expected = .002885, tolerance = 1e-4)
})


test_that("kappam_fleiss", {
  diag_kf <- kappam_fleiss(ratings = diagnoses, variant = "fleiss")
  # actually, Conger does not fit here to the diagnoses dataset
  #+as we have no select set of six fixed raters!
  diag_kc <- kappam_fleiss(ratings = diagnoses, variant = "conger")
  
  K_NAMES <- c("method", "subjects", "raters", "categories", "ratings", "balanced",
               "agreem", "value0", "value", "se_j", "stat.name", "statistic", "p.value")
  expect_type(diag_kf, "list")
  expect_named(diag_kf, expected = c(K_NAMES, "se0", "statistic0", "p.value0"))
  expect_equal(diag_kf$value0, expected = 0.43024, tolerance = 1e-4) # from irr::kappam.fleiss
  
  expect_type(diag_kc, "list")
  expect_named(diag_kc, expected = K_NAMES)
  # diagnoses with un-permuated rater1-rater6: 0.44181
  expect_equal(diag_kc$value0, expected = 0.43112, tolerance = 1e-4) # from irr::kappam.fleiss
  
  # Fleiss, "Stat. methods for rates and prop.", 3rd ed, section 18.3, table 18.8:
  # five ratings on each of ten subjects into one of three categories
  # we assume that we have identical raters
  fleissmEx18.8 <- list(
    nSj = 10,
    nRt = 5)
  
  fleissmEx18.8$tab_cnt_sj <- matrix(c(1, 4, 0,
                                   2, 0, 3,
                                   0, 0, 5,
                                   4, 0, 1,
                                   3, 0, 2,
                                   1, 4, 0,
                                   5, 0, 0,
                                   0, 4, 1,
                                   1, 0, 4,
                                   3, 0, 2),
                                 nrow = fleissmEx18.8$nSj, byrow = TRUE,
                                 dimnames = list(subj = paste0("S", seq_len(fleissmEx18.8$nSj)),
                                                 cat = LETTERS[1:3]))
  fleissmEx18.8$prob_tab <- colSums(fleissmEx18.8$tab_cnt_sj) / sum(fleissmEx18.8$tab_cnt_sj)
  
  
  fleissmEx18.8$rat <- fleissmEx18.8$tab_cnt_sj |> 
    tibble::as_tibble(rownames = "sjID") |> 
    tidyr::pivot_longer(cols = c(A, B, C), names_to = "Category", values_to = "Cnt") |> 
    tidyr::uncount(weights = Cnt) |> 
    # optionally, permute rows within sjID to have exchangeable raters
    #+with no permutation it is ordered and first rater tends to give A etc
    dplyr::mutate(rID = dplyr::row_number(), .by = sjID) |> 
    tidyr::pivot_wider(id_cols = sjID,
                       names_from = rID, names_prefix = "R",
                       values_from = Category) |> 
    tibble::column_to_rownames(var = "sjID")
  
  fleissmEx18.8$kappa <- kappam_fleiss(ratings = fleissmEx18.8$rat, detail = TRUE)
  # 
  with(fleissmEx18.8, {
    expect_identical(kappa$categories, 3L)
    expect_equal(kappa$value0, 0.42, tolerance = .01) #value given by Fleiss
    expect_equal(kappa$value0,
                 1 - (nSj * nRt^2 - sum(tab_cnt_sj^2)) / (nSj * nRt * (nRt-1) * crossprod(prob_tab, 1-prob_tab)[1]))
    expect_equal(kappa$se0, 0.072, tolerance = .01)
    expect_equal(kappa$statistic0, 5.83, tolerance = .01)
    expect_identical(colnames(kappa$detail), expected = c("kappa_j", "se0_j", "z_j", "p.value_j"))
    expect_equal(kappa$detail[, "kappa_j"], c(A=.29, B=.67, C=.35), tolerance = .01)
    expect_equal(kappa$detail[, "se0_j"], rep(0.10, 3), ignore_attr = TRUE)
  })
  
  # cf Fleiss (2003), sec 18.3, table 18.7
  fleissmEx18.7 <- matrix(c(2, 2,
                            2, 0,
                            3, 2,
                            4, 3,
                            3, 3,
                            4, 1,
                            3, 0,
                            5, 0,
                            2, 0,
                            4, 4,
                            5, 5,
                            3, 3,
                            4, 4,
                            4, 3,
                            2, 0,
                            2, 2,
                            3, 1,
                            2, 1,
                            4, 1,
                            5, 4,
                            3, 2,
                            4, 0,
                            3, 0,
                            3, 3,
                            2, 2), ncol = 2, byrow = TRUE,
                          dimnames = list(subj = paste0("S", seq_len(25)),
                                          cnt = c("m_i", "x_i"))) |> 
    tibble::as_tibble(rownames = "sjID") |> 
    dplyr::rename(A = x_i) |> 
    dplyr::mutate(B = m_i - A, m_i = NULL) |> 
    tidyr::pivot_longer(cols = c(A, B), names_to = "Category", values_to = "Cnt") |> 
    tidyr::uncount(weights = Cnt) |> 
    # optionally, permute rows within sjID to have exchangeable raters
    #+with no permutation it is ordered and first rater tends to give A etc
    dplyr::mutate(rID = dplyr::row_number(), .by = sjID) |> 
    tidyr::pivot_wider(id_cols = sjID,
                       names_from = rID, names_prefix = "R",
                       values_from = Category) |> 
    tibble::column_to_rownames(var = "sjID")
  
  kfEx18.7 <- kappam_fleiss(ratings = fleissmEx18.7)
  kfEx18.7_cong <- kappam_fleiss(ratings = fleissmEx18.7, variant = "cong")
  kfEx18.7_rob <- kappam_fleiss(ratings = fleissmEx18.7, variant = "rob")
  
  expect_identical(kfEx18.7$subjects, 25L)
  expect_identical(kfEx18.7$raters, 5L)
  expect_identical(kfEx18.7$ratings, 81L)
  expect_identical(kfEx18.7$categories, 2L)
  expect_false(kfEx18.7$balanced)
  expect_equal(kfEx18.7$value0, 0.54, tolerance = .05) #cf sect 18.3, p. 613
  expect_equal(kfEx18.7$se0, 0.103, tolerance = .01)
  expect_equal(kfEx18.7$statistic0, 5.24, tolerance = .05)
  
  # conger reduces the chance agreement a little bit, hence kappa increases
  expect_gt(kfEx18.7_cong$value0, kfEx18.7$value0)
  # robust assumes smallest possible chance agreement, hence kappa increases
  expect_gt(kfEx18.7_rob$value0, kfEx18.7$value0)
})


test_that("kappam_gold", {
  
  # take first rater as reference
  kg_anx <- kappam_gold(ratings = anxiety)
  
  expect_identical(kg_anx$subjects, expected = NROW(anxiety))
  expect_identical(kg_anx$raters, expected = NCOL(anxiety)-1)
  # kappam_gold as average of pairwise Cohen's kappa
  expect_equal(kg_anx$value0,
               expected = mean(purrr::map_dbl(2L:NCOL(anxiety), ~ kappa2(ratings = anxiety[, c(1L,.x)])$value)))
  
  
  kg_diag <- kappam_gold(ratings = diagnoses)
  expect_identical(kg_diag$subjects, expected = NROW(diagnoses))
  expect_identical(kg_diag$raters, expected = NCOL(diagnoses)-1)
  # kappam_gold as average of pairwise Cohen's kappa
  expect_equal(kg_diag$value0,
               expected = mean(purrr::map_dbl(2:NCOL(diagnoses),
                                              ~ kappa2(ratings = diagnoses[, c(1L, .x)])$value)))
  
  # copes with NAs
  expect_no_error({kg_uicc <- kappam_gold(ratings = uicc)})
  expect_type(kg_uicc, type = "list")
  expect_named(kg_uicc, expected = c("method",
                                     "subjects", "raters", "categories",
                                     "agreem",
                                     "value0", "value",
                                     "se_j", "conf.level", "ci.lo", "ci.hi", "ci.width"))
  expect_identical(kg_uicc$subjects, NROW(uicc))
  expect_identical(kg_uicc$raters, NCOL(uicc)-1)
  # pairwise Cohen's kappa with gold standard
  k2_uicc <- purrr::map(.x = 2L:NCOL(uicc), .f = ~ kappa2(ratings = uicc[, c(1L, .x)]))
  expect_identical(purrr::map_dbl(k2_uicc, "subjects"), expected = colSums(!is.na(uicc[,-1L])),
                   ignore_attr = "names")
  expect_equal(kg_uicc$value0,
               expected = weighted.mean(x = purrr::map_dbl(k2_uicc, "value"),
                                        w = purrr::map_dbl(k2_uicc, "subjects")))
})


test_that("kappa homogeneity test", {
  diagnoses1 <- diagnoses[1:10,]
  diagnoses2 <- diagnoses[-(1:10),]
  
  
  # homogeneity test
  kt <- kappa_test(kappas = list(kappa1=kappam_fleiss(diagnoses1),
                                 kappa2=kappam_fleiss(diagnoses2)),
                   val = "value0", se = "se0")
  
  expect_named(kt, expected = c("method", "data.name", "estimate", "statistic", "p.value", 
                                "alternative", "parameter", "conf.int"))
  expect_gt(kt$p.value, expected = 0.1)
})


test_that("kappa simulation", {
  simK <- simulKappa(nRater = 8, cats = 3, nSubj = 11, mcSim = 9,
             # assumed prob for classification by raters
             probs = matrix(c(.6, .2, .1, # subjects of cat 1
                              .3, .4, .3, # subjects of cat 2
                              .1, .4, .5  # subjects of cat 3
             ), nrow = 3, byrow = TRUE))
  
  expect_type(simK, "list")
  expect_s3_class(simK, class = "data.frame")
  
  expect_named(simK,
               expected = c("nRater", "nSubjTotal", "kappam_gold", "ci_lo", "ci_hi", "ci_halfwidth"))
  expect_identical(simK$nRater, expected = rep.int(8, times = 9))
  
})