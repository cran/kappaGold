TOL <- 1.5e-8

#' Cohen's kappa for nominal data
#' 
#' The data of ratings must be stored in a two column object,
#' each rater is a columns and the subjects are in the rows.
#' 
#' Every rating category is used and the levels are sorted.
#' Weighting is currently not implemented.
#' 
#' @examples
#' # 2 raters have assessed 4 subjects into categories "A", "B" or "C"
#' # organize ratings as two column matrix, one row per subject rated
#' m <- rbind(sj1 = c("A", "A"),
#'            sj2 = c("C", "B"),
#'            sj3 = c("B", "C"),
#'            sj4 = c("C", "C"))
#'            
#' # Cohen's kappa -----
#' kappa2(ratings = m)
#' 
#' # robust variant ---------
#' kappa2(ratings = m, robust = TRUE)
#' 
#' @param ratings matrix (dimension nx2), containing the ratings as subjects by raters
#' @param robust flag. Use robust estimate for random chance of agreement by Brennan-Prediger?
#' @param ratingScale Possible levels for the rating. Or `NULL`.
#' @returns list containing Cohen's kappa agreement measure (value) or `NULL` if no valid subjects
#' @seealso [irr::kappa2()]
#' @export
kappa2 <- function (ratings, robust = FALSE, ratingScale = NULL) {
  ratings <- as.matrix(ratings)
  # handle ratingScale first, before dropping incomplete observations
  if (is.null(ratingScale)) {
    # sort also drops NAs (due to na.last=NA)
    ratingScale <- sort(as.character(unique(c(ratings[,1L], ratings[,2L]))), na.last = NA)
  } else {
    if (anyDuplicated(ratingScale)) {
      stop("Duplicated entries in rating scale!", call. = FALSE)
    }#fi
    ratingScale <- as.character(ratingScale)
    
    if (!all(as.character(ratings[!is.na(ratings)]) %in% ratingScale)) {
      stop("Ratings ", paste0(unique(as.character(ratings)), collpase="*"),
           " do not match provided rating scale ", paste0(ratingScale, collapse="*"), call. = FALSE)
    }#fi
  }#esle
  
  nCat <- length(ratingScale)
  if (nCat <= 1L) {
    stop("Rating scale needs **at least two** levels!", call. = FALSE)
  }#fi
  
  # check assumptions
  if (NCOL(ratings) != 2L) {
    stop("Please provide **exactly two** raters!", call. = FALSE)
  }
  
  # remove subjects with NA ratings
  ratings <- stats::na.omit(ratings)
  
  nSj <- NROW(ratings) # nbr subjects
  if (nSj < 1) {
    return(NULL)
  }
  
  
  # cross-table of two raters
  rtab <- table(factor(ratings[, 1L], levels = ratingScale),
                factor(ratings[, 2L], levels = ratingScale))
  
  rtab.rs <- .rowSums(rtab, m = nCat, n = nCat)
  rtab.cs <- .colSums(rtab, m = nCat, n = nCat)
  
  # agreement as found on diagonal elements
  agreeP <- sum(diag(rtab)) / nSj 
  chanceP <- if (robust) {
    # uniform
    1/nCat
  } else {
    # marginal proportions
    crossprod(rtab.rs, rtab.cs)[1L] / nSj^2
  }
  
  # Cohen's kappa for 2 raters
  k2 <- agreeP - chanceP
  # normalization if numerator and denominator != 0
  #>agreeP = chanceP is always kappa=0 (no need to normalize)
  #>chanceP=1 => denominator 0
  if (abs(k2) > TOL && chanceP < 1) {
    k2 <- k2 / (1L - chanceP)
  }
  
  # standard error
  # cf Large Sample Standard Errors of Kappa ... (Fleiss, 1969)
  pRS <- rtab.rs / nSj
  pCS <- rtab.cs / nSj
  ptab0 <- rtab / nSj; diag(ptab0) <- 0
  dsTerm <- 0
  for (i in seq_along(ratingScale)) {
    for (j in seq_along(ratingScale)) {
      dsTerm <- dsTerm + ptab0[i,j] * (pCS[[i]] + pRS[[j]])^2
    }#rof
  }#rof
  
  varK2 <- 1/(nSj*(1-chanceP)^4) * (
    crossprod(diag(rtab)/nSj, ((1-chanceP) - (pCS + pRS) * (1-agreeP))^2)[1L] +
      (1-agreeP)^2 * dsTerm - (agreeP*chanceP - 2 * chanceP + agreeP)^2)
  
  # return
  list(
    method = "Cohen's Kappa for two Raters",
    subjects = nSj, raters = 2, categories = nCat,
    robust = robust,
    agreem = agreeP, value = k2,
    #XXX currently, SE only implemented for standard Cohen's Kappa
    SE = if (!robust && is.finite(varK2) && varK2 >= 0) sqrt(varK2) else NA_real_
  )
}


#' Fleiss' kappa for multiple nominal-scale raters
#' 
#' When multiple raters judge subjects on a nominal scale we can assess their agreement with Fleiss' kappa.
#' It is a generalization of Cohen's Kappa for two raters and there are different variants how to assess chance agreement.
#' 
#' Different **variants** of Fleiss' kappa are implemented.
#' By default (`variant="fleiss"`), the original Fleiss Kappa (1971) is calculated, together with an asymptotic standard error and test for kappa=0.
#' It assumes that the raters involved are not assumed to be the same (one-way ANOVA setting).
#' The marginal category proportions determine the chance agreement.
#' Setting `variant="conger"` gives the variant of Conger (1980) that reduces to Cohen's kappa when m=2 raters. 
#' It assumes identical raters for the different subjects (two-way ANOVA setting).
#' The chance agreement is based on the category proportions of each rater separately.
#' Typically, the Conger variant yields slightly higher values than Fleiss kappa.
#' `variant="robust"` assumes a chance agreement of two raters to be simply 1/q, where q is the number of categories (uniform model).
#' 
#' @examples
#' # 4 subjects were rated by 3 raters in categories "1", "2" or "3"
#' # organize ratings as matrix with subjects in rows and raters in columns
#' m <- matrix(c("3", "2", "3",
#'               "2", "2", "1",
#'               "1", "3", "1",
#'               "2", "2", "3"), ncol = 3, byrow = TRUE)
#' kappam_fleiss(m)
#'
#' # show category-wise kappas -----
#' kappam_fleiss(m, detail = TRUE)
#'
#' @param ratings matrix (subjects by raters), containing the ratings
#' @param variant Which variant of kappa? Default is Fleiss (1971). Other options are Conger (1980) or robust variant.
#' @param detail Should category-wise Kappas be computed? Only available for the Fleiss (1971) variant.
#' @param ratingScale Specify possible levels for the rating. Default `NULL` means to use all unique levels from the sample.
#' @returns list containing Fleiss's kappa agreement measure (value) or `NULL` if no subjects
#' @seealso [irr::kappam.fleiss()]
#' @export
kappam_fleiss <- function (ratings, variant = c("fleiss", "conger", "robust", "uniform"),
                           detail = FALSE, ratingScale = NULL) {
  variant <- match.arg(variant)
  
  ratings <- as.matrix(ratings)
  
  # handle ratingScale first, before dropping incomplete observations
  if (is.null(ratingScale)) {
    # sort also drops NAs (due to na.last=NA)
    ratingScale <- sort(unique(as.character(ratings)), na.last = NA)
  } else {
    if (anyDuplicated(ratingScale)) {
      stop("Duplicated entries in rating scale!", call. = FALSE)
    }
    ratingScale <- as.character(ratingScale)
    
    if (!all(as.character(ratings[!is.na(ratings)]) %in% ratingScale)) {
      stop("Ratings ", paste0(unique(as.character(ratings)), collpase="*"),
           " do not match provided rating scale ", paste0(ratingScale, collapse="*"), call. = FALSE)
    }
  }
  nCat <- length(ratingScale)
  if (nCat <= 1L) {
    stop("Rating scale needs **at least two** levels!", call. = FALSE)
  }
  
  # drop raters with only NA-ratings
  # using `drop = FALSE` important to stay matrix,
  #+even in edge cases with single col/row
  ratings <- ratings[, colSums(!is.na(ratings)) >= 1L, drop = FALSE]
  # drop subjects that are rated not at all or only once
  ratings <- ratings[rowSums(!is.na(ratings)) >= 2L, , drop = FALSE]
  
  nSj <- NROW(ratings)
  nr <- NCOL(ratings)
  
  
  if (nSj <= 0) {
    #warning("No subjects left!", call. = FALSE)
    return(NULL)
  }
  
  if (nr <= 0) {
    #warning("No raters left!", call. = FALSE)
    return(NULL)
  }
  
  method <- switch(variant, 
                   fleiss = "Fleiss' Kappa for m Raters",
                   conger = "Fleiss' Kappa for m Raters (Conger variant)",
                   uniform =,
                   robust = "Fleiss' Kappa for m Raters (robust variant)",
                   stop("Unknown variant of Fleiss kappa!", call. = FALSE)
  )
  
  
  # count nbr of rated categories per subject (nSj x nCat)
  tab_cnt_sj <- t(apply(ratings, 1,
                        FUN = function(ro) tabulate(factor(ro, levels = ratingScale), nbins = nCat)))
  rt_cnt_sj <- .rowSums(tab_cnt_sj, m = nSj, n = nCat)
  
  # build return value
  rval <- list(method = method,
               subjects = nSj, raters = nr, categories = nCat,
               ratings = sum(tab_cnt_sj), balanced = isTRUE(stats::sd(rt_cnt_sj) == 0))
  
  # @param idx index vector of rows to use
  kappam_fleiss0 <- function(idx) {
    idxl <- length(idx)
    nRatings <- sum(tab_cnt_sj[idx,])
    cat_cnt <- .colSums(tab_cnt_sj[idx,], m = idxl, n = nCat)
    cat_ssq <- crossprod(cat_cnt)[1L]
    # prop of concordant pairs per subject:
    P_sj <- 1/(rt_cnt_sj[idx] - 1) * (.rowSums(tab_cnt_sj[idx,]^2, m = idxl, n = nCat)/rt_cnt_sj[idx] - 1)
    agreeP <- stats::weighted.mean(P_sj, w = rt_cnt_sj[idx])
    #XXX weighting? # or simply mean(P_sj)
    #+P_sj are proportions, based on denominator rt_cnt_sj * (rt_cnt_sj-1)
    #+
    #+balanced data used for agreeP:
    #sum((colSums(sj_cnt_tab^2) - nr) / (nr * (nr - 1) * nSj)) #or
    #mean((colSums(sj_cnt_tab^2) - nr)) / (nr * (nr - 1))
    
    chanceP <- switch(variant, 
                      fleiss = cat_ssq / nRatings^2,
                      conger = local({
                        # counts of rated categories per rater (nr x nCat)
                        tab_cnt_rt <- t(apply(ratings[idx,], 2,
                                              FUN = function(co) tabulate(factor(co, levels = ratingScale), nbins = nCat)))
                        tab_prop_rt <- proportions(tab_cnt_rt, margin = 1)
                        # Conger (1980), p. 325 divides the correction term by (nr-1) (not by nr)
                        #+but (nr-1) leads to different results when I compare with his derivation in an balanced example
                        #+Formula 1:
                        #cat_ssq / nRatings^2 - sum(apply(tab_prop_rt, 2, stats::var)) / nr
                        
                        # For Conger, chanceP is "average proportion of raters who agreed on the classification of each subject"
                        #+He compares all pairs of raters, -- which is well-founded in his balanced setting as each rater rates all subjects
                        #+For the unbalanced setting, it might be that two raters do not share any subject.
                        #+But that might not be a problem if we assume exchangeable subjects (and raters)
                        #+Formula 2:
                        #pagr_mat <- tcrossprod(tab_prop_rt)
                        #sum(pagr_mat[lower.tri(pagr_mat)]) / sum(lower.tri(pagr_mat))
                        
                        #+Formulas 1 & 2 agree for balanced setting but are slightly different in unbalanced settings
                        #+I use Formula 2 because it seems more robust
                        pagr_mat <- tcrossprod(tab_prop_rt)
                        sum(pagr_mat[lower.tri(pagr_mat)]) / sum(lower.tri(pagr_mat))
                      }),
                      uniform =,
                      robust = {
                        # prop. of two raters being concordant
                        1 / nCat # XXX is it that simple?
                        # XXX could use observed nbr of ratings per subject (in particular when number of raters varies per subj)
                        #+and use same procedure as for P_sj
                        #+only using expected r_ij under 1/q assumption.
                        #+In my examples, it was close (but little lower) than 1/q
                      },
                      stop("Unknown variant of Fleiss kappa!", call. = FALSE)
    )
    
    rval0 <- (agreeP - chanceP) / (1 - chanceP)
    attr(rval0, "agreeP") <- agreeP

    rval0    
  }#fn0
  
  
  val0 <- kappam_fleiss0(idx = seq_len(nSj))
  kappa_j <- victorinox(est = kappam_fleiss0, idx = seq_len(nSj))
  
  # bias corrected estimate
  val <- as.numeric(val0 - kappa_j$bias_j)
  # standard error
  se_j <- kappa_j$se_j
  
  u <- val / se_j
  p.val <- 2 * stats::pnorm(q = abs(u), lower.tail = FALSE)
  # add results to return value
  rval <- c(rval, list(
    agreem = attr(val0, "agreeP"), value0 = as.numeric(val0),
    value = val, se_j = se_j, stat.name = "z", statistic = u, p.value = p.val)
  )
  
  
  # SE0 and detail
  if (variant == "fleiss" && (nCat == 2 || rval$balanced)) {
    
    # avg number of ratings given per subject
    rt_cnt <- rval$ratings / nSj # == mean(rt_cnt_sj)
    pj <- .colSums(tab_cnt_sj, m = nSj, n = nCat) / rval$ratings #prop. of categories overall
    qj <- 1 - pj
    pj.qj <- crossprod(pj, qj)[1L]
    
    
    if (nCat == 2) {
      
      # cf Fleiss (2003), 18.3, p. 613 eq (18.46) for SE0 (for kappa=0)
      rt_cnt_H <- 1 / mean(1 / rt_cnt_sj) #harmonic mean
      pq <- pj[1L] * qj[1L] #== prod(pj)
      SEkappa0 <- sqrt(2 * (rt_cnt_H-1) + ((rt_cnt - rt_cnt_H) * (1 - 4 * pq)) / (rt_cnt * pq)) /
        ((rt_cnt - 1) * sqrt(nSj * rt_cnt_H))
      
      # Lipsitz, Laird and Brennan ("Simple moment estimates ..", 1994) propose SE for any kappa (not only SE0),
      #+using a estimating equation approach. See also "Estimating the K-coefficient from a select sample" (2001)
      #+They claim their approach is easy to enhance to nCat > 2
      
    } else {
      
      stopifnot(rval$balanced)
      # alternative kappa formula in this setting, cf Fleiss (2003), p.614, eq (18.50)
      # 1 - (nSj * rt_cnt^2 - sum(tab_sj_cnt^2)) / (nSj * rt_cnt * (rt_cnt - 1) * pj.qj)
      
      # cf Fleiss (2003) 18.3, p. 616
      SEkappa0 <- sqrt(2) / (pj.qj * sqrt(nSj * rt_cnt * (rt_cnt - 1))) * 
        sqrt(pj.qj^2 - crossprod(pj * qj, qj - pj)[1L]) #SE0
      
      
      # kappa for each category level dichotomized (level vs rest)?
      if (detail) {
        # XXX think of calculating detail also for un-balanced case and with other variants than Fleiss?!
        #cf Fleiss (2003) 18.3, p.614, eq (18.50)
        stopifnot(rval$balanced)
        value_j <- 1 - diag(crossprod(tab_cnt_sj, rt_cnt - tab_cnt_sj)) / (nSj * rt_cnt * (rt_cnt - 1) * pj * qj)
        # NaN to NA
        is.na(value_j) <- !is.finite(value_j)
        SEkappa0_j <- sqrt(2 / (nSj * rt_cnt * (rt_cnt - 1))) #SE0
        u_j <- value_j / SEkappa0_j
        #SEkappa0_j <- rep_len(SEkappa0_j, length(value_j))
        #is.na(SEkappa0_j) <- !is.finite(value_j)
        p.value_j <- 2 * (1 - stats::pnorm(abs(u_j)))
        detail_j <- cbind(kappa_j = value_j,
                          se0_j = SEkappa0_j,
                          z_j = u_j,
                          p.value_j = p.value_j)
        rownames(detail_j) <- ratingScale
        
        rval <- c(rval, list(detail = detail_j))
      } #fi detail
    } #esle
    
    u0 <- rval$value0 / SEkappa0
    p.value0 <- 2 * stats::pnorm(abs(u), lower.tail = FALSE) #two-sided
    
    rval <- c(rval,
              se0 = SEkappa0, statistic0 = u0, p.value0 = p.value0)
  }#fi fleiss-variant (2cat or balanced)
  
  rval
}



#' Agreement of a group of nominal-scale raters with a gold standard
#'
#' First, Cohen's kappa is calculated between each rater against the gold
#' standard which is taken from the 1st column. The average of these kappas is
#' returned as 'kappam_gold0'. The variant setting (`robust=`) is forwarded to
#' Cohen's kappa. A bias-corrected version 'kappam_gold' and a corresponding
#' confidence interval are provided as well via the jackknife method.
#'
#' @examples
#' # matrix with subjects in rows and raters in columns.
#' # 1st column is taken as goldstandard
#' m <- matrix(c("O", "G", "O",
#'               "G", "G", "R",
#'               "R", "R", "R",
#'               "G", "G", "O"), ncol = 3, byrow = TRUE)
#' kappam_gold(m)
#'
#' @param ratings matrix subjects by raters
#' @param robust flag. Use robust estimate for random chance of agreement by
#'   Brennan-Prediger?
#' @param ratingScale Possible levels for the rating. Or `NULL`.
#' @param conf.level confidence level for confidence interval
#' @returns list. agreement measures (raw and bias-corrected) kappa with
#'   confidence interval. Entry `raters` refers to the number of tested raters,
#'   not counting the reference rater
#' @export
kappam_gold <- function(ratings, robust = FALSE, ratingScale = NULL, conf.level = .95) {
  ratings <- as.matrix(ratings)
  
  # handle ratingScale first, before dropping incomplete observations
  if (is.null(ratingScale)) {
    # sort also drops NAs (due to na.last=NA)
    ratingScale <- sort(unique(as.character(ratings)), na.last = NA)
  } else {
    if (anyDuplicated(ratingScale)) {
      stop("Duplicated entries in rating scale!", call. = FALSE)
    }
    ratingScale <- as.character(ratingScale)
    
    if (!all(as.character(ratings[!is.na(ratings)]) %in% ratingScale)) {
      stop("Ratings ", paste0(unique(as.character(ratings)), collpase="*"),
           " do not match provided rating scale ", paste0(ratingScale, collapse="*"), call. = FALSE)
    }
  }#esle
  nCat <- length(ratingScale)
  if (nCat <= 1L) {
    stop("Rating scale needs **at least two** levels!", call. = FALSE)
  }
  
  # gold standard
  subjGoldIdx <- which(!is.na(ratings[,1L]))
  if (!length(subjGoldIdx)) {
    stop("No subject with gold standard rating!", call. = FALSE)
  }
  
  # keep only subjects where gold-standard is given
  ratings <- ratings[subjGoldIdx,]
  # drop raters with only NA-ratings
  ratings <- ratings[, .colSums(!is.na(ratings),
                                m = NROW(ratings), n = NCOL(ratings)) >= 1L]
  
  nSj <- NROW(ratings)
  nRaters <- NCOL(ratings) # each rater is in a column
  stopifnot(nSj >= 1L, nRaters >= 2L)
  
  # raw kappa gold. 
  # @param idx row index of data to use (used for jackknifing)
  # @param what character. Which quantity from kappa2-list to work with? Default "value" is Cohen's kappa.
  kappam_gold0 <- function(idx, what = "value") {
    # Cohen's kappa for all pairwise ratings
    k2L <- purrr::map(.x = 2L:nRaters,
                      .f = ~ kappa2(ratings[idx, c(1L, .x), drop = FALSE],
                                    robust = robust, ratingScale = ratingScale))
    # drop invalid cases (no valid subjects)
    k2L <- purrr::compact(k2L)
    
    # build weighted average depending on available data
    #XXX better weighting would be via inverse of squared SE of kappa2!
    stats::weighted.mean(x = purrr::map_dbl(k2L, what),
                         w = purrr::map_dbl(k2L, "subjects"))
  }#fn
  
  # raw agreement proportion (has no bias in example runs)
  agreem <- kappam_gold0(idx = seq_len(nSj), what = "agreem")
  
  value0 <- kappam_gold0(idx = seq_len(nSj), what = "value")
  kgold_j <- victorinox(est = kappam_gold0, idx = seq_len(nSj))
  
  # bias corrected estimate
  value <- value0 - kgold_j$bias_j
  # standard error
  se_j <- kgold_j$se_j
  
  
  # 95% CI for bias-corrected estimate
  # stats::qt(1 - (1-conf.level)/2, df = max(1, nSubj-1)) * se
  ci <- value + c(-1, 1) * stats::qnorm(1 - (1-conf.level)/2) * se_j
  
  # return:
  list(
    method = "Averaged Cohen's Kappa with gold standard",
    subjects = nSj, raters = nRaters-1, categories = nCat,
    agreem = agreem,
    value0 = value0, value = value,
    se_j = se_j, conf.level = conf.level,
    ci.lo = ci[[1]], ci.hi = ci[[2]], ci.width = diff(ci)
  )
}

#' Significance test for homogeneity of kappa coefficients
#'
#' When groups of different subjects are rated on a nominal scale. Assuming
#' independence of subjects and their ratings between groups a chi-squared test
#' for equality of kappa between these groups is performed. The test requires
#' estimates of kappa and its standard error per group.
#'
#' A common overall kappa coefficient across groups is estimated. The test
#' statistic assesses the weighted squared deviance of the individual kappas
#' from the overall kappa estimate. The weights depend on the provided standard
#' errors.
#'
#' @examples
#' # script concordance test on 34 clinical situations,
#' # rated by 39 students and 11 experts
#' kappa_stud <- kappam_fleiss(SC_test[, 1:39])
#' kappa_expert <- kappam_fleiss(SC_test[, 40:50])
#'
#' # compare student and expert agreement
#' kappa_test(kappas = list(kappa_stud, kappa_expert))
#'
#'
#' @param kappas list of kappas from different groups. It uses the kappa
#'   estimate and its standard error.
#' @param val character. Name of field to extract kappa coefficient estimate.
#' @param se character. Name of field to extract standard error of kappa.
#' @param conf.level numeric. confidence level of confidence interval for
#'   overall kappa
#' @returns list containing the test results, including the entries `statistic`
#'   and `p.value` (class `htest`)
#' @references Joseph L. Fleiss, Statistical Methods for Rates and Proportions,
#'   3rd ed., 2003, section 18.1
#' @export
kappa_test <- function(kappas, val = "value0", se = "se0", conf.level = 0.95) {
  # check inputs
  isOK <- TRUE
  stopifnot(is.list(kappas), length(kappas) >= 2, nzchar(val), nzchar(se))
  purrr::iwalk(kappas, function(k, i) {
    if (!is.list(k) || !all(c(val, se) %in% names(k))) {
      isOK <- FALSE
      warning("Kappa provided as ",i, " is not a list or lacks requested entries!",
              call. = FALSE)
    }#fi
  })
  stopifnot(isOK)
  
  k <- purrr::map_dbl(kappas, val)
  w <- purrr::map_dbl(kappas, se)^-2  # inverse variance weighting
  kappa_ov <- c(`overall kappa` = stats::weighted.mean(x = k, w = w))
  
  chisq_df <- c(df = length(kappas)-1)
  
  teststat <- c(`X^2`=sum(w * (k - kappa_ov)^2))
  pval <- stats::pchisq(q = teststat, df = chisq_df, lower.tail = FALSE)
  ci_ov <- kappa_ov + c(-1, 1) * stats::qnorm(p=1-(1-conf.level)/2) * sum(w)^-.5
  attr(ci_ov, which = "conf.level") <- conf.level
  
  structure(list(method = "Test for Equality of Kappa in different groups",
                 data.name = paste(deparse1(substitute(kappas)), "(containing", length(kappas), "kappa values)"),
                 estimate = kappa_ov,
                 statistic = teststat, p.value = pval, 
                 alternative = "two.sided", parameter = chisq_df,
                 conf.int = ci_ov),
            class = "htest")
}


#' Simulate rating data and calculate agreement with gold standard
#' 
#' The function generates simulation data according to given categories and probabilities.
#' and can repeatedly apply function [kappam_gold()].
#' Currently, there is no variation in probabilities from rater to rater,
#' only sampling variability from multinomial distribution is at work.
#' 
#' 
#' This function is future-aware for the repeated evaluation of [kappam_gold()]
#' that is triggered by this function.
#' 
#' @examples
#' # repeatedly estimate agreement with goldstandard for simulated data
#' simulKappa(nRater = 8, cats = 3, nSubj = 11,
#'            # assumed prob for classification by raters
#'            probs = matrix(c(.6, .2, .1, # subjects of cat 1
#'                             .3, .4, .3, # subjects of cat 2
#'                             .1, .4, .5  # subjects of cat 3
#'            ), nrow = 3, byrow = TRUE))
#' 
#' 
#' @param nRater numeric. number of raters.
#' @param cats categories specified either as character vector or just the
#'   numbers of categories.
#' @param nSubj numeric. number of subjects per gold standard category. Either a
#'   single number or as vector of numbers per category, e.g. for non-balanced
#'   situation.
#' @param probs numeric square matrix (nCat x nCat) with classification
#'   probabilities. Row `i` has probabilities of rater categorization for
#'   subjects of category `i` (gold standard).
#' @param mcSim numeric. Number of Monte-Carlo simulations.
#' @param simOnly logical. Need only simulation data? Default is `FALSE`.
#' @returns dataframe of kappa-gold on the simulated datasets or (when
#'   `simOnly=TRUE`) list of length `mcSim` with each element a simulated data
#'   set with goldrating in first column and then the raters.
#' @export
simulKappa <- function(nRater, cats, nSubj, probs, mcSim = 10, simOnly=FALSE) {
  # check input
  if (is.numeric(cats)) {
    if (length(cats) != 1 || cats <= 2 || cats > 26^2) {
      stop("The number of categories given is invalid!", call. = FALSE)
    }#fi
    cats <- if (cats <= 26) {
      LETTERS[seq_len(cats)]
    } else {
      with(tidyr::expand_grid(c1 = LETTERS, c2 = LETTERS),
           paste0(c1, c2))[seq_len(cats)]
    }
  }#fi
  
  stopifnot(is.character(cats) || is.factor(cats))
  cats <- unique(as.character(cats))
  ncats <- length(cats)
  if (ncats <= 1) {
    stop("We need at least 2 categories!", call. = FALSE)
  }
  
  if (length(nSubj) == 1L) {
    nSubj <- rep.int(nSubj, times = ncats)
  }#fi
  
  stopifnot(is.numeric(nSubj), length(nSubj) == ncats,
            all(is.finite(nSubj)), all(nSubj >= 1L))
  
  if (missing(probs)) {
    stop("Please specify the assumed probabilities how subjects are rated.",
         call. = FALSE)
  }#fi
  stopifnot(is.matrix(probs), is.numeric(probs), NROW(probs) == ncats,
            NCOL(probs) == NROW(probs))
  nSubjTotal <- sum(nSubj)
  
  
  classif_of_rater <- vector(mode = 'list', length = ncats)
  classif_per_rater <- replicate(n = nRater, expr = {
    # currently, no variation in pp from rater to rater. could be added here.
    #+only sampling variability from multinomial distribution is at work
    for (ctgry in seq_along(cats)) {
      
      # raters act independently of each other,
      #+propensity of a subject to a category is not modelled individually
      # the order of ratings does not play a role here as we only look to the gold-standard
      # Within raters, the lexical ordering would increase the agreement
      classif_of_rater[[ctgry]] <- apply(
        X = stats::rmultinom(n=mcSim, size = nSubj[[ctgry]], prob = probs[ctgry,]),
        MARGIN = 2,
        FUN = function(x) rep.int(cats, times = x),
        simplify = FALSE)
    }#rof ctgry
    
    # per rater, combine the classifications for the subjects of the different categories
    purrr::pmap(classif_of_rater, .f = c)
  }, simplify = FALSE)
  rm(classif_of_rater)
  
  names(classif_per_rater) <- paste0("R", seq_len(nRater))
  
  
  goldStdOutcome <- rep.int(cats, times = nSubj)
  # build simulated data set per simulation round
  simData <- purrr::pmap(.l = classif_per_rater,
                         .f = function(...) # pass on arguments to tibble (in order to become columns)
                           tibble(gold = goldStdOutcome, # gold standard as first column
                                  ...,
                                  .rows = nSubjTotal))
  
  if (isTRUE(simOnly)) return(simData)
  
  resL <- future.apply::future_lapply(X = simData, FUN = kappam_gold)
  
  # return
  tibble(nRater = nRater, nSubjTotal = nSubjTotal,
         kappam_gold = purrr::map_dbl(.x = resL, "value"),
         ci_lo = purrr::map_dbl(.x = resL, "ci.lo"),
         ci_hi = purrr::map_dbl(.x = resL, "ci.hi"),
         ci_halfwidth = purrr::map_dbl(.x = resL, "ci.width") / 2L)
}

