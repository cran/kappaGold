# inference on kappa coefficients



#' Significance test for homogeneity of kappa coefficients in independent groups
#'
#' The null hypothesis states that the kappas for all involved groups are the
#' same ("homogeneous"). A prerequisite is that the groups are independent of
#' each other, this means the groups are comprised of different subjects and
#' each group has different raters. Each rater employs a nominal scale. The test
#' requires estimates of kappa and its standard error per group.
#'
#' A common overall kappa coefficient across groups is estimated. The test
#' statistic assesses the weighted squared deviance of the individual kappas
#' from the overall kappa estimate. The weights depend on the provided standard
#' errors. Under H0, the test statistics is chi-square distributed.
#'
#' @examples
#' # three independent agreement studies (different raters, different subjects)
#' # each study involves two raters that employ a binary rating scale
#' k2_studies <- lapply(agreem_binary, kappa2)
#'
#' # combined estimate and test for homogeneity of kappa
#' kappa_test(kappas = k2_studies, val = "value", se = "se")
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
  stopifnot(is.numeric(conf.level), length(conf.level) == 1,
            conf.level > 0, conf.level < 1)
  purrr::iwalk(kappas, function(k, i) {
    if (!is.list(k) || !all(c(val, se) %in% names(k))) {
      isOK <- FALSE
      warning("Kappa provided in position ", i,
              " is not a list or lacks requested entries!",
              call. = FALSE)
    }#fi
  })
  
  if (!isOK) {
    return(invisible(NULL))
  }
  
  k <- purrr::map_dbl(kappas, val)
  v <- purrr::map_dbl(kappas, se)^2
  # inverse variance weighting
  kappa_ov <- c(`overall kappa` = stats::weighted.mean(x = k, w = 1/v))
  
  chisq_df <- c(df = length(kappas)-1)
  
  teststat <- c(`X^2`=sum((k - kappa_ov)^2 / v))
  pval <- stats::pchisq(q = teststat, df = chisq_df, lower.tail = FALSE)
  ci_ov <- kappa_ov + c(-1, 1) * stats::qnorm(p=1-(1-conf.level)/2) * sum(1/v)^-.5
  ci_ov[2] <- min(ci_ov[2], 1)
  attr(ci_ov, which = "conf.level") <- conf.level
  
  structure(list(method = "Test for Equality of Kappa in independent groups",
                 data.name = paste(deparse1(substitute(kappas)), "(containing",
                                   length(kappas), "kappa values)"),
                 estimate = kappa_ov,
                 statistic = teststat, p.value = pval, 
                 alternative = "two.sided", parameter = chisq_df,
                 conf.int = ci_ov),
            class = "htest")
}


#' Test for homogeneity of kappa in correlated groups
#'
#' Bootstrap test on kappa based on data with common subjects. The differences
#' in kappa between all groups (but first) relative to first group (e.g., Group
#' 2 - Group 1) are considered.
#'
#' # Note
#' Due to limitations of the `htest` print method the confidence interval shown
#' by `print` refers to the 1st difference `k1-k2`. If there are more than 2
#' groups access all confidence intervals via entry `conf.int`.
#'
#' @param ratings matrix. ratings as sbj x raters, including the multiple groups
#'   to be tested
#' @param grpIdx list. Comprises numeric index vectors per group. Each group is
#'   defined as set of raters (i.e., columns)
#' @param kappaF function or list of functions. kappa function to apply on each group.
#' @param kappaF_args list. Further arguments for the kappa function. By
#'   default, these settings apply to all groups, but the settings can be
#'   specified per group (as list of lists).
#' @param B numeric. number of bootstrap samples. At least 1000 are recommended
#'   for stable results.
#' @param alternative character. Direction of alternative. Currently only
#'   `'two.sided'` is supported.
#' @param conf.level numeric. confidence level for confidence intervals
#' @returns list. test results as class `htest`. The confidence interval shown
#'   by `print` refers to the 1st difference `k1-k2`.
#' @export
#'
#' @examples
#' # Compare Fleiss kappa between students and expert raters
#' # For real analyses use more bootstrap samples (B >= 1000)
#' kappa_test_corr(ratings = SC_test, grpIdx = list(S=1:39, E=40:50), B = 125,
#'                 kappaF = kappam_fleiss,
#'                 kappaF_args = list(variant = "fleiss", ratingScale=-2:2))
#' 
kappa_test_corr <- function(ratings, grpIdx,
                            kappaF, kappaF_args = list(),
                            B = 100, alternative = "two.sided", conf.level = .95) {
  
  stopifnot(is.matrix(ratings), NCOL(ratings) > 2, NROW(ratings) > 1)
  stopifnot(is.list(grpIdx), length(grpIdx) >= 2, all(lengths(grpIdx) >= 2),
            all(vapply(grpIdx, FUN = is.numeric, FUN.VALUE = TRUE)))
  nGrps <- length(grpIdx)
  
  stopifnot(is.numeric(conf.level), length(conf.level) == 1L,
            conf.level > 0, conf.level < 1)
  stopifnot(is.character(alternative), alternative == "two.sided")
  
  stopifnot(is.numeric(B), length(B) == 1L, B > 13)
  B <- as.integer(B)
  
  stopifnot(is.function(kappaF) || is.list(kappaF))
  
  if (is.function(kappaF)) {
    if (!rlang::has_name(formals(kappaF), name = "ratingScale")) {
      warning("Provided function in kappaF= does not allow to set `ratingScale=`",
              call. = FALSE)
    }#fi
  } else {
    stopifnot(is.list(kappaF))
    stopifnot(all(purrr::map_lgl(kappaF, .f = is.function)))
    
    if (length(kappaF) != nGrps) {
      stop("kappaF= Please provide a list of kappa-functions, one per group!",
           call. = FALSE)
    }
    if (!all(purrr::map_lgl(.x = kappaF,
                            .f = \(kF) rlang::has_name(formals(kF), name = "ratingScale")))) {
      warning("Not all functions in kappaF= do allow to set `ratingScale=`",
              call. = FALSE)
    }
  }
  
  stopifnot(is.list(kappaF_args))
  # check if we have separate kappaF-arguments per group?
  kappaF_args_depth <- purrr::pluck_depth(kappaF_args)
  stopifnot(kappaF_args_depth > 1)
  if (kappaF_args_depth == 2L) {
    stopifnot(rlang::is_named2(kappaF_args))
    if (!"ratingScale" %in% names(kappaF_args)) {
      warning("Providing ratingScale= is recommended!", call. = FALSE)
    }#fi
    # rep keeps the names (rep_len does not!)
    kappaF_args <- rep(list(kappaF_args), times = nGrps)
  } else if (kappaF_args_depth == 3L) {
    if (length(kappaF_args) != nGrps) {
      stop("Please provide kappa-arguments for each group as list!", call. = FALSE)
    }
    if (!all(purrr::map_lgl(.x = kappaF_args, .f = \(kFa) rlang::has_name(kFa, name = "ratingScale")))) {
      warning("Providing ratingScale= is recommended for every group!", call. = FALSE)
    }
  }
  
  # contract: have arguments for kappaF for each group!
  stopifnot(length(kappaF_args) == nGrps)
  
  # get kappas per group on bootstrapped sample
  # argument i (index) is ignored: we use sample() to draw bootstrap samples
  kStarFun <- function (idx) {
    # draw bootstrap sample of subjects (=rows)
    sjIdx <- sample(NROW(ratings), size = NROW(ratings), replace = TRUE)
    vapply(seq_len(nGrps), FUN = function (igr) {
      kF <- if (is.list(kappaF)) kappaF[[igr]] else kappaF
      do.call(kF, args = c(list(ratings = ratings[sjIdx, grpIdx[[igr]], drop = FALSE]),
                               kappaF_args[[igr]]))$value
    }, FUN.VALUE = 9.9)
    
  }#fn
  
  kppsStar <- future.apply::future_vapply(X = seq_len(B),
                                          FUN = kStarFun,
                                          FUN.VALUE = rep.int(0, nGrps),
                                          future.seed = TRUE,
                                          future.packages = c("kappaGold"))
  
  
  M <- .rowMeans(kppsStar, m = nGrps, n = B)
  V <- stats::var(t(kppsStar)) #variance of bootstrapped kappas
  # contrast matrix
  C <- cbind(1, diag(-1, nrow = nGrps-1))
  Ck <- C %*% M
  # Hotelling's T-square
  T2 <- as.vector(crossprod(Ck, solve(tcrossprod(C %*% V, C)) %*% Ck))
  
  # confidence intervals for k1 - kj
  ciMat <- as.vector(Ck) + tcrossprod(sqrt(stats::qchisq(p = conf.level, df = nGrps-1) * diag(C %*% V %*% t(C))), c(-1,1))
  dimnames(ciMat) <- list(paste0("k1-k", 2:nGrps), c("L", "U"))
  attr(ciMat, which = "conf.level") <- conf.level
  
  # return value
  rval <- lst(estimate = as.vector(Ck),
              null.value = if (nGrps > 2) {
                c(`k1 - kj (for all j=2...)` = 0)
              } else {
                c(`k1 - k2` = 0)
              },
              statistic = c(T2=T2),
              parameter = list(k1 = M[[1]], se_k1 = sqrt(V[[1]]), B = B),
              data.name = paste(NROW(ratings), "subjects rated by", nGrps, "correlated groups of raters"),
              alternative = alternative,
              method = "2-sided Hotelling's T-square bootstrap test on kappa in correlated groups",
              conf.int = t(ciMat),
              p.value = stats::pchisq(q = T2, df = nGrps-1, lower.tail = FALSE)
  )
  class(rval) <- "htest"
  
  rval
}
