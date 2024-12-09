#' Agreement between two groups of raters
#'
#' This function expands upon Cohen's and Fleiss' Kappa as measures for
#' interrater agreement while taking into account the heterogeneity within each
#' group.
#'
#' Data need to be stored with raters in columns.
#' 
#' 
#' @param ratings matrix of subjects x raters for both groups of raters
#' @param refIdx numeric. indices of raters that constitute the reference group.
#'   Can also be all negative to define rater group by exclusion.
#' @param ratingScale character vector of the levels for the rating. Or `NULL`.
#' @param weights optional weighting schemes: `"unweighted"`,
#'   `"linear"`,`"quadratic"`
#' @param conf.level confidence level for interval estimation
#' @returns list. kappa agreement between two groups of raters
#' @references Vanbelle, S., Albert, A. Agreement between Two Independent Groups
#'   of Raters. Psychometrika 74, 477–491 (2009).
#'   \doi{10.1007/s11336-009-9116-1}
#' @examples
#' # compare student ratings with ratings of 11 experts
#' kappam_vanbelle(SC_test, refIdx = 40:50)
#' 
#' @export
kappam_vanbelle <- function(ratings, refIdx, ratingScale = NULL,
                          weights = c("unweighted", "linear", "quadratic"),
                          conf.level=.95) {
  
  stopifnot(is.numeric(refIdx), length(refIdx) >= 1,
            min(abs(refIdx)) >= 1, max(abs(refIdx)) <= NCOL(ratings))
  weights <- match.arg(weights)
  stopifnot(is.numeric(conf.level), length(conf.level) == 1L,
            conf.level > 0, conf.level < 1)
  
  ratings <- stats::na.omit(as.matrix(ratings))
  stopifnot(!is.null(ratings), is.matrix(ratings), NROW(ratings) >= 1)
  
  ratingsGr1 <- ratings[, refIdx, drop = FALSE]
  ratingsGr2 <- ratings[, -refIdx, drop = FALSE]
  nRat1 <- NCOL(ratingsGr1)
  nRat2 <- NCOL(ratingsGr2)
  
  nSubj <- NROW(ratings)
  
  if (is.null(ratingScale)) {
    ratingScale <- sort(unique(as.vector(ratings)))
  } else {
    if (anyDuplicated(ratingScale)) {
      stop("Duplicated entries in rating scale!", call. = FALSE)
    }#fi
    ratingScale <- as.character(ratingScale)
    
    if (!all(as.character(ratings) %in% ratingScale)) {
      stop("Ratings and provided rating scale do not match!", call. = FALSE)
    }#fi
  }
  
  # number of categories
  nCat <- length(ratingScale)
  
  # count ratings (columns) per subject (rows)
  nij1 <- t(apply(ratingsGr1, MARGIN = 1L,
                  FUN = function(subjR) tabulate(factor(subjR, levels = ratingScale), nbins = nCat)))
  nij2 <- t(apply(ratingsGr2, MARGIN = 1L,
                  FUN = function(subjR) tabulate(factor(subjR, levels = ratingScale), nbins = nCat)))
  colnames(nij1) <- colnames(nij2) <- ratingScale
  
  
  # weights in matrix format: "agreement" form with 1s on the diagonal
  w <- diag(nrow = nCat)
  w <- switch(weights,
              linear = 1L - abs(row(w) - col(w)) / (nCat-1L),
              quadratic = 1L - ((row(w) - col(w)) / (nCat-1L))^2L,
              unweighted = w,
              stop("Unknown weights given: ", sQuote(weights), call. = FALSE)
  )
  
  
  # inner helper function (also usable with jackknife)
  kappam_vanbelle0 <- function(idx = seq_len(nSubj)) {
    
    nSubj0 <- length(idx)
    nij1 <- nij1[idx,, drop = FALSE]
    nij2 <- nij2[idx,, drop = FALSE]
    
    pij1 <- nij1 / nRat1
    pij2 <- nij2 / nRat2
    kxk <- crossprod(pij1, pij2) / nSubj0
    
    stopifnot(NCOL(kxk) == nCat)
    
    # observed agreement  
    agreeP <- sum(w * kxk)
    
    #Hint:
    #+You can use apply for pij1 and pij2 and then take the mean over the result of pmax
    # entry1 <- max(sum(outer(pij1[1,], pij1[1,])*w), sum(outer(pij2[1,], pij2[1,])*w))
    maxAgreeP0 <- rep_len(NA_real_, length.out = nSubj0)
    
    for (i in seq_len(nSubj0)) {
      maxAgreeP0[i] <- max(sum(outer(pij1[i,], pij1[i,])*w), sum(outer(pij2[i,], pij2[i,])*w))
    }
    maxAgreeP <- mean(maxAgreeP0)
    
    # outer product of sums (margins) of group 1 & 2: relative frequencies of each category
    eij <- outer(rowSums(kxk), colSums(kxk))
    chanceP <- sum(eij * w) 
    
    # return
    (agreeP - chanceP) / (maxAgreeP - chanceP)
  }
  
  
  # jackknife -----
  value0 <- kappam_vanbelle0(idx = seq_len(nSubj))
  kgold_jk <- victorinox(est = kappam_vanbelle0, idx = seq_len(nSubj))
  
  # bias corrected estimate
  value <- value0 - kgold_jk$bias_j
  se <- kgold_jk$se_j
  
  # 95% CI for bias-corrected estimate
  ci <- value + c(-1L, 1L) * stats::qt(1L - (1L-conf.level)/2L, df = nSubj-1L) * se
  
  # return ----
  list(method = "Vanbelle's Kappa for two groups of raters (with bias-correction, 95% CI and optional weighting)",
       subjects = nSubj, raters = NCOL(ratings), categories = nCat,
       weights = weights,
       value0 = value0, value = value,
       se = se,
       conf.level = conf.level, ci.lo = ci[[1L]], 
       ci.hi = ci[[2L]], ci.width = ci[[2L]] - ci[[1L]])
}

