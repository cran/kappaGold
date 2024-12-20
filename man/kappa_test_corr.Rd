% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kappa_inference.R
\name{kappa_test_corr}
\alias{kappa_test_corr}
\title{Test for homogeneity of kappa in correlated groups}
\usage{
kappa_test_corr(
  ratings,
  grpIdx,
  kappaF,
  kappaF_args = list(),
  B = 100,
  alternative = "two.sided",
  conf.level = 0.95
)
}
\arguments{
\item{ratings}{matrix. ratings as sbj x raters, including the multiple groups
to be tested}

\item{grpIdx}{list. Comprises numeric index vectors per group. Each group is
defined as set of raters (i.e., columns)}

\item{kappaF}{function or list of functions. kappa function to apply on each group.}

\item{kappaF_args}{list. Further arguments for the kappa function. By
default, these settings apply to all groups, but the settings can be
specified per group (as list of lists).}

\item{B}{numeric. number of bootstrap samples. At least 1000 are recommended
for stable results.}

\item{alternative}{character. Direction of alternative. Currently only
\code{'two.sided'} is supported.}

\item{conf.level}{numeric. confidence level for confidence intervals}
}
\value{
list. test results as class \code{htest}. The confidence interval shown
by \code{print} refers to the 1st difference \code{k1-k2}.
}
\description{
Bootstrap test on kappa based on data with common subjects. The differences
in kappa between all groups (but first) relative to first group (e.g., Group
2 - Group 1) are considered.
}
\section{Note}{
Due to limitations of the \code{htest} print method the confidence interval shown
by \code{print} refers to the 1st difference \code{k1-k2}. If there are more than 2
groups access all confidence intervals via entry \code{conf.int}.
}

\examples{
# Compare Fleiss kappa between students and expert raters
# For real analyses use more bootstrap samples (B >= 1000)
kappa_test_corr(ratings = SC_test, grpIdx = list(S=1:39, E=40:50), B = 125,
                kappaF = kappam_fleiss,
                kappaF_args = list(variant = "fleiss", ratingScale=-2:2))

}
