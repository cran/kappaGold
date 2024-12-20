% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kappa.R
\name{kappam_fleiss}
\alias{kappam_fleiss}
\title{Fleiss' kappa for multiple nominal-scale raters}
\usage{
kappam_fleiss(
  ratings,
  variant = c("fleiss", "conger", "robust", "uniform"),
  detail = FALSE,
  ratingScale = NULL
)
}
\arguments{
\item{ratings}{matrix (subjects by raters), containing the ratings}

\item{variant}{Which variant of kappa? Default is Fleiss (1971). Other options are Conger (1980) or robust variant.}

\item{detail}{Should category-wise Kappas be computed? Only available for the Fleiss (1971) variant.}

\item{ratingScale}{Specify possible levels for the rating. Default \code{NULL} means to use all unique levels from the sample.}
}
\value{
list containing Fleiss's kappa agreement measure (value) or \code{NULL} if no subjects
}
\description{
When multiple raters judge subjects on a nominal scale we can assess their agreement with Fleiss' kappa.
It is a generalization of Cohen's Kappa for two raters and there are different variants how to assess chance agreement.
}
\details{
Different \strong{variants} of Fleiss' kappa are implemented.
By default (\code{variant="fleiss"}), the original Fleiss Kappa (1971) is calculated, together with an asymptotic standard error and test for kappa=0.
It assumes that the raters involved are not assumed to be the same (one-way ANOVA setting).
The marginal category proportions determine the chance agreement.
Setting \code{variant="conger"} gives the variant of Conger (1980) that reduces to Cohen's kappa when m=2 raters.
It assumes identical raters for the different subjects (two-way ANOVA setting).
The chance agreement is based on the category proportions of each rater separately.
Typically, the Conger variant yields slightly higher values than Fleiss kappa.
\code{variant="robust"} assumes a chance agreement of two raters to be simply 1/q, where q is the number of categories (uniform model).
}
\examples{
# 4 subjects were rated by 3 raters in categories "1", "2" or "3"
# organize ratings as matrix with subjects in rows and raters in columns
m <- matrix(c("3", "2", "3",
              "2", "2", "1",
              "1", "3", "1",
              "2", "2", "3"), ncol = 3, byrow = TRUE)
kappam_fleiss(m)

# show category-wise kappas -----
kappam_fleiss(m, detail = TRUE)

}
\seealso{
\code{\link[irr:kappam.fleiss]{irr::kappam.fleiss()}}
}
