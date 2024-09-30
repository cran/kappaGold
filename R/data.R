#' Script concordance test (SCT).
#'
#' In medical education, the script concordance test (SCT) (Charlin, Gagnon,
#' Sibert, & Van der Vleuten, 2002) is used to score physicians or medical
#' students in their ability to solve clinical situations as compared to answers
#' given by experts. The test consists of a number of items to be evaluated on a
#' 5-point Likert scale.
#'
#' Each item represents a clinical situation (called an 'assumption') likely to
#' be encountered in the physician’s practice. The situation has to be unclear,
#' even for an expert. The task of the subjects being evaluated is to consider
#' the effect of new information on the assumption to solve the situation. The
#' data incorporates 50 raters, 39 students and 11 experts.
#' 
#' Each rater judges the same 34 assumptions.
#'
#' @format A matrix with 34 rows and 50 columns. Columns 1 to 39 are student
#'   raters, columns 40 to 50 are experts. Each rater applies to each clinical
#'   situation one of five levels ranging from -2 to 2 with the following
#'   meaning:
#' \describe{
#'   \item{-2}{The assumption is practically eliminated;}
#'   \item{-1}{The assumption becomes less likely;}
#'   \item{0}{The information has no effect on the assumption;}
#'   \item{+1}{The assumption becomes more likely;}
#'   \item{+2}{The assumption is virtually the only possible one.}
#' }
#' @source Sophie Vanbelle (personal communication, 2021)
#' @references Vanbelle, S., Albert, A. Agreement between Two Independent Groups
#'   of Raters. Psychometrika 74, 477–491 (2009).
#'   \doi{10.1007/s11336-009-9116-1}
"SC_test"


#' Staging of colorectal carcinoma
#'
#' Staging of carcinoma is done by different medical professions. Gold standard
#' is the (histo-)pathological rating of a tissue sample but this information
#' typically only becomes available late, after surgery. However prior to
#' surgery the carcinoma is also staged by radiologists in the clinical setting
#' on the basis of MRI scans.
#'
#' These fictitious data were inspired by the OCUM trial. The simulation uses
#' the following two assumptions: over-staging occurs more frequently than
#' under-staging and an error by two categories is less likely than an error by
#' only one category.
#' 
#' Stages conform to the UICC classification according to the TNM
#' classification. Note that cases in stage IV do not appear in this data set
#' and that the following description of stages is simplified.
#' 
#' 1. **I** Until T2, N0, M0
#' 2. **II** From T3, N0, M0
#' 3. **III** Any T, N1/N2, M0
#' 
#' @format
#' A data frame with 21 observations and 6 variables:
#' \describe{
#'   \item{patho}{the (histo-)pathological staging (gold standard) with categories I, II or III}
#'   \item{rad1}{the clinical staging with categories I, II or III by radiologist 1}
#'   \item{rad2}{the clinical staging with categories I, II or III by radiologist 2}
#'   \item{rad3}{the clinical staging with categories I, II or III by radiologist 3}
#'   \item{rad4}{the clinical staging with categories I, II or III by radiologist 4}
#'   \item{rad5}{the clinical staging with categories I, II or III by radiologist 5}
#' }
#' 
#' 
#' @source simulated data
#' @references Kreis, M. E. et al., MRI-Based Use of Neoadjuvant
#'   Chemoradiotherapy in Rectal Carcinoma: Surgical Quality and
#'   Histopathological Outcome of the OCUM Trial
#'   \doi{10.1245/s10434-019-07696-y}
"stagingData"


#' Psychiatric diagnoses
#'
#' N = 30 patients were given one of k = 5 diagnoses by some n = 6
#' psychiatrists out of 43 psychiatrists in total. The diagnoses are
#' 1. Depression
#' 2. PD (=Personality Disorder)
#' 3. Schizophrenia
#' 4. Neurosis
#' 5. Other
#' 
#' A total of 43 psychiatrists provided diagnoses. In the actual
#' study (Sandifer, Hordern, Timbury, & Green, 1968), between 6 and 10
#' psychiatrists from the pool of 43 were unsystematically selected to diagnose
#' a subject. Fleiss randomly selected six diagnoses per subject to bring the
#' number of assignments per patient down to a constant of six.
#'
#'
#'
#' As there is not a fixed set of six raters the ratings from the same column
#' are not related to each other. Therefore, compared to the dataset with the
#' same name in package `irr`, we applied a permutation of the six ratings.
#'
#' @format ## `diagnoses`
#' A matrix with 30 rows and 6 columns:
#' \describe{
#'   \item{rater1}{1st rating of some six raters}
#'   \item{rater2}{2nd rating of some six raters}
#'   \item{rater3}{3rd rating of some six raters}
#'   \item{rater4}{4th rating of some six raters}
#'   \item{rater5}{5th rating of some six raters}
#'   \item{rater6}{6th rating of some six raters}
#' }
#' 
#' @references Sandifer, M. G., Hordern, A., Timbury, G. C., & Green, L. M.
#'   Psychiatric diagnosis: A comparative study in North Carolina, London and
#'   Glasgow. British Journal of Psychiatry, 1968, 114, 1-9.
#' @references Fleiss, J. L. Measuring nominal scale agreement among many
#'   raters. Psychological Bulletin, 1971, 76(5), 378–382.
#'   \doi{10.1037/h0031619}
#' @seealso This dataset is also available as `diagnoses` in the irr-package on
#'   CRAN.
"diagnoses"
