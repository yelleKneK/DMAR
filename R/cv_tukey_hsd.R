# Critical value for the Tukey Honestly Significant Difference (HSD) test
#' Provides the Critical Value for the Tukey Honestly Significant Difference (HSD) Test
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate). For the Tukey HSD test,
#'   the full \code{alpha_level} applies to the upper tail of the Studentized range distribution
#'   (see Details).
#' @param df The error degrees of freedom from the ANOVA (a positive number; typically
#'   \code{N - groups}).
#' @param groups The number of groups whose means are being compared (an integer of at
#'   least 2).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value in a output style (a \code{data.frame} with one
#'   row per critical value, following the format used by \code{\link{cv_t}} and
#'   \code{\link{cv_z}}).
#'
#' @details The Tukey HSD test compares all pairs of group means using the Studentized
#'   range distribution (\code{\link[stats]{qtukey}}). Because the Studentized range is the
#'   absolute difference between the largest and smallest sample means (scaled by a standard
#'   error), it is non-negative and the associated distribution has support on
#'   \eqn{[0, \infty)}. As a consequence, the Type I error rate \code{alpha_level} is \emph{not}
#'   split between two tails in the way it is for the (symmetric) \emph{t}- and
#'   \emph{z}-distributions in \code{\link{cv_t}} and \code{\link{cv_z}}; rather, the full
#'   \code{alpha_level} applies to the upper tail.
#'
#'   The reported critical value is on the scale used for pairwise comparisons of group
#'   means, i.e., \eqn{q_{1-\alpha, k, df} / \sqrt{2}}, where \eqn{k} is the number of
#'   groups. A pair of means is declared significantly different when the absolute
#'   standardized difference between them exceeds this critical value.
#'
#'   The Tukey HSD critical value is a quantile of the Studentized range
#'   distribution, which base R supplies through \code{\link[stats]{qtukey}},
#'   so unlike \code{\link{cv_dunnett}} and \code{\link{cv_smm}} this function
#'   needs no multivariate distribution machinery and does not require the
#'   \pkg{mvtnorm} package. Maxwell, Delaney, and Kelley (2027, Chapter 5)
#'   develop the Tukey method as the procedure for all-pairwise comparisons
#'   within the multiple-comparisons problem.
#'
#' @references
#' Tukey, J. W. (1953). \emph{The problem of multiple comparisons}.
#'   Unpublished manuscript, Princeton University.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem,
#'   where the Tukey method for all-pairwise comparisons is developed.)
#'
#' @examples
#' # Following the all-pairwise comparisons setting of Maxwell, Delaney, and
#' # Kelley (2027, Chapter 5): every pair of k = 3 group means is compared,
#' # with 27 error degrees of freedom, holding the family-wise error rate at
#' # .05.
#' cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
#'
#' # Using DMAR's test_market data (6 marketing panels, N = 24).
#' fit <- aov(brand_movement ~ panel, data = test_market)
#' cv_tukey_hsd(
#'   alpha_level  = .05,
#'   df     = df.residual(fit),
#'   groups = nlevels(test_market$panel)
#' )
#'
#' # A more stringent alpha with a simple (non-verbose) result.
#' cv_tukey_hsd(alpha_level = .01, df = 27, groups = 3, verbose = FALSE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_t}}, \code{\link{cv_z}}, \code{\link[stats]{TukeyHSD}},
#'   \code{\link[stats]{qtukey}}
#'
#' @keywords htest design
#'
#' @family critical values
#'
#' @export
#' @import stats

cv_tukey_hsd <- function(alpha_level, df, groups, verbose = TRUE) {
  # Degrees of freedom checks (mirrors cv_t).
  if (missing(df)) stop("You must specify the degrees of freedom (i.e., 'df'), which must be a positive number.")
  if (!df > 0) stop("You must specify a positive value for the degrees of freedom.")

  # Groups check (specific to Tukey HSD).
  if (missing(groups)) stop("You must specify the number of groups (i.e., 'groups').")
  if (groups < 2) stop("This function is appropriate for situations in which there are two or more groups.")

  # Alpha check. For the Studentized range distribution, alpha_level is not split between tails
  # (see Details), so only a single alpha_level argument is needed.
  if (missing(alpha_level)) stop("You must specify 'alpha_level'.")
  if (alpha_level <= 0 || alpha_level >= 1) stop("Specify 'alpha_level' to be greater than zero and less than 1.")

  # Critical value on the Studentized range scale.
  q_critical <- qtukey(p = 1 - alpha_level, nmeans = groups, df = df)

  # Rescale to the scale used for pairwise mean comparisons.
  value <- q_critical / sqrt(2)

  # Areas on the Studentized range scale (probabilities are invariant under the monotonic
  # rescaling by 1/sqrt(2), so these also correspond to the reported value).
  area_less    <- ptukey(q_critical, nmeans = groups, df = df, lower.tail = TRUE)
  area_greater <- ptukey(q_critical, nmeans = groups, df = df, lower.tail = FALSE)

  term <- "upper_cv"

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term = term, value = value, area_less = area_less, area_greater = area_greater)))
  }
  if (verbose == FALSE) {
    return(.as_dmar_tbl(data.frame(term = term, value = value)))
  }
}
