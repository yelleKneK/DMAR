# Critical value for the Tukey Honestly Significant Difference (HSD) test
#' Provides the Critical Value for the Tukey Honestly Significant Difference (HSD) Test
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate). For the Tukey HSD test,
#'   the full \code{alpha_level} applies to the upper tail of the Studentized range distribution
#'   (see Details).
#' @param df The error degrees of freedom from the ANOVA (a positive number; typically
#'   \code{N - n_groups}).
#' @param n_groups The number of groups whose means are being compared (an integer of at
#'   least 2).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value in a output style (a \code{data.frame} with one
#'   row per critical value, following the format used by \code{\link{cv_t}} and
#'   \code{\link{cv_z}}).
#'
#' @details The Tukey HSD test compares all pairs of group means using the Studentized
#'   range distribution. Because the Studentized range is the
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
#'   distribution, so unlike \code{\link{cv_dunnett}} and \code{\link{cv_smm}}
#'   this function needs no multivariate distribution machinery and does not
#'   require the \pkg{mvtnorm} package. Maxwell, Delaney, and Kelley (2027,
#'   Chapter 5) develop the Tukey method as the procedure for all-pairwise
#'   comparisons within the multiple-comparisons problem.
#'
#'   \strong{How the quantile is computed.} Over most of its range the critical
#'   value is base R's \code{\link[stats]{qtukey}}, which R documents as accurate
#'   to the fourth decimal place (it inverts \code{ptukey} by a secant search),
#'   and the areas are \code{\link[stats]{ptukey}}, which evaluates its integral
#'   with a 16-point Legendre formula and for which R states no accuracy.
#'   Compared with an independent computation, the two fall short in five
#'   places. (1) At small error df, increasingly so as the number of groups and
#'   the confidence level grow: with 20 groups and 2 error df,
#'   \code{qtukey(.99, 20, 2)} is 44.526 where the exact quantile is 37.943, so
#'   the critical value would be 31.485 rather than 26.830; below 2 df both
#'   functions return \code{NaN}. (2) At large error df: \code{qtukey(.999, 3,
#'   18000)} is \eqn{1.2 \times 10^{-5}} too large, and beyond 25000 df both
#'   return the infinite-df distribution. (3) Far in either tail:
#'   \code{qtukey(.9999, 10, 11)} is 12.008462 where the exact quantile is
#'   12.008580, so \code{cv_tukey_hsd(1e-4, df = 11, n_groups = 10)} would be
#'   8.491265 rather than 8.491348; at .99999 \code{qtukey} can be off by
#'   several units (with 2000 groups and 1000 error df it is 12.930 where the
#'   exact quantile is 9.872); in the lower tail it was more than
#'   \eqn{10^{-5}} off at .0001 (with 75 groups) and at .00001 can be far off
#'   (\code{qtukey(.00001, 3, 150)} is 0.000092 where the exact quantile is
#'   0.006023) or \code{NaN}. (4) With more than 50 groups \code{ptukey} is
#'   off by more than \eqn{10^{-6}} below its .75 point (by up to
#'   \eqn{8.3 \times 10^{-6}} with 2000 groups, and with 1500 or more groups
#'   it can return 0 where the probability is .001), and with thousands of
#'   groups at very large error df it is off just above that point as well
#'   (\code{ptukey(7.6250574, 5000, Inf)} is .7500010 where the exact
#'   probability is .7500000). (5) Anywhere, the secant
#'   search can stop short of the
#'   root of \code{ptukey}: \code{qtukey(.975, 10, 12)} is 6.003766 where the
#'   exact quantile is 6.003807, and \code{qtukey(.9975, 5, 23)} is 6.018665
#'   where it is 6.018546.
#'
#'   In these places DMAR computes the distribution function exactly, as the
#'   probability integral of the range,
#'   \eqn{P(W \le w) = k \int \phi(z) [\Phi(z) - \Phi(z - w)]^{k - 1} dz} for
#'   the range \eqn{W} of \eqn{k} standard normal variates
#'   (\eqn{k} the number of groups), integrated against the chi square
#'   distribution of the error variance by adaptive quadrature, and takes the
#'   critical value as the root of that function. It does so (a) below an error
#'   df that grows with the number of groups: fewer than 9 error df for up to 6
#'   groups, fewer than 11 for up to 10, 15 for up to 20, 21 for up to 50, 31
#'   for up to 100, 41 for up to 250, 47 for up to 400, 57 for up to 700, 71 for
#'   up to 1000, 81 for up to 1500, 91 for up to 2000, 101 for up to 3000, and
#'   111 for more; (b) above 10000 and up to \eqn{10^8} error df; (c) for every
#'   quantile at a probability above .999 (an \code{alpha_level} below .001)
#'   or below .001, at any error df; (d) with more than 50 groups,
#'   for every probability below .75 and every quantile below the .75 point
#'   (below the .80 point with more than 1000 groups), at any error df; and
#'   (e) wherever \code{qtukey} fails, or its value is more than
#'   \eqn{2 \times 10^{-6}} from the root of \code{ptukey} (the residual of
#'   \code{ptukey} there divided by its slope). These rules were set by
#'   comparison with two independent computations of the studentized range, as
#'   the cells in which \code{qtukey} was more than \eqn{10^{-5}} from the exact
#'   quantile or \code{ptukey} more than \eqn{10^{-6}} from the exact
#'   probability, with a margin of at least one measured step. R's values were
#'   compared with 2 to 10000 groups at error df up to 20000, with 2 to 20000
#'   groups at infinite df (which R uses for every df above 25000), and with 3
#'   to 10000 groups at \eqn{10^8 + 1}, \eqn{10^9}, and \eqn{10^{10}} df, at
#'   probabilities from .00001 to .99999 (52,892 cells for the second
#'   computation alone). Elsewhere the
#'   values are those of \code{qtukey} and \code{ptukey}; at the measured cells
#'   there they were within \eqn{9.6 \times 10^{-6}} of the exact quantile and
#'   \eqn{9.1 \times 10^{-7}} of the exact probability. The largest quantile
#'   error (700 groups, 57 df, at .9965) came from a dense probe
#'   between the grid's cells near the df bounds in (a), where the errors come
#'   closest to the \eqn{10^{-5}} rule; at the grid's own cells the largest
#'   was \eqn{8.3 \times 10^{-6}}. Between the measured cells larger errors
#'   cannot be ruled out, and beyond 10000 groups (20000 at infinite df) the
#'   df bound in (a) and the cuts in (d) are extrapolated, not measured. Above
#'   \eqn{10^8} error df the distribution is within
#'   \eqn{2.3 \times 10^{-7}} in probability (\eqn{1.4 \times 10^{-6}} in the
#'   quantile) of the infinite-df one, and the values are those of
#'   \code{qtukey} and \code{ptukey} at infinite df, except under (c), (d), and
#'   (e), where the exact infinite-df distribution is used; at the three df
#'   measured there DMAR's values were within
#'   \eqn{1.6 \times 10^{-6}} of the exact quantile and
#'   \eqn{9.1 \times 10^{-7}} of the exact probability. The exact computation
#'   agrees with the independent computations to \eqn{2 \times 10^{-8}} in the
#'   quantile at probabilities from .01 to .99999, and less closely far in the
#'   lower tail, where the quantile is small (to \eqn{4.6 \times 10^{-8}} at
#'   .001, \eqn{1.1 \times 10^{-7}} at .0001, and \eqn{4.8 \times 10^{-7}} at
#'   .00001); with two groups it reproduces the closed form
#'   \eqn{q = \sqrt{2}\, t_{1 - \alpha/2, \mathit{df}}}.
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
#' cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3)
#'
#' # Using DMAR's test_market data (6 marketing panels, N = 24).
#' fit <- aov(brand_movement ~ panel, data = test_market)
#' cv_tukey_hsd(
#'   alpha_level  = .05,
#'   df     = df.residual(fit),
#'   n_groups = nlevels(test_market$panel)
#' )
#'
#' # A more stringent alpha with a simple (non-verbose) result.
#' cv_tukey_hsd(alpha_level = .01, df = 27, n_groups = 3, verbose = FALSE)
#'
#' # With 20 groups and only 2 error degrees of freedom the critical value
#' # comes from the exact distribution, 26.830. Base R's qtukey() is not
#' # accurate this far out; its value would give 31.485:
#' cv_tukey_hsd(alpha_level = .01, df = 2, n_groups = 20)
#' qtukey(.99, nmeans = 20, df = 2) / sqrt(2)
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

cv_tukey_hsd <- function(alpha_level, df, n_groups, verbose = TRUE) {
  if (!missing(alpha_level)) .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("cv_tukey_hsd", c("df", "n_groups"), required = c("df", "n_groups"))
  # Degrees of freedom checks (mirrors cv_t).
  if (missing(df)) .dmar_refuse_argument(
    "You must specify the degrees of freedom (i.e., 'df'), which must be a positive number.",
    arg = "df")
  if (!df > 0) .dmar_refuse_argument(
    "You must specify a positive value for the degrees of freedom.", arg = NULL)

  # Groups check (specific to Tukey HSD).
  if (missing(n_groups)) .dmar_refuse_argument(
    "You must specify the number of groups (i.e., 'n_groups').",
    arg = "n_groups")
  if (n_groups < 2) .dmar_refuse_argument(
    "This function is appropriate for situations in which there are two or more groups.",
    arg = NULL)

  # Alpha check. For the Studentized range distribution, alpha_level is not split between tails
  # (see Details), so only a single alpha_level argument is needed.
  if (missing(alpha_level)) .dmar_refuse_argument(
    "You must specify 'alpha_level'.", arg = "alpha_level")
  if (alpha_level <= 0 || alpha_level >= 1) .dmar_refuse_argument(
    "Specify 'alpha_level' to be greater than zero and less than 1.",
    arg = "alpha_level")

  # Critical value on the Studentized range scale: stats::qtukey() where it is
  # accurate, the exact quantile elsewhere (R/studentized_range_internals.R).
  q_critical <- .sr_quantile(1 - alpha_level, n_groups, df)

  # Rescale to the scale used for pairwise mean comparisons.
  value <- q_critical / sqrt(2)

  # Areas on the Studentized range scale (probabilities are invariant under the monotonic
  # rescaling by 1/sqrt(2), so these also correspond to the reported value).
  area_less    <- .sr_cdf(q_critical, n_groups, df, lower_tail = TRUE)
  area_greater <- .sr_cdf(q_critical, n_groups, df, lower_tail = FALSE)

  term <- "upper_cv"

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term = term, value = value, area_less = area_less, area_greater = area_greater)))
  }
  if (verbose == FALSE) {
    return(.as_dmar_tbl(data.frame(term = term, value = value)))
  }
}
