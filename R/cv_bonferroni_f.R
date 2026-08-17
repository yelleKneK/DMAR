# Critical value for the Bonferroni-adjusted F test of a set of contrasts.
#' Provides the Bonferroni-Adjusted Critical Value for an \emph{F} Test of One of Several Contrasts
#'
#' @param alpha_level The family-wise Type I error rate (i.e., the rate for the set
#'   of \code{n_comparisons} tests taken together). Default \code{0.05}, the
#'   level Maxwell, Delaney, and Kelley (2027) tabulate.
#' @param df_denominator The denominator (error) degrees of freedom (a
#'   positive number). In a one-way design with \eqn{N} observations and
#'   \eqn{a} groups this is \eqn{N - a}.
#' @param n_comparisons The number of comparisons in the family, \eqn{C} (a
#'   positive integer).
#' @param df_numerator The numerator degrees of freedom. Default \code{1},
#'   because a contrast carries a single degree of freedom, which is the case
#'   the Appendix table covers.
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value in a output style (a \code{data.frame}
#'   following the format used by \code{\link{cv_t}}). When
#'   \code{verbose = TRUE} the \code{area_greater} column reports the
#'   per-comparison error rate \eqn{\alpha/C} that the adjustment spends on
#'   each test.
#'
#' @details The Bonferroni adjustment tests each of \eqn{C} contrasts at
#'   \eqn{\alpha/C} rather than \eqn{\alpha}, which holds the family-wise
#'   error rate at or below \eqn{\alpha} whatever the contrasts are and
#'   however they are correlated. The critical value is therefore an ordinary
#'   upper-tail \emph{F} quantile read at the smaller per-comparison rate,
#'   \deqn{F_{\alpha/C;\,\mathrm{df_{num}},\,\mathrm{df_{den}}},}
#'   which is what \code{\link{cv_f}} would return if handed \code{alpha / C}.
#'   This function exists because the adjustment is worth naming: the whole of
#'   it is the division, and seeing \eqn{\alpha/C} reported back in
#'   \code{area_greater} is the point. Maxwell, Delaney, and Kelley (2027)
#'   tabulate these values for one numerator degree of freedom and a
#'   family-wise alpha of .05 in their Appendix Table A.3.
#'
#'   The default \code{df_numerator = 1} covers the case the table addresses
#'   and the one that arises in practice, since a contrast among means is a
#'   single-degree-of-freedom question. Supply a larger value to Bonferroni
#'   adjust a family of multiple-degree-of-freedom model comparisons.
#'
#'   The procedure is often called Dunn's, after Dunn (1961), who first
#'   applied the Bonferroni inequality to multiple contrasts. It is not the
#'   rank-sum procedure of \code{\link{dunn_test}}, which the same author
#'   published three years later.
#'
#'   \strong{When something else is better.} Bonferroni makes no use of the
#'   structure of the family, so a procedure built for a particular structure
#'   beats it there: \code{\link{cv_tukey_hsd}} is more powerful for all
#'   pairwise comparisons, and \code{\link{cv_dunnett}} is more powerful for
#'   comparing several treatments to one control. Bonferroni's advantage is
#'   generality; it applies to any set of contrasts chosen in advance, and it
#'   can beat Tukey's method when only a few of the pairwise comparisons were
#'   planned. Because it is conservative, a step-down variant such as Holm's
#'   is uniformly more powerful while controlling the same rate, and is
#'   available through the \code{method} argument of
#'   \code{\link{contrast_adjusted}} and \code{\link[stats]{p.adjust}}.
#'
#' @references
#' Dunn, O. J. (1961). Multiple comparisons among means. \emph{Journal of the
#'   American Statistical Association, 56}(293), 52--64.
#'   \doi{10.2307/2282330}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 5 on the multiple-comparisons
#'   problem; Appendix Table A.3 reports these critical values.)
#'
#' @examples
#' # Five planned contrasts with 20 error degrees of freedom, holding the
#' # family-wise error rate at .05. The area_greater column shows the .01
#' # per-comparison rate the adjustment spends on each test.
#' cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5)
#'
#' # It is the F critical value read at alpha / C.
#' cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5,
#'                 verbose = FALSE)$value
#' cv_f(alpha_level = .05 / 5, df_numerator = 1, df_denominator = 20,
#'      verbose = FALSE)$value[2]
#'
#' # With one comparison there is nothing to adjust.
#' cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 1,
#'                 verbose = FALSE)$value
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_f}}, \code{\link{cv_tukey_hsd}},
#'   \code{\link{cv_dunnett}}, \code{\link{cv_scheffe}},
#'   \code{\link{contrast_adjusted}}
#'
#' @keywords distribution htest design
#'
#' @family critical values
#'
#' @export
#' @import stats

cv_bonferroni_f <- function(alpha_level = 0.05, df_denominator, n_comparisons,
                            df_numerator = 1, verbose = TRUE) {
  if (missing(df_denominator)) stop("You must specify the denominator degrees of freedom (i.e., 'df_denominator'), which must be a positive number.")
  if (missing(n_comparisons)) stop("You must specify the number of comparisons in the family (i.e., 'n_comparisons').")
  if (!df_denominator > 0) stop("You must specify a positive value for 'df_denominator'.")
  if (!df_numerator > 0) stop("You must specify a positive value for 'df_numerator'.")
  if (n_comparisons < 1 || n_comparisons != round(n_comparisons)) stop("'n_comparisons' must be a positive integer.")
  if (alpha_level <= 0 || alpha_level >= 1) stop("Specify 'alpha_level' to be greater than zero and less than 1.")

  # The whole of the Bonferroni adjustment: spend alpha_level/C on each comparison.
  alpha_per_comparison <- alpha_level / n_comparisons

  term <- "upper_cv"
  value <- qf(p = 1 - alpha_per_comparison, df1 = df_numerator, df2 = df_denominator)
  area_less <- pf(value, df1 = df_numerator, df2 = df_denominator, lower.tail = TRUE)
  area_greater <- pf(value, df1 = df_numerator, df2 = df_denominator, lower.tail = FALSE)

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term, value, area_less, area_greater)))
  }
  .as_dmar_tbl(data.frame(term, value))
}
