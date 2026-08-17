# Critical value for Dunnett's procedure for comparing treatments to a control.
# The numeric engine (.dunnett_upper_cv / .dunnett_cdf) lives in
# R/dunnett_internals.R, shared with ci_dunnett().

#' Provides the Critical Value for Dunnett's Many-to-One Comparisons Procedure
#'
#' @param alpha_level Type I error rate (i.e., the family-wise false-positive rate).
#' @param df Error degrees of freedom (typically \eqn{N - k}, where \eqn{k} is
#'   the total number of groups including the control). May be \code{Inf}, in
#'   which case the known-variance (normal) limit is used.
#' @param n_comparisons The number of treatment-versus-control comparisons
#'   (i.e., \eqn{k - 1} where \eqn{k} is the total number of groups).
#' @param alternative The form of the alternative hypothesis: one of
#'   \code{"not_equal"} (two-sided; default), \code{"greater"} (treatments
#'   greater than control), or \code{"less"} (treatments less than control).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value in a output style (a
#'   \code{data.frame} following the format used by \code{\link{cv_t}}).
#'
#' @details Dunnett's procedure controls the family-wise error rate for the
#'   special case of comparing each of \eqn{k - 1} treatments to a single
#'   control (the many-to-one comparisons setting), using a multivariate
#'   \emph{t} reference with constant pairwise correlation \eqn{1/2}. That
#'   correlation is exact for a balanced design: each comparison is
#'   \eqn{(\bar Y_i - \bar Y_0)}, all sharing the one control mean \eqn{\bar Y_0}
#'   of variance \eqn{\sigma^2/n}, while each comparison has variance
#'   \eqn{2\sigma^2/n}, so any two comparisons correlate
#'   \eqn{(\sigma^2/n)/(2\sigma^2/n) = 1/2}. Maxwell, Delaney, and Kelley (2027,
#'   Chapter 5) develop the many-to-one comparisons problem and tabulate these
#'   critical values.
#'
#'   \strong{How it is computed.} The critical value is a quantile of a
#'   \eqn{(k-1)}-dimensional multivariate \emph{t} distribution whose
#'   correlation matrix has a unit diagonal and off-diagonal entries of 1/2.
#'   Because that correlation is a single common value, the comparisons have the
#'   one-factor representation \eqn{Z_i = \sqrt{1/2}\,W + \sqrt{1/2}\,U_i} with a
#'   shared factor \eqn{W} and independent \eqn{U_i}, all standard normal.
#'   Conditioning on \eqn{W} and on the common scale estimate
#'   \eqn{S = \sqrt{\chi^2_{df}/df}} makes the comparisons independent, so the
#'   \eqn{(k-1)}-dimensional integral collapses to two nested one-dimensional
#'   integrals:
#'   \deqn{P\!\left(\max_i T_i \le d\right) = \int_0^\infty\!\!\int_{-\infty}^{\infty}
#'         \Bigl[\Phi\bigl((d\,s - \sqrt{1/2}\,w)/\sqrt{1/2}\bigr)\Bigr]^{k-1}
#'         \phi(w)\,dw\; f_S(s)\,ds}
#'   for the one-sided value (the two-sided value replaces the bracket with the
#'   probability that \eqn{|T_i| \le d}). This package evaluates those integrals
#'   with \code{\link[stats]{integrate}} and inverts with
#'   \code{\link[stats]{uniroot}}, so the returned value is deterministic and
#'   accurate to the solver tolerance, with no Monte Carlo simulation and no
#'   random seed. When \eqn{k - 1 = 1} the value is the ordinary one- or
#'   two-sided \emph{t} critical value.
#'
#' @references
#' Dunnett, C. W. (1955). A multiple comparison procedure for comparing
#'   several treatments with a control. \emph{Journal of the American
#'   Statistical Association, 50}(272), 1096--1121.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 5 on the
#'   multiple-comparisons problem; Dunnett's many-to-one comparisons procedure
#'   is developed there, and Appendix Tables A.6 and A.7 report its critical
#'   values for two- and one-tailed tests.)
#'
#' @examples
#' # Following the many-to-one comparisons setting of Maxwell, Delaney, and
#' # Kelley (2027, Chapter 5): three treatment groups each compared to a single
#' # control (k = 4 groups, so 3 comparisons) with 36 error degrees of freedom.
#' # The two-sided critical value controls the family-wise error rate at .05.
#' cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3)
#'
#' # When the treatments are expected only to exceed the control, a one-sided
#' # ("treatments greater than control") critical value is smaller.
#' cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3, alternative = "greater")
#'
#' @note The constant-correlation assumption holds for balanced designs (equal
#'   n per group); for severely unbalanced designs use
#'   \code{\link[multcomp]{glht}} instead.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_t}}, \code{\link{cv_tukey_hsd}},
#'   \code{\link{cv_smm}}, \code{\link{cv_scheffe}},
#'   \code{\link{contrast_test}}
#'
#' @keywords design htest
#'
#' @family critical values
#'
#' @export
#' @import stats
cv_dunnett <- function(alpha_level, df, n_comparisons,
                       alternative = "not_equal", verbose = TRUE) {
  if (missing(alpha_level))         stop("You must specify 'alpha_level'.")
  if (missing(df))            stop("You must specify the error degrees of freedom 'df'.")
  if (missing(n_comparisons)) stop("You must specify 'n_comparisons' (= number of treatments compared to control).")
  if (alpha_level <= 0 || alpha_level >= 1) stop("'alpha_level' must be in (0, 1).")
  if (df <= 0)                stop("'df' must be positive.")
  if (n_comparisons < 1 || n_comparisons != round(n_comparisons)) {
    stop("'n_comparisons' must be a positive integer.")
  }

  alt <- tolower(alternative)
  alt_map <- list(
    "not_equal"   = c("ne", "2s", "not equal", "two.sided", "two sided",
                       "two-sided", "!=", "not_equal"),
    "greater"     = c("greater than", "greater-than", "greater",
                       "greater.than", "gt", "g", ">", ">="),
    "less"        = c("less than", "less-than", "less", "lesser",
                       "less.than", "lt", "l", "<", "<=")
  )
  alt_canonical <- NA_character_
  for (key in names(alt_map)) {
    if (alt %in% alt_map[[key]]) {
      alt_canonical <- key
      break
    }
  }
  if (is.na(alt_canonical)) {
    stop("'alternative' must be one of \"not.equal\", \"greater\", or \"less\".")
  }

  two_sided <- alt_canonical == "not_equal"
  value <- .dunnett_upper_cv(alpha_level, df, n_comparisons, two_sided)

  if (two_sided) {
    term <- "upper_cv"
  } else {
    if (alt_canonical == "less") value <- -value
    term <- if (alt_canonical == "greater") "upper_cv" else "lower_cv"
  }
  area_less    <- 1 - alpha_level
  area_greater <- alpha_level

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term = term, value = value,
                      area_less = area_less, area_greater = area_greater)))
  }
  .as_dmar_tbl(data.frame(term = term, value = value))
}
