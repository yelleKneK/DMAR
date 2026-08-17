# Critical value(s) for the Scheffe procedure for any-contrast simultaneous tests.
#' Provides the Critical Value for the Scheffé Procedure
#'
#' @param alpha_level Type I error rate (i.e., the family-wise false-positive rate).
#' @param df_numerator The numerator degrees of freedom (typically the number
#'   of groups minus 1, \eqn{k - 1}, in a one-way ANOVA).
#' @param df_denominator The denominator (error) degrees of freedom (typically
#'   \eqn{N - k} in a one-way ANOVA).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value in a output style (a
#'   \code{data.frame} with one row per critical value, following the format
#'   used by \code{\link{cv_t}} and \code{\link{cv_tukey_hsd}}).
#'
#' @details The Scheffé critical value protects the family-wise error rate
#'   for the simultaneous test of \emph{any} contrast (or family of
#'   contrasts) in a fixed-effects ANOVA, including data-driven contrasts
#'   selected after looking at the data. It is therefore the most
#'   conservative of the standard procedures.
#'
#'   The critical value, on the scale of a \emph{t}-statistic, is
#'   \deqn{t_{\mathrm{crit}}^{\mathrm{Scheffe}} = \sqrt{(k-1)\, F_{1-\alpha,\,k-1,\,df_{\mathrm{denominator}}}},}
#'   so that a contrast \eqn{\hat\psi} with standard error \eqn{\mathit{SE}_{\hat\psi}} is
#'   declared significant when
#'   \eqn{|\hat\psi/\mathit{SE}_{\hat\psi}| > t_{\mathrm{crit}}^{\mathrm{Scheffe}}}.
#'   The corresponding simultaneous confidence interval is
#'   \eqn{\hat\psi \pm t_{\mathrm{crit}}^{\mathrm{Scheffe}} \cdot \mathit{SE}_{\hat\psi}}.
#'
#'   Like the Studentized range distribution underlying Tukey HSD, the
#'   Scheffé reference is one-sided (the underlying \emph{F} statistic is
#'   non-negative), so \code{alpha_level} is \emph{not} split between two tails.
#'
#'   The Scheffé critical value is a function of a univariate \emph{F}
#'   quantile, which base R supplies through \code{\link[stats]{qf}}, so unlike
#'   \code{\link{cv_dunnett}} and \code{\link{cv_smm}} this function needs no
#'   multivariate distribution machinery. Scheffé's procedure earns its
#'   simultaneous protection over the infinite family of all possible
#'   contrasts by projecting onto the overall \emph{F} test rather than by
#'   integrating a multivariate \emph{t} density; that is why a single
#'   univariate quantile suffices and the \pkg{mvtnorm} package is not
#'   required here. Maxwell, Delaney, and Kelley (2027, Chapter 5) develop the
#'   Scheffé method as the procedure for arbitrary contrasts within the
#'   multiple-comparisons problem.
#'
#' @references
#' Scheffe, H. (1953). A method for judging all contrasts in the analysis
#'   of variance. \emph{Biometrika, 40}, 87--104.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 5 on the
#'   multiple-comparisons problem, where the Scheffé method for arbitrary
#'   contrasts is developed.)
#'
#' @examples
#' # Following the arbitrary-contrasts setting of Maxwell, Delaney, and Kelley
#' # (2027, Chapter 5): a one-way ANOVA with k = 4 groups and 36 error degrees
#' # of freedom (e.g., n = 10 per group). The Scheffé critical value protects
#' # the family-wise error rate over any contrast, including contrasts chosen
#' # after looking at the data.
#' cv_scheffe(alpha_level = .05, df_numerator = 3, df_denominator = 36)
#'
#' # The price of that protection is a larger multiplier than the unadjusted
#' # t critical value at the same error degrees of freedom:
#' cv_t(alpha_level = .05, df = 36)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_t}}, \code{\link{cv_tukey_hsd}},
#'   \code{\link{contrast_test}}
#'
#' @keywords design htest
#'
#' @family critical values
#'
#' @export
#' @import stats
cv_scheffe <- function(alpha_level, df_numerator, df_denominator, verbose = TRUE) {
  if (missing(alpha_level)) stop("You must specify 'alpha_level'.")
  if (alpha_level <= 0 || alpha_level >= 1) stop("'alpha_level' must be in (0, 1).")
  if (missing(df_numerator) || !is.numeric(df_numerator) ||
      df_numerator < 1) {
    stop("'df_numerator' must be a positive integer (e.g., k - 1 in a one-way ANOVA).")
  }
  if (missing(df_denominator) || !is.numeric(df_denominator) ||
      df_denominator <= 0) {
    stop("'df_denominator' must be a positive number.")
  }

  f_critical <- stats::qf(1 - alpha_level, df1 = df_numerator, df2 = df_denominator)
  value      <- sqrt(df_numerator * f_critical)
  area_less    <- stats::pf(f_critical, df1 = df_numerator,
                            df2 = df_denominator, lower.tail = TRUE)
  area_greater <- stats::pf(f_critical, df1 = df_numerator,
                            df2 = df_denominator, lower.tail = FALSE)

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term = "upper_cv", value = value,
                      area_less = area_less, area_greater = area_greater)))
  }
  .as_dmar_tbl(data.frame(term = "upper_cv", value = value))
}
