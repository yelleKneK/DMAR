# Critical value for the Studentized Maximum Modulus (SMM) distribution.

# Internal: the upper-alpha_level SMM critical value. The m independent |t_i| share a
# single chi/df scale S, so the joint CDF is the one-dimensional integral
#   P(max|Z_i|/S <= c) = int_0^inf (2*pnorm(c*s) - 1)^m f_S(s) ds,
# evaluated over u ~ chi^2_df (with S = sqrt(u/df)) for numerical stability and
# inverted by uniroot.  Deterministic; no Monte Carlo.
.smm_upper_cv <- function(alpha_level, df, m) {
  if (m == 1)
    return(if (is.infinite(df)) stats::qnorm(1 - alpha_level / 2)
           else                 stats::qt(1 - alpha_level / 2, df))
  if (is.infinite(df))
    return(stats::qnorm((1 + (1 - alpha_level)^(1 / m)) / 2))
  target <- 1 - alpha_level
  # Integrate over the chi squared error variate. The limits are the extreme
  # quantiles of that variate rather than (0, Inf): for a large df the density
  # is a narrow bump far from the origin, and an adaptive rule asked to cover
  # the whole half-line samples its way past the mass and returns zero (at
  # df = 200 it does exactly that), which leaves the root finder with no
  # bracket. Placing the limits on the mass makes the rule reliable at every df.
  u_lo <- stats::qchisq(1e-12, df)
  u_hi <- stats::qchisq(1e-12, df, lower.tail = FALSE)
  cdf <- function(c) stats::integrate(function(u)
      (2 * stats::pnorm(c * sqrt(u / df)) - 1)^m * stats::dchisq(u, df),
      lower = u_lo, upper = u_hi, rel.tol = 1e-10)$value
  # The SMM value is at least the two-sided single-t value (equality at m = 1);
  # use that as a lower bracket and grow the upper end until the CDF clears 1-alpha_level.
  lo <- stats::qt(1 - alpha_level / 2, df)
  hi <- lo
  while (cdf(hi) < target && hi < 1e4) hi <- hi * 1.5
  stats::uniroot(function(c) cdf(c) - target, c(lo, hi), tol = 1e-10)$root
}

#' Provides the Critical Value of the Studentized Maximum Modulus Distribution
#'
#' @param alpha_level Type I error rate (i.e., the family-wise false-positive rate).
#' @param df Error degrees of freedom (typically \eqn{N - k} in a one-way ANOVA).
#'   May be \code{Inf}, in which case the known-variance (normal) limit is used.
#' @param n_comparisons The number of simultaneous comparisons (\eqn{m}) for
#'   which a family-wise critical value is desired.
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value as a \code{data.frame}, following
#'   the format used by \code{\link{cv_t}}.
#'
#' @details The Studentized maximum modulus (SMM) distribution is the
#'   distribution of \deqn{\max_{i=1,\ldots,m}\, |Z_i| \,/\, S,}
#'   where \eqn{Z_1, \ldots, Z_m} are independent standard normal variates
#'   and \eqn{S} is independent of the \eqn{Z_i} and equal to
#'   \eqn{\sqrt{\chi^2_{df} / df}}.
#'
#'   \strong{What the SMM is used for.} The SMM critical value is the
#'   multiplier that turns a set of \eqn{m} individual estimates into a family
#'   of simultaneous confidence intervals, or equivalently a family of tests,
#'   that jointly control the family-wise error rate at level \code{alpha_level}.
#'   Because the modulus is the largest of \eqn{m} standardized statistics in
#'   absolute value, requiring that maximum to clear the critical value bounds
#'   the chance of any one of the \eqn{m} intervals failing to cover (or any
#'   one of the \eqn{m} tests producing a false positive). Maxwell, Delaney,
#'   and Kelley (2027, Chapter 5) describe this use in the context of the
#'   multiple-comparisons problem: when a researcher forms several means or
#'   contrasts and wants the stated coverage to hold across the whole set
#'   rather than one interval at a time, the SMM supplies the simultaneous
#'   critical value. A pair-by-pair construction, applied to \eqn{m}
#'   comparisons, is \eqn{\hat\psi_i \pm c_{\alpha;m,df}\,\mathit{SE}_{\hat\psi_i}}
#'   with \eqn{c_{\alpha;m,df}} the SMM critical value returned here. When
#'   \eqn{m = 1} it reduces to the ordinary two-sided \emph{t} critical value.
#'
#'   \strong{How it is computed.} The \eqn{m} statistics are independent given
#'   the common scale estimate \eqn{S}, so the joint distribution factorizes
#'   after conditioning on \eqn{S} and the \eqn{m}-dimensional integral that
#'   defines the maximum modulus collapses to a single one-dimensional integral,
#'   \deqn{P\!\left(\max_i |Z_i|/S \le c\right) =
#'         \int_0^\infty \bigl[\,2\,\Phi(c\,s) - 1\,\bigr]^{m}\, f_S(s)\; ds,}
#'   where \eqn{f_S} is the density of \eqn{S = \sqrt{\chi^2_{df}/df}}. This
#'   package evaluates that integral with \code{\link[stats]{integrate}} and
#'   inverts it with \code{\link[stats]{uniroot}}, so the returned value is
#'   deterministic and accurate to the solver tolerance (there is no Monte
#'   Carlo simulation and no random seed). When \code{df = Inf} the scale is
#'   degenerate at 1 and the closed form
#'   \eqn{c = \Phi^{-1}\!\bigl((1 + (1-\alpha)^{1/m})/2\bigr)} is returned; when
#'   \eqn{m = 1} the value is exactly \eqn{t_{1-\alpha/2,\,df}}.
#'
#' @references
#' Stoline, M. R., & Ury, H. K. (1979). Tables of the Studentized
#'   maximum modulus distribution and an application to multiple
#'   comparisons among means. \emph{Technometrics, 21}(1), 87--93.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 5 on the
#'   multiple-comparisons problem, where simultaneous confidence intervals for
#'   several means or contrasts are developed; Appendix Table A.5 reports SMM
#'   critical values.)
#'
#' @examples
#' # Following the simultaneous-intervals setting of Maxwell, Delaney, and
#' # Kelley (2027, Chapter 5): a researcher forms m = 5 comparisons and wants
#' # all five confidence intervals to hold simultaneously at the .05 level,
#' # with 36 error degrees of freedom.
#' cv_smm(alpha_level = .05, df = 36, n_comparisons = 5)
#'
#' # When m = 1, the SMM critical value reduces to the two-sided
#' # t critical value:
#' cv_smm(alpha_level = .05, df = 36, n_comparisons = 1)$value
#' cv_t(alpha_level = .05, df = 36)$value[2]   # upper_cv from cv_t
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_t}}, \code{\link{cv_tukey_hsd}},
#'   \code{\link{cv_dunnett}}, \code{\link{cv_scheffe}}
#'
#' @keywords design htest
#'
#' @family critical values
#'
#' @export
#' @import stats
cv_smm <- function(alpha_level, df, n_comparisons, verbose = TRUE) {
  if (missing(alpha_level))         stop("You must specify 'alpha_level'.")
  if (missing(df))            stop("You must specify the error degrees of freedom 'df'.")
  if (missing(n_comparisons)) stop("You must specify the number of comparisons 'n_comparisons'.")
  if (alpha_level <= 0 || alpha_level >= 1) stop("'alpha_level' must be in (0, 1).")
  if (df <= 0)                stop("'df' must be positive.")
  if (n_comparisons < 1 || n_comparisons != round(n_comparisons)) {
    stop("'n_comparisons' must be a positive integer.")
  }

  value <- .smm_upper_cv(alpha_level, df, n_comparisons)

  area_less    <- 1 - alpha_level
  area_greater <- alpha_level

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term = "upper_cv", value = value,
                      area_less = area_less, area_greater = area_greater)))
  }
  .as_dmar_tbl(data.frame(term = "upper_cv", value = value))
}
