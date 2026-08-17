# Confidence interval for partial eta squared (effect size for ANOVA).
#' Confidence Interval for Partial Eta Squared (Effect Size for ANOVA)
#'
#' Computes the point estimate and an exact, noncentrality-based confidence
#' interval for the population \emph{partial} eta squared (\eqn{\eta^2_p}).
#' Accepts either the raw ANOVA summary (\emph{F}, effect df, error df, total
#' \emph{N}) or a fitted \code{aov}/\code{lm}/\code{aovlist} object, in which
#' case the function returns one row per effect (with \code{stratum}
#' identification for within-subjects fits).
#'
#' This is the explicitly-named counterpart of \code{\link{ci_eta_squared}}.
#' The two share point-estimate and CI machinery: in a one-way ANOVA partial
#' \eqn{\eta^2} coincides with \eqn{\eta^2}; in a factorial or within-subjects
#' ANOVA both functions return the per-effect \emph{partial} value computed
#' against that effect's own error stratum. Use \code{ci_eta_squared_partial}
#' when you want the function name to make the partial interpretation explicit.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or \code{aovlist}.
#' @param F_value Observed \emph{F}-value (ignored if \code{object} is supplied).
#' @param df_effect Numerator degrees of freedom for the effect
#'   (ignored if \code{object} is supplied).
#' @param df_error Error (residual) degrees of freedom
#'   (ignored if \code{object} is supplied).
#' @param N Total sample size (ignored if \code{object} is supplied).
#' @param conf_level Desired confidence coverage; default \code{0.95}.
#' @param alpha_lower,alpha_upper Optional Type I error on the lower and upper
#'   side.
#'
#' @return A \code{data.frame} with one row per effect. Single-stratum
#'   fits and the raw interface return columns \code{effect},
#'   \code{eta_squared_partial}, \code{lower_limit}, \code{upper_limit},
#'   \code{F_value}, \code{df_effect}, \code{df_error}, \code{N}.
#'   \code{aovlist} fits additionally include a \code{stratum} column. With
#'   the raw-argument interface \code{effect} is \code{"overall"}.
#'
#' @references
#' Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor ANOVA
#'   designs. \emph{Educational and Psychological Measurement, 33}(1), 107--112.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Smithson, M. (2001). Correct confidence intervals for various regression
#'   effect sizes and parameters: The importance of noncentral distributions in
#'   computing intervals. \emph{Educational and Psychological Measurement, 61},
#'   605--632. \doi{10.1177/00131640121971392}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{eta_squared_partial}}, \code{\link{ci_eta_squared}}
#'
#' @examples
#' # Raw-argument interface.
#' ci_eta_squared_partial(F_value = 11.221, df_effect = 4,
#'                        df_error = 50, N = 55)
#'
#' # Two-factor ANOVA: per-effect partial eta squared with CI for the
#' # manipulated expectancy treatment and the measured grade
#' # classification (pygmalion data, N = 310). The treatment by grade
#' # interaction is weak here (F = 1.19), so the additive model is used.
#' fit <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
#' ci_eta_squared_partial(fit)
#'
#' # Within-subjects ANOVA.
#' set.seed(113)
#' n <- 20
#' rm_data <- data.frame(
#'   subject = factor(rep(seq_len(n), each = 3)),
#'   time    = factor(rep(c("Pre", "Mid", "Post"), n),
#'                    levels = c("Pre", "Mid", "Post")),
#'   y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
#'             0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
#' )
#' fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
#' ci_eta_squared_partial(fit_rm)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_eta_squared_partial <- function(
  object      = NULL,
  F_value     = NULL,
  df_effect   = NULL,
  df_error    = NULL,
  N           = NULL,
  conf_level  = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL
) {
  alphas <- .ci_eta_squared_alphas(conf_level, alpha_lower, alpha_upper)
  alpha_lower <- alphas$alpha_lower
  alpha_upper <- alphas$alpha_upper

  if (!is.null(object)) {
    # One clamp warning per call, not one per clamped effect.
    return(.warn_ncf_clamp_once(
      .ci_eta_squared_from_model(object, alpha_lower, alpha_upper,
                                    value_name = "eta_squared_partial")))
  }
  if (is.null(F_value) || is.null(df_effect) || is.null(df_error) || is.null(N)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', 'df_error', and 'N'.")
  }
  .ci_eta_squared_one(F_value, df_effect, df_error, N,
                         alpha_lower, alpha_upper,
                         effect_label = "overall",
                         value_name   = "eta_squared_partial")
}
