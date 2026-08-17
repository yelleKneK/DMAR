# Confidence interval for omega squared (ANOVA effect size).
#' Confidence Interval for Omega Squared (Effect Size for ANOVA)
#'
#' Computes the point estimate and an exact, noncentrality-based confidence
#' interval for the population omega squared (\eqn{\omega^2}), the proportion of
#' variance in the dependent variable accounted for by a fixed effect. Accepts
#' either the raw ANOVA summary (\emph{F}, effect df, error df, total \emph{N})
#' or a fitted \code{\link[stats]{aov}} / \code{\link[stats]{lm}} object, in
#' which case the function returns a row per effect (partial \eqn{\omega^2} in
#' factorial designs).
#'
#' @param object Optional. A fitted \code{\link[stats]{aov}} or
#'   \code{\link[stats]{lm}} object. When supplied, the function loops over the
#'   non-\code{Residuals} rows of \code{\link[stats]{anova}(object)} and returns
#'   one row per effect.
#' @param F_value Observed \emph{F}-value from the fixed-effects ANOVA
#'   (ignored if \code{object} is supplied).
#' @param df_effect Numerator degrees of freedom for the effect
#'   (ignored if \code{object} is supplied).
#' @param df_error Error (residual) degrees of freedom
#'   (ignored if \code{object} is supplied).
#' @param N Total sample size (ignored if \code{object} is supplied;
#'   \code{\link[stats]{nobs}(object)} is used instead).
#' @param conf_level Desired confidence coverage; default \code{0.95}. Used only
#'   when \code{alpha_lower} and \code{alpha_upper} are both \code{NULL}.
#' @param alpha_lower,alpha_upper Optional Type I error on the lower and upper
#'   side. If both are \code{NULL}, a symmetric interval at \code{conf_level} is
#'   used. If both are supplied, \code{conf_level} is recomputed as
#'   \code{1 - alpha_lower - alpha_upper}.
#'
#' @return A \code{data.frame} with one row per effect and the columns
#'   \code{effect}, \code{omega_squared} (point estimate), \code{lower_limit},
#'   \code{upper_limit}, \code{F_value}, \code{df_effect}, \code{df_error}, and
#'   \code{N}. When the raw-argument interface is used, \code{effect} is
#'   \code{"overall"}.
#'
#' @details
#' \strong{Point estimate.} The function reports the usual sample omega squared,
#' which for a one-way design can be written as
#' \deqn{\hat{\omega}^2 = \frac{\mathit{SS}_{\text{effect}} - df_{\text{effect}} \cdot \mathit{MS}_{\text{error}}}{\mathit{SS}_{\text{total}} + \mathit{MS}_{\text{error}}} = \frac{df_{\text{effect}} (F - 1)}{df_{\text{effect}} (F - 1) + N}}
#' (Hays, 1994; Keppel, 1991). For factorial designs the same formula applied per
#' effect yields \emph{partial} omega squared (Olejnik & Algina, 2003); values
#' below zero are truncated to zero.
#'
#' \strong{Confidence interval.} The CI is constructed by Steiger's (2004,
#' Proposition 1) confidence interval transformation principle: a CI for the
#' noncentrality parameter \eqn{\lambda} of the \emph{F} distribution is
#' obtained (via \code{\link{conf_limits_ncf}}) and then mapped through
#' \deqn{\omega^2_{\text{bound}} = \frac{\lambda_{\text{bound}}}{\lambda_{\text{bound}} + N}.}
#' When the lower CI on \eqn{\lambda} is not identified (i.e., the observed
#' \emph{F} is below the one-sided critical value), the lower limit on
#' \eqn{\omega^2} is set to 0, matching the convention used in
#' \code{\link{ci_pvaf}}. In a one-way design, the
#' interval produced here is identical to the CI for \eqn{\eta^2} from
#' \code{\link{ci_pvaf}}; the two estimands coincide in the population and
#' differ only in their \emph{sample} estimators (an implication of the
#' confidence interval transformation principle of Steiger, 2004).
#'
#' \strong{Sums of squares in factorial designs.} When a fitted model is
#' supplied, the function reads the \emph{F}-values from \code{anova()}, which
#' in base \R uses Type I (sequential) sums of squares. For balanced designs,
#' Types I, II, and III give identical \emph{F}-values; for unbalanced designs
#' they differ. If Type II or III \emph{F}-values are required, compute them
#' with e.g.\ \code{car::Anova(object, type = 3)} and pass the relevant
#' \emph{F} / df into the raw-argument interface.
#'
#' @references
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Keppel, G. (1991). \emph{Design and analysis: A researcher's handbook} (3rd
#'   ed.). Prentice Hall.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
#'   statistics: Measures of effect size for some common research designs.
#'   \emph{Psychological Methods, 8}(4), 434--447.
#'   \doi{10.1037/1082-989X.8.4.434}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{ci_pvaf}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' # 1. Raw-argument interface. Bargman's (1970) example, also used in
#' #        Venables (1975), Fleishman (1980), and Steiger (2004): a 5-group
#' #        one-way ANOVA with 11 subjects per group, observed F = 11.221.
#' ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#'
#' # Same example with a 90% confidence interval.
#' ci_omega_squared(
#'   F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
#'   conf_level = 0.90
#' )
#'
#' # 2. One way ANOVA from a fitted model: mean IQ gain differs across
#' #        the six grades of the pygmalion data (N = 310).
#' fit_one <- aov(iq_gain ~ factor(grade), data = pygmalion)
#' ci_omega_squared(fit_one)
#'
#' # 3. Two-factor ANOVA: partial omega squared per effect for the
#' #        manipulated expectancy treatment and the measured grade
#' #        classification (pygmalion data, N = 310). The treatment by
#' #        grade interaction is weak here (F = 1.19), so the additive
#' #        model is used.
#' fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
#' ci_omega_squared(fit_additive)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family confidence intervals for effect sizes
#'
#' @export
#' @import stats

ci_omega_squared <- function(
  object      = NULL,
  F_value     = NULL,
  df_effect   = NULL,
  df_error    = NULL,
  N           = NULL,
  conf_level  = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL
) {
  # Alpha / confidence-level plumbing (parallels ci_pvaf).
  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  } else if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_lower - alpha_upper
  } else if (!is.null(alpha_lower) && is.null(alpha_upper)) {
    stop("'alpha_lower' specified but not 'alpha_upper'.")
  } else {
    stop("'alpha_upper' specified but not 'alpha_lower'.")
  }

  if (alpha_lower < 0 || alpha_lower > 0.5) stop("'alpha_lower' must be in [0, 0.5].")
  if (alpha_upper < 0 || alpha_upper > 0.5) stop("'alpha_upper' must be in [0, 0.5].")
  if (conf_level <= 0 || conf_level >= 1)   stop("'conf_level' must be strictly between 0 and 1.")

  # Dispatch on input style.
  if (!is.null(object)) {
    # One clamp warning per call, not one per clamped effect.
    return(.warn_ncf_clamp_once(
      .ci_omega_squared_from_model(object, alpha_lower, alpha_upper)))
  }

  if (is.null(F_value) || is.null(df_effect) || is.null(df_error) || is.null(N)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', 'df_error', and 'N'.")
  }

  .ci_omega_squared_one(
    F_value, df_effect, df_error, N, alpha_lower, alpha_upper,
    effect_label = "overall"
  )
}


# Single-effect worker: validates inputs, computes the point estimate, pulls
# the NCP CI, and maps it through Steiger's (2004) transformation.
.ci_omega_squared_one <- function(F_value, df_effect, df_error, N,
                                  alpha_lower, alpha_upper, effect_label) {
  if (!is.numeric(F_value) || F_value <= 0)
    stop("'F_value' must be a single positive number.")
  if (!is.numeric(df_effect) || df_effect <= 0)
    stop("'df_effect' must be a positive number.")
  if (!is.numeric(df_error)  || df_error  <= 0)
    stop("'df_error' must be a positive number.")
  if (!is.numeric(N) || N <= df_effect + df_error)
    stop("'N' must be greater than df_effect + df_error.")

  # Point estimate. Bounded below at 0, matching Hays (1994) and common practice.
  num <- df_effect * (F_value - 1)
  omega_sq <- max(0, num / (num + N))

  # Noncentrality-parameter CI.
  ncp_lims <- .conf_limits_ncf_for(
    caller      = "ci_omega_squared",
    quantity    = "omega squared",
    F_value     = F_value,
    df_1        = df_effect,
    df_2        = df_error,
    alpha_lower = alpha_lower,
    alpha_upper = alpha_upper,
    conf_level  = NULL
  )

  lo_ncp <- ncp_lims$value[ncp_lims$term == "lower_limit"]
  up_ncp <- ncp_lims$value[ncp_lims$term == "upper_limit"]

  # Steiger (2004) Eq. 16; his Eq. 17 is this map applied to the worked example.
  lower_limit <- lo_ncp / (lo_ncp + N)
  upper_limit <- if (is.infinite(up_ncp)) 1 else up_ncp / (up_ncp + N)

  out <- data.frame(
    effect        = effect_label,
    omega_squared = omega_sq,
    lower_limit   = lower_limit,
    upper_limit   = upper_limit,
    F_value       = F_value,
    df_effect     = df_effect,
    df_error      = df_error,
    N             = N,
    stringsAsFactors = FALSE,
    row.names     = NULL
  )
  class(out) <- c("dmar_ci_anova", "dmar_tbl", "data.frame")
  out
}


# Model-based worker: loops over the non-Residuals rows of anova(object).
.ci_omega_squared_from_model <- function(object, alpha_lower, alpha_upper) {
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be an aov or lm fit.")
  }

  tbl <- stats::anova(object)

  resid_row <- which(rownames(tbl) == "Residuals")
  if (length(resid_row) == 0L) {
    stop("Could not find a 'Residuals' row in anova(object); is this a standard fixed-effects ANOVA?")
  }

  df_error <- tbl[resid_row, "Df"]
  N <- stats::nobs(object)

  effect_rows <- setdiff(rownames(tbl), "Residuals")

  results <- lapply(effect_rows, function(effect) {
    F_val <- tbl[effect, "F value"]
    df_e  <- tbl[effect, "Df"]
    if (is.na(F_val) || F_val <= 0) return(NULL)
    .ci_omega_squared_one(
      F_val, df_e, df_error, N, alpha_lower, alpha_upper,
      effect_label = effect
    )
  })
  results <- Filter(Negate(is.null), results)

  if (length(results) == 0L) {
    stop("No effects with a computable F statistic were found in the model.")
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  class(out) <- c("dmar_ci_anova", "dmar_tbl", "data.frame")
  out
}
