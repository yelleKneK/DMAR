# Confidence interval for eta squared (effect size for ANOVA).
#' Confidence Interval for Eta Squared (Effect Size for ANOVA)
#'
#' Computes the point estimate and an exact, noncentrality-based confidence
#' interval for the population eta squared (\eqn{\eta^2}), the proportion of
#' variance in the dependent variable accounted for by a fixed effect. Accepts
#' either the raw ANOVA summary (\emph{F}, effect df, error df, total \emph{N})
#' or a fitted model object. Supports both between-subjects designs
#' (\code{\link[stats]{aov}} / \code{\link[stats]{lm}}) and within-subjects /
#' mixed designs (\code{aovlist} fits with an \code{Error()} term in the
#' formula). For factorial and within-subjects designs the function returns
#' one row per effect with the CI for \emph{partial} \eqn{\eta^2} computed
#' against that effect's own error stratum.
#'
#' @param object Optional. A fitted model object of class
#'   \code{\link[stats]{aov}}, \code{\link[stats]{lm}}, or \code{aovlist}
#'   (multi-stratum aov fit, e.g.\
#'   \code{aov(y ~ A + Error(subject/A), data = d)}). For multi-stratum fits
#'   the function walks every error stratum and returns one row per
#'   non-\code{Residuals} effect, identifying the stratum used.
#' @param F_value Observed \emph{F}-value (ignored if \code{object} is supplied).
#' @param df_effect Numerator degrees of freedom for the effect
#'   (ignored if \code{object} is supplied).
#' @param df_error Error (residual) degrees of freedom
#'   (ignored if \code{object} is supplied).
#' @param N Total sample size, the total number of observations (ignored
#'   if \code{object} is supplied; derived automatically from a fitted
#'   model, via \code{\link[stats]{nobs}(object)} for single-stratum
#'   fits, or, for \code{aovlist} fits, recovered as one more than the
#'   sum of the effect and residual degrees of freedom across every
#'   stratum, which is the total number of observations).
#' @param conf_level Desired confidence coverage; default \code{0.95}. Used only
#'   when \code{alpha_lower} and \code{alpha_upper} are both \code{NULL}.
#' @param alpha_lower,alpha_upper Optional Type I error on the lower and upper
#'   side. If both are \code{NULL}, a symmetric interval at \code{conf_level} is
#'   used. If both are supplied, \code{conf_level} is recomputed as
#'   \code{1 - alpha_lower - alpha_upper}.
#'
#' @return A \code{data.frame} with one row per effect. Single-stratum
#'   fits and the raw interface return columns \code{effect},
#'   \code{eta_squared}, \code{lower_limit}, \code{upper_limit},
#'   \code{F_value}, \code{df_effect}, \code{df_error}, \code{N}.
#'   \code{aovlist} (within-subjects / mixed) fits additionally include a
#'   \code{stratum} column. With the raw-argument interface \code{effect}
#'   is \code{"overall"}.
#'
#' @details
#' \strong{Point estimate.} \eqn{\hat{\eta}^2 = df_{\text{effect}} \cdot F /
#' (df_{\text{effect}} \cdot F + df_{\text{error}})}, which equals
#' \eqn{\mathit{SS}_{\text{effect}}/(\mathit{SS}_{\text{effect}} +
#' \mathit{SS}_{\text{error}})}. In a one-way ANOVA this is also
#' \eqn{\mathit{SS}_{\text{effect}}/\mathit{SS}_{\text{total}}}. In a factorial
#' design the same expression gives the per-effect \emph{partial} \eqn{\eta^2}.
#'
#' \strong{Confidence interval.} The CI is constructed by Steiger's (2004)
#' confidence interval transformation principle: a CI for the noncentrality
#' parameter \eqn{\lambda} of the \emph{F} distribution is obtained (via
#' \code{\link{conf_limits_ncf}}) and then mapped through
#' \deqn{\eta^2_{\text{bound}} = \frac{\lambda_{\text{bound}}}{\lambda_{\text{bound}} + N}.}
#' This is the same transformation used by \code{\link{ci_pvaf}} and
#' \code{\link{ci_omega_squared}}; the three functions share CI machinery and
#' differ only in their sample point estimators. When the lower CI on
#' \eqn{\lambda} is not identified (i.e., the observed \emph{F} is below the
#' one-sided critical value), the lower limit on \eqn{\eta^2} is set to 0.
#'
#' \strong{Designs supported.}
#' \itemize{
#'   \item \emph{Between-subjects ANOVA}: fitted \code{aov}/\code{lm}, or
#'     raw \emph{F}/df/N.
#'   \item \emph{Within-subjects or mixed ANOVA}: fitted \code{aovlist}.
#'     Each effect's CI is built from its own stratum's \emph{F} and
#'     residual \emph{df}; \code{N} is the total number of observations
#'     across all strata. The reported \code{stratum} column identifies
#'     which error term each row used.
#' }
#'
#' \strong{Sums of squares in factorial designs.} When a fitted model is
#' supplied, \emph{F}-values are read from \code{anova()} (single-stratum)
#' or \code{summary()} (multi-stratum), both of which use Type I
#' (sequential) sums of squares in base \R. For balanced designs Types
#' I, II, and III agree; for unbalanced designs they differ. If Type II
#' or III \emph{F}-values are required, compute them with e.g.\
#' \code{car::Anova(object, type = 3)} and pass the relevant \emph{F},
#' degrees of freedom, and \emph{N} into the raw-argument interface.
#'
#' @references
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
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
#' @seealso \code{\link{eta_squared}}, \code{\link{ci_eta_squared_partial}},
#'   \code{\link{ci_omega_squared}}, \code{\link{ci_pvaf}},
#'   \code{\link{conf_limits_ncf}}
#'
#' @examples
#' # 1. Raw-argument interface. Bargman's (1970) example.
#' ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#'
#' # Same example with a 90% confidence interval.
#' ci_eta_squared(
#'   F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
#'   conf_level = 0.90
#' )
#'
#' # 2. One way ANOVA from a fitted model: mean IQ gain differs across
#' #        the six grades of the pygmalion data (N = 310).
#' fit_one <- aov(iq_gain ~ factor(grade), data = pygmalion)
#' ci_eta_squared(fit_one)
#'
#' # 3. Two-factor ANOVA: partial eta squared per effect for the
#' #        manipulated expectancy treatment and the measured grade
#' #        classification (pygmalion data, N = 310). The treatment by
#' #        grade interaction is weak here (F = 1.19), so the additive
#' #        model is used.
#' fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
#' ci_eta_squared(fit_additive)
#'
#' # 4. Within-subjects ANOVA. CI computed against the within-subjects
#' #        error stratum (the row reports which one via 'stratum').
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
#' ci_eta_squared(fit_rm)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family confidence intervals for effect sizes
#'
#' @export
#' @import stats

ci_eta_squared <- function(
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
                                    value_name = "eta_squared")))
  }
  if (is.null(F_value) || is.null(df_effect) || is.null(df_error) || is.null(N)) {
    stop("Either provide 'object' (a fitted aov or lm) or all of 'F_value', 'df_effect', 'df_error', and 'N'.")
  }
  .ci_eta_squared_one(F_value, df_effect, df_error, N,
                         alpha_lower, alpha_upper,
                         effect_label = "overall",
                         value_name   = "eta_squared")
}


# Alpha / confidence-level plumbing (shared with the partial variant).
.ci_eta_squared_alphas <- function(conf_level, alpha_lower, alpha_upper) {
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
  list(alpha_lower = alpha_lower, alpha_upper = alpha_upper)
}


# Single-effect worker (point estimate + CI).
.ci_eta_squared_one <- function(F_value, df_effect, df_error, N,
                                   alpha_lower, alpha_upper,
                                   effect_label, value_name) {
  if (!is.numeric(F_value) || length(F_value) != 1L || F_value < 0)
    stop("'F_value' must be a single non-negative number.")
  if (!is.numeric(df_effect) || length(df_effect) != 1L || df_effect <= 0)
    stop("'df_effect' must be a single positive number.")
  if (!is.numeric(df_error) || length(df_error) != 1L || df_error <= 0)
    stop("'df_error' must be a single positive number.")
  if (!is.numeric(N) || length(N) != 1L || N <= df_effect + df_error)
    stop("'N' must be greater than df_effect + df_error.")

  num <- df_effect * F_value
  eta_sq <- num / (num + df_error)

  partial <- value_name == "eta_squared_partial"
  ncp_lims <- .conf_limits_ncf_for(
    caller      = if (partial) "ci_eta_squared_partial" else "ci_eta_squared",
    quantity    = if (partial) "partial eta squared" else "eta squared",
    F_value     = F_value,
    df_1        = df_effect,
    df_2        = df_error,
    alpha_lower = alpha_lower,
    alpha_upper = alpha_upper,
    conf_level  = NULL
  )

  lo_ncp <- ncp_lims$value[ncp_lims$term == "lower_limit"]
  up_ncp <- ncp_lims$value[ncp_lims$term == "upper_limit"]

  # Steiger (2004) Eq. 37; matches ci_pvaf / ci_omega_squared.
  lower_limit <- lo_ncp / (lo_ncp + N)
  upper_limit <- if (is.infinite(up_ncp)) 1 else up_ncp / (up_ncp + N)

  out <- data.frame(
    effect      = effect_label,
    placeholder = eta_sq,
    lower_limit = lower_limit,
    upper_limit = upper_limit,
    F_value     = F_value,
    df_effect   = df_effect,
    df_error    = df_error,
    N           = N,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )
  names(out)[2] <- value_name
  class(out) <- c("dmar_ci_anova", "dmar_tbl", "data.frame")
  out
}


# Model-based worker: dispatches on aovlist (multi-stratum) vs aov/lm.
.ci_eta_squared_from_model <- function(object, alpha_lower, alpha_upper, value_name) {
  if (inherits(object, "aovlist")) {
    return(.ci_eta_squared_from_aovlist(object, alpha_lower, alpha_upper, value_name))
  }
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be an aov, lm, or aovlist fit.")
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
    if (is.na(F_val) || F_val < 0) return(NULL)
    .ci_eta_squared_one(F_val, df_e, df_error, N,
                           alpha_lower, alpha_upper,
                           effect_label = effect,
                           value_name   = value_name)
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


# Within-subjects (multi-stratum) worker: each effect uses its own stratum's
# error df / SS in the Steiger transformation, so the CI is for partial
# eta squared per effect.
.ci_eta_squared_from_aovlist <- function(object, alpha_lower, alpha_upper, value_name) {
  effects <- .aovlist_effects_table(object)
  N <- .aovlist_nobs(object)
  results <- lapply(seq_len(nrow(effects)), function(i) {
    row <- effects[i, , drop = FALSE]
    if (is.na(row$F_value) || row$F_value < 0) return(NULL)
    out <- .ci_eta_squared_one(row$F_value, row$df_effect, row$df_error, N,
                               alpha_lower, alpha_upper,
                               effect_label = row$effect,
                               value_name   = value_name)
    out$stratum <- row$stratum
    class(out) <- c("dmar_ci_anova", "dmar_tbl", "data.frame")
    out
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("No effects with a computable F statistic were found in the aovlist.")
  }
  out <- do.call(rbind, results)
  out <- out[, c("effect", value_name, "lower_limit", "upper_limit",
                 "stratum", "F_value", "df_effect", "df_error", "N")]
  rownames(out) <- NULL
  class(out) <- c("dmar_ci_anova", "dmar_tbl", "data.frame")
  out
}
