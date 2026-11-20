# =============================================================================
#  Internal helpers for the meta-analysis family
# =============================================================================
#
# Shared by the exported synthesis functions:
#
#   meta_es()   -- pool any effect sizes given (yi, vi)
#   meta_smd()  -- pool standardized mean differences
#   meta_r()    -- pool correlations on the Fisher's Z scale
#
# Nothing in this file is exported. The numerics follow the standard
# random effects model y_i = mu + u_i + e_i with u_i ~ N(0, tau^2) and
# e_i ~ N(0, v_i), v_i known.
#
# Layout:
#   (A) tau^2 estimators: DerSimonian-Laird, Paule-Mandel, REML
#   (B) Q-profile confidence interval for tau^2 (Viechtbauer, 2007)
#   (C) The assembled fit: .meta_fit() returns everything the public
#       functions report (estimate, SE, Hartung-Knapp, heterogeneity,
#       prediction interval)
#   (D) Correlations corrected for attenuation: the Fisher Z effect sizes
#       and sampling variances meta_r() pools

# --- (A) tau^2 estimators ----------------------------------------------------

# Cochran's Q at given weights.
.meta_Q <- function(yi, vi, tau2 = 0) {
  w <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  sum(w * (yi - mu)^2)
}

# DerSimonian and Laird (1986) moment estimator.
.meta_tau2_dl <- function(yi, vi) {
  k <- length(yi)
  w <- 1 / vi
  Q <- .meta_Q(yi, vi)
  C <- sum(w) - sum(w^2) / sum(w)
  max(0, (Q - (k - 1)) / C)
}

# The tau^2 at which the generalized Q statistic equals `target`, and 0 when
# Q at tau^2 = 0 is already at or below the target. The Paule-Mandel
# estimate and both limits of the Q-profile interval are values of this
# kind.
#
# The generalized Q is a number without units, and it decreases in tau^2:
# with w_i = 1 / (v_i + tau^2) and mu the weighted mean, its derivative is
# -sum(w_i^2 (y_i - mu)^2), the term from mu vanishing because
# sum(w_i (y_i - mu)) = 0. The root is therefore unique. The search is made
# on u = tau^2 / s, with s = max(v_i) + var(y_i) as in .meta_tau2_reml(), and
# Q is computed from y_i / sqrt(s) and v_i / s. A search on tau^2 itself,
# stopped at an absolute distance from the root, stopped within the same
# distance whatever the units of the effect sizes: for effect sizes in units
# of 1e-6 that distance, 1.5e-8 at the square root of the machine epsilon,
# was larger than the whole bracket, and an end of the bracket came back as
# the root. On u the root for (c y_i, c^2 v_i) is c^2 times the root for
# (y_i, v_i) up to rounding error.
#
# The test against the boundary is made here, on the same scaled Q the
# search uses. Made on the unscaled data, as it once was, it could disagree
# by rounding error with the scaled Q when Q at 0 sat within a few units of
# rounding of the target, and uniroot() then stopped with "f() values at end
# points not of opposite sign". Now Q at 0 at or below the target gives 0
# and anything above it a root that is positive and, at such an input,
# within rounding error of 0, so the result moves continuously across the
# boundary in every unit.
#
# uniroot() stops when the bracket around the root is narrower than
# 4 eps |u| + tol, eps the machine epsilon (zeroin's test), so a tolerance
# fixed in u, as the machine epsilon to the power 0.75 once was, stops at a
# distance relative to s and not to the root. When one sampling variance
# dwarfs the others, s is far above tau^2 and the root lost accuracy
# relative to itself: an error of 1.7e-2 of the lower limit for twelve
# studies, eleven with variances near 1e-3 and one with variance 1e8, and
# tau^2 near .01. The tolerance is now eps times min(v_i / s), about the
# smallest change in u that can change some v_i / s + u, and with it Q, in
# floating point, so the stopping distance is relative to the root down to
# that resolution. On 2,000 random sets of 2 to 30 studies this took at
# most 17 iterations (14 at the earlier tolerance).
#
# The upper end of the bracket needs no search: with SS the sum of squares
# of y_i / sqrt(s) about their mean, Q at u is at most SS / (min(v_i / s) +
# u), the weighted sum of squares about the weighted mean being at most the
# one about the unweighted mean, and so below the target at u = SS / target.
# The doubling loop guards only against Q there rounding up to the target.
.meta_Q_root <- function(yi, vi, target) {
  s  <- max(vi) + stats::var(yi)
  ys <- yi / sqrt(s)
  vs <- vi / s
  f  <- function(u) .meta_Q(ys, vs, u) - target
  f_lower <- f(0)
  if (f_lower <= 0) return(0)
  upper <- sum((ys - mean(ys))^2) / target
  f_upper <- f(upper)
  while (f_upper > 0) {
    upper <- 2 * upper
    f_upper <- f(upper)
  }
  u <- stats::uniroot(f, c(0, upper), f.lower = f_lower, f.upper = f_upper,
                      tol = .Machine$double.eps * min(vs),
                      maxiter = 1000L)$root
  u * s
}

# Paule and Mandel (1982) estimator: the tau^2 at which the generalized
# Q statistic equals its expectation k - 1, and 0 when Q at tau^2 = 0 is
# already at or below k - 1 (the test is made in .meta_Q_root()).
.meta_tau2_pm <- function(yi, vi) {
  .meta_Q_root(yi, vi, length(yi) - 1)
}

# Restricted maximum likelihood, by direct maximization of the restricted
# log-likelihood profile in tau^2 (one-dimensional, smooth).
#
# The search is made on u = tau^2 / s, with s = max(v_i) + var(y_i), the
# scale of the data in squared units of the effect size, over 0 <= u <= 10.
# optimize() stops at an absolute step in its argument (about a third of
# `tol`), so a search on tau^2 itself stops within the same absolute
# distance of the maximum whatever the units: 3e-05 at the default
# tolerance, the fourth root of the machine epsilon, and 5e-09 at its
# square root. That is small for standardized mean differences and as
# large as tau^2 itself for effect sizes in small units (risk differences
# with variances near 1e-5). On u the stopping distance is relative to s,
# and the estimate for (c y_i, c^2 v_i) is c^2 times the estimate for
# (y_i, v_i). The log-likelihood is written in the scaled quantities, less
# the terms in log(s) that do not depend on u, so the comparison with the
# boundary value is unit free as well.
#
# optimize() stops within 2 (sqrt(eps) |u| + tol / 3) of the maximum, eps
# the machine epsilon. A tolerance fixed in u, as the square root of eps
# once was, therefore stops at a distance relative to s, about 1e-8 of it,
# and not relative to the estimate. When tau^2 is far below s the estimate
# lost accuracy relative to itself: errors up to 8e-4 of tau^2 on 150
# random sets with sampling variances from 1e-3 to 3, and, with one
# sampling variance far above the others (twelve studies, eleven with
# variances near 1e-3, tau^2 near .01), 6e-3 when that variance was 1e4,
# 0.27 at 1e6, and 186 times tau^2 at 1e8, against the root of the REML
# score equation in 200-bit arithmetic. As in .meta_Q_root(), the
# tolerance is now eps times min(v_i / s), so the stopping distance is
# relative to the estimate; on the same sets the largest error became
# 1.4e-6 and 1e-7, the remainder being the flatness of the likelihood at
# its maximum, which limits how closely any search on the likelihood
# itself can place it.
.meta_tau2_reml <- function(yi, vi) {
  s <- max(vi) + stats::var(yi)
  vs <- vi / s
  ll <- function(u) {
    w  <- 1 / (vs + u)
    mu <- sum(w * yi) / sum(w)
    -0.5 * sum(log(vs + u)) - 0.5 * log(sum(w)) -
      0.5 * sum(w * (yi - mu)^2) / s
  }
  opt <- stats::optimize(ll, c(0, 10), maximum = TRUE,
                         tol = .Machine$double.eps * min(vs))
  # The maximum can sit at the boundary; prefer 0 when it fits as well.
  if (ll(0) >= opt$objective - 1e-10) 0 else opt$maximum * s
}

# --- (B) Q-profile confidence interval for tau^2 -----------------------------
#
# Invert the generalized Q statistic (Viechtbauer, 2007): the limits are the
# tau^2 values at which Q_gen(tau2) equals the chi square quantiles with
# k - 1 df. The lower limit is 0 when even tau^2 = 0 fails to push Q above
# the upper quantile, and the upper limit is 0 when Q at tau^2 = 0 is
# already at or below the lower quantile. Both limits, boundary tests
# included, are found on the scale of the data by .meta_Q_root(), so the
# limits for (c y_i, c^2 v_i) are c^2 times the limits for (y_i, v_i), and
# the I^2 limits, which map them through the typical within-study variance,
# are the same.
.meta_tau2_ci <- function(yi, vi, conf_level = 0.95) {
  k <- length(yi)
  alpha <- 1 - conf_level
  q_hi <- stats::qchisq(1 - alpha / 2, df = k - 1)
  q_lo <- stats::qchisq(alpha / 2, df = k - 1)
  c(.meta_Q_root(yi, vi, q_hi), .meta_Q_root(yi, vi, q_lo))
}

# --- (C) the assembled fit ---------------------------------------------------
#
# Returns the named list every public meta_* function formats. tau^2 method
# is one of "reml" (default), "pm", "dl", or "fe" (fixed effect: tau^2 = 0
# by assumption, no prediction interval). hartung_knapp applies the
# Hartung-Knapp-Sidik-Jonkman small-sample adjustment: the pooled SE is
# rescaled from the weighted residuals and the reference distribution is
# t with k - 1 df. The prediction interval follows Higgins, Thompson, and
# Spiegelhalter (2009), with t on k - 2 df, and uses the same SE.
#
# When the effect sizes are identical, exactly or up to rounding error, the
# Hartung-Knapp SE is zero or numerically zero. The numbers are reported as
# computed (as metafor's test = "knha" reports them) with one warning of
# class dmar_meta_hk_zero_se, which states the conventional standard error
# and the normal interval it gives and carries them as fields. The check is
# on the effect sizes themselves, at 8 units of relative rounding error in
# their spread, so that effects that differ in a late digit (0.5 against
# 0.5 + 1e-9) do not raise it.
#
# Under the adjustment the fit also returns q, the weighted residual mean
# square on k - 1 df, and se_conventional, (sum w)^(-1/2); the adjusted SE
# is sqrt(q) times the second. The public functions record both as
# attributes, so that an interval much narrower than the conventional one
# can be seen for what it is.
#
# `transform` and `scale_note` serve the warning only. meta_r() fits on the
# Fisher Z scale and reports the estimate and its limits as correlations,
# so it passes tanh and " (on the Fisher Z scale)": the estimate and the
# two intervals in the message and in the fields are then correlations,
# and the conventional standard error is labeled with the scale it is on.
.meta_fit <- function(yi, vi, method = "reml", hartung_knapp = TRUE,
                      conf_level = 0.95, transform = identity,
                      scale_note = "") {
  k <- length(yi)
  tau2 <- switch(method,
                 fe   = 0,
                 dl   = .meta_tau2_dl(yi, vi),
                 pm   = .meta_tau2_pm(yi, vi),
                 reml = .meta_tau2_reml(yi, vi))
  w   <- 1 / (vi + tau2)
  est <- sum(w * yi) / sum(w)
  q       <- NA_real_
  se_conv <- sqrt(1 / sum(w))

  if (method != "fe" && hartung_knapp) {
    # The unmodified HKSJ variance, as metafor's test = "knha": it is
    # deliberately not bounded below by the conventional variance, so the
    # interval is narrower than the normal one when the weighted residual
    # mean square q is below (z / t)^2, the squared ratio of the normal
    # quantile to the t quantile.
    q    <- sum(w * (yi - est)^2) / (k - 1)
    se   <- sqrt(sum(w * (yi - est)^2) / ((k - 1) * sum(w)))
    crit <- stats::qt(1 - (1 - conf_level) / 2, df = k - 1)
    stat <- est / se
    if (se == 0 ||
        diff(range(yi)) <= 8 * .Machine$double.eps * max(abs(yi))) {
      .meta_warn_hk_zero_se(k, est, se, stat, crit, se_conv, conf_level,
                            has_prediction = k >= 3, transform = transform,
                            scale_note = scale_note)
    }
    p    <- 2 * stats::pt(-abs(stat), df = k - 1)
    test_distribution <- "t"
    test_df <- k - 1
  } else {
    se   <- se_conv
    crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    stat <- est / se
    p    <- 2 * stats::pnorm(-abs(stat))
    test_distribution <- "normal"
    test_df <- NA_real_
  }

  Q    <- .meta_Q(yi, vi)
  Q_df <- k - 1
  Q_p  <- stats::pchisq(Q, df = Q_df, lower.tail = FALSE)

  # Typical within-study variance (Higgins & Thompson, 2002), the bridge
  # between tau^2 and I^2 / H^2; mapping the tau^2 interval through it
  # gives intervals for I^2 and H^2.
  wfix <- 1 / vi
  s2_typ <- (k - 1) * sum(wfix) / (sum(wfix)^2 - sum(wfix^2))
  i2_of <- function(t2) 100 * t2 / (t2 + s2_typ)
  tau2_ci <- if (method == "fe") c(NA_real_, NA_real_) else
    .meta_tau2_ci(yi, vi, conf_level)

  pred <- c(NA_real_, NA_real_)
  if (method != "fe" && k >= 3) {
    t_pred <- stats::qt(1 - (1 - conf_level) / 2, df = k - 2)
    half   <- t_pred * sqrt(tau2 + se^2)
    pred   <- c(est - half, est + half)
  }

  list(k = k, method = method, hartung_knapp = (method != "fe" && hartung_knapp),
       estimate = est, se = se, stat = stat,
       test_distribution = test_distribution, test_df = test_df, p = p,
       lower = est - crit * se, upper = est + crit * se,
       tau2 = tau2, tau2_lower = tau2_ci[1], tau2_upper = tau2_ci[2],
       I2 = i2_of(tau2), I2_lower = i2_of(tau2_ci[1]),
       I2_upper = i2_of(tau2_ci[2]),
       H2 = (tau2 + s2_typ) / s2_typ,
       Q = Q, Q_df = Q_df, Q_p = Q_p,
       prediction_lower = pred[1], prediction_upper = pred[2],
       hartung_knapp_q = q, se_conventional = se_conv)
}

# The warning for identical effect sizes under the Hartung-Knapp
# adjustment. It says what the adjustment computes and why the result is
# zero, shows the degenerate limits, says what the result is not evidence
# of, and states the conventional standard error with the normal interval
# it gives, which the condition also carries (n_studies, estimate, se_conventional,
# lower_conventional, upper_conventional). Numbers are written to four
# significant digits. The t statistic is infinite when the standard error
# is exactly 0 and the estimate is not; 0 / 0, and so undefined, when both
# are 0; and the estimate over a rounding error when the standard error is
# only numerically zero.
.meta_warn_hk_zero_se <- function(k, est, se, stat, crit, se_conv,
                                  conf_level, has_prediction, transform,
                                  scale_note) {
  num <- function(x) format(signif(x, 4), trim = TRUE)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  lower_conv <- transform(est - z * se_conv)
  upper_conv <- transform(est + z * se_conv)
  statistic_text <- if (is.nan(stat))
    "the t statistic is undefined, since the common effect size is 0"
  else if (is.infinite(stat))
    "the t statistic is infinite"
  else paste0("the t statistic, ", num(stat), ", is the estimate divided ",
              "by a rounding error")
  .dmar_warn(c(
    "The ", k, " effect sizes are identical, exactly or up to rounding ",
    "error, so the Hartung-Knapp standard error is zero or numerically ",
    "zero. That adjustment multiplies the conventional standard error by ",
    "the square root of q, the weighted mean square of the effect sizes ",
    "about the pooled estimate, and with no differences among the effect ",
    "sizes q is 0. The confidence interval therefore has no width (",
    num(transform(est - crit * se)), ", ", num(transform(est + crit * se)),
    ")", if (has_prediction) ", nor has the prediction interval",
    ", and ", statistic_text, ". This is not evidence that the effect is ",
    "known exactly: each effect size has a sampling variance, and ",
    "agreement among ", k, " estimates, which coarse rounding or chance ",
    "can produce, does not remove it. The conventional standard error, ",
    "from the sampling variances alone, is ", num(se_conv), scale_note,
    ", with the ", format(100 * conf_level, trim = TRUE),
    "% normal interval ", num(lower_conv), " to ", num(upper_conv),
    ". The numbers are reported as computed. To report the conventional ",
    "interval, set hartung_knapp = FALSE; and check whether the effect ",
    "sizes were rounded or entered more than once."),
    "meta_hk_zero_se",
    n_studies = k, estimate = transform(est), se_conventional = se_conv,
    lower_conventional = lower_conv, upper_conventional = upper_conv)
}

# Format a .meta_fit() list as the standard table all meta_* functions
# return. The estimate scale rows can be transformed (meta_r back-transforms
# from Fisher's Z) through `transform`.
#
# The ratio of the estimate to its standard error is the row
# test_statistic under either reference distribution: Student's t on
# n_studies - 1 degrees of freedom under the Hartung-Knapp adjustment, the
# standard normal without it. DMAR 1.0.x named the row t or z by the
# reference; the NEWS rule of 1.1.0 names a statistic whose reference
# depends on an option test_statistic, and the reference is in the
# attributes test_distribution ("t" or "normal") and test_df (NA for the
# normal), which the line under the printed table states.
.meta_table <- function(fit, conf_level, transform = identity) {
  tr <- transform
  out <- data.frame(
    term = c("estimate", "se", "test_statistic", "p_value",
             "ci_lower", "ci_upper",
             "prediction_lower", "prediction_upper",
             "tau2", "tau2_ci_lower", "tau2_ci_upper", "tau",
             "I2", "I2_ci_lower", "I2_ci_upper", "H2",
             "Q", "Q_df", "Q_p", "n_studies"),
    value = c(tr(fit$estimate), fit$se, fit$stat, fit$p,
              tr(fit$lower), tr(fit$upper),
              tr(fit$prediction_lower), tr(fit$prediction_upper),
              fit$tau2, fit$tau2_lower, fit$tau2_upper, sqrt(fit$tau2),
              fit$I2, fit$I2_lower, fit$I2_upper, fit$H2,
              fit$Q, fit$Q_df, fit$Q_p, fit$k),
    stringsAsFactors = FALSE
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      p_terms = c("p_value", "Q_p"), subclass = "dmar_meta")
  attr(out, "method") <- fit$method
  attr(out, "hartung_knapp") <- fit$hartung_knapp
  attr(out, "test_distribution") <- fit$test_distribution
  attr(out, "test_df") <- fit$test_df
  # q and the conventional standard error, on the scale the model is fit
  # on, whenever the adjustment is applied.
  if (fit$hartung_knapp) {
    attr(out, "hartung_knapp_q") <- fit$hartung_knapp_q
    attr(out, "se_conventional") <- fit$se_conventional
  }
  out
}

# The line a meta-analysis table adds beneath itself: the reference
# distribution of its test_statistic row, which DMAR 1.0.x showed in the
# name of the row (t or z) and the table now states in its attributes
# test_distribution and test_df.
#' @exportS3Method .dmar_tbl_footer dmar_meta
#' @noRd
.dmar_tbl_footer.dmar_meta <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .meta_test_line(x), "test_distribution")
}

# The words of that line, or NULL for a table that records no reference.
.meta_test_line <- function(x) {
  dist <- attr(x, "test_distribution", exact = TRUE)
  if (identical(dist, "t")) {
    df <- attr(x, "test_df", exact = TRUE)
    sprintf(paste0("test_statistic: referred to Student's t on %s df ",
                   "(Hartung-Knapp adjustment)"),
            format(df, scientific = FALSE))
  } else if (identical(dist, "normal")) {
    "test_statistic: referred to the standard normal"
  } else {
    NULL
  }
}

# Shared input validation for (yi, vi) pairs.
.meta_check_yivi <- function(yi, vi) {
  if (!is.numeric(yi) || length(yi) < 2L || anyNA(yi)) {
    .dmar_refuse_argument(
      "Supply two or more effect sizes with no missing values.", arg = NULL)
  }
  if (!all(is.finite(yi))) {
    .dmar_refuse_argument(c(
      "Supply finite effect sizes: an infinite effect size cannot be ",
      "pooled."), arg = NULL)
  }
  if (!is.numeric(vi) || length(vi) != length(yi) || anyNA(vi) ||
      any(vi <= 0)) {
    .dmar_refuse_argument(
      "Supply a positive sampling variance for each effect size.", arg = NULL)
  }
  if (!all(is.finite(vi))) {
    .dmar_refuse_argument(
      "Supply a finite sampling variance for each effect size.", arg = NULL)
  }
  invisible(TRUE)
}

# --- (D) correlations corrected for attenuation ------------------------------
#
# The Fisher Z effect size and its sampling variance for meta_r(). With
# a_i = rxx_i * ryy_i the product of the two reliabilities, taken as known
# constants, the corrected correlation is r_c = r / sqrt(a), and
# y_i = atanh(r_c). Since atanh(r) has sampling variance 1 / (n - 3), the
# delta method gives y_i the variance
#
#   1 / (n - 3) * g^2,  g = d atanh(r_c) / d atanh(r)
#                         = (1 - r^2) / (sqrt(a) * (1 - r_c^2)),
#
# that is, (1 - r^2)^2 / (a (1 - r_c^2)^2 (n - 3)). The correction divides
# the sampling error of r by sqrt(a) along with r itself, and the factor
# 1 / (1 - r_c^2) grows as the corrected correlation nears 1 in magnitude,
# so g exceeds 1 whenever a < 1. g is evaluated at each study's observed r,
# so v_i depends on study i alone, as the model of .meta_fit() takes it.
# With a = 1 (no correction) r_c is r, the numerator and denominator of g
# are the same number, g is exactly 1, and y_i and v_i are exactly
# atanh(r) and 1 / (n - 3), so an analysis without reliabilities is the
# analysis of uncorrected correlations to the last bit. Uncertainty in the
# reliabilities themselves is not in v_i.
.meta_r_yivi <- function(r, n, rxx, ryy) {
  a   <- rxx * ryy
  r_c <- r / sqrt(a)
  g   <- (1 - r^2) / (sqrt(a) * (1 - r_c^2))
  list(yi = atanh(r_c), vi = g^2 / (n - 3), r_c = r_c)
}
