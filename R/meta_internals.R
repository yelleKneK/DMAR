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

# Paule and Mandel (1982) estimator: the tau^2 at which the generalized
# Q statistic equals its expectation k - 1.
.meta_tau2_pm <- function(yi, vi) {
  k <- length(yi)
  f <- function(t2) .meta_Q(yi, vi, t2) - (k - 1)
  if (f(0) <= 0) return(0)
  upper <- max(vi) * 10 + stats::var(yi) * 10
  while (f(upper) > 0) upper <- upper * 4
  stats::uniroot(f, c(0, upper), tol = .Machine$double.eps^0.5)$root
}

# Restricted maximum likelihood, by direct maximization of the restricted
# log-likelihood profile in tau^2 (one-dimensional, smooth).
.meta_tau2_reml <- function(yi, vi) {
  ll <- function(t2) {
    w  <- 1 / (vi + t2)
    mu <- sum(w * yi) / sum(w)
    -0.5 * sum(log(vi + t2)) - 0.5 * log(sum(w)) -
      0.5 * sum(w * (yi - mu)^2)
  }
  upper <- max(vi) * 10 + stats::var(yi) * 10
  opt <- stats::optimize(ll, c(0, upper), maximum = TRUE,
                         tol = .Machine$double.eps^0.25)
  # The maximum can sit at the boundary; prefer 0 when it fits as well.
  if (ll(0) >= opt$objective - 1e-10) 0 else opt$maximum
}

# --- (B) Q-profile confidence interval for tau^2 -----------------------------
#
# Invert the generalized Q statistic (Viechtbauer, 2007): the limits are the
# tau^2 values at which Q_gen(tau2) equals the chi-square quantiles with
# k - 1 df. The lower limit is 0 when even tau^2 = 0 fails to push Q above
# the upper quantile.
.meta_tau2_ci <- function(yi, vi, conf_level = 0.95) {
  k <- length(yi)
  alpha <- 1 - conf_level
  q_hi <- stats::qchisq(1 - alpha / 2, df = k - 1)
  q_lo <- stats::qchisq(alpha / 2, df = k - 1)
  upper_bound <- max(vi) * 100 + stats::var(yi) * 100

  lo <- if (.meta_Q(yi, vi, 0) <= q_hi) 0 else {
    stats::uniroot(function(t2) .meta_Q(yi, vi, t2) - q_hi,
                   c(0, upper_bound), tol = .Machine$double.eps^0.5)$root
  }
  hi <- if (.meta_Q(yi, vi, 0) <= q_lo) 0 else {
    f <- function(t2) .meta_Q(yi, vi, t2) - q_lo
    ub <- upper_bound
    while (f(ub) > 0) ub <- ub * 4
    stats::uniroot(f, c(0, ub), tol = .Machine$double.eps^0.5)$root
  }
  c(lo, hi)
}

# --- (C) the assembled fit ---------------------------------------------------
#
# Returns the named list every public meta_* function formats. tau^2 method
# is one of "reml" (default), "pm", "dl", or "fe" (fixed effect: tau^2 = 0
# by assumption, no prediction interval). hartung_knapp applies the
# Hartung-Knapp-Sidik-Jonkman small-sample adjustment: the pooled SE is
# rescaled from the weighted residuals and the reference distribution is
# t with k - 1 df. The prediction interval follows Higgins, Thompson, and
# Spiegelhalter (2009), with t on k - 2 df.
.meta_fit <- function(yi, vi, method = "reml", hartung_knapp = TRUE,
                      conf_level = 0.95) {
  k <- length(yi)
  tau2 <- switch(method,
                 fe   = 0,
                 dl   = .meta_tau2_dl(yi, vi),
                 pm   = .meta_tau2_pm(yi, vi),
                 reml = .meta_tau2_reml(yi, vi))
  w   <- 1 / (vi + tau2)
  est <- sum(w * yi) / sum(w)

  if (method != "fe" && hartung_knapp) {
    # HKSJ variance; floored at the conventional Wald variance is not
    # applied (we report the unmodified HKSJ, as metafor's test="knha").
    se   <- sqrt(sum(w * (yi - est)^2) / ((k - 1) * sum(w)))
    crit <- stats::qt(1 - (1 - conf_level) / 2, df = k - 1)
    stat <- est / se
    p    <- 2 * stats::pt(-abs(stat), df = k - 1)
    stat_label <- "t"
  } else {
    se   <- sqrt(1 / sum(w))
    crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    stat <- est / se
    p    <- 2 * stats::pnorm(-abs(stat))
    stat_label <- "z"
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
       estimate = est, se = se, stat = stat, stat_label = stat_label, p = p,
       lower = est - crit * se, upper = est + crit * se,
       tau2 = tau2, tau2_lower = tau2_ci[1], tau2_upper = tau2_ci[2],
       I2 = i2_of(tau2), I2_lower = i2_of(tau2_ci[1]),
       I2_upper = i2_of(tau2_ci[2]),
       H2 = (tau2 + s2_typ) / s2_typ,
       Q = Q, Q_df = Q_df, Q_p = Q_p,
       prediction_lower = pred[1], prediction_upper = pred[2])
}

# Format a .meta_fit() list as the standard table all meta_* functions
# return. The estimate scale rows can be transformed (meta_r back-transforms
# from Fisher's Z) through `transform`.
.meta_table <- function(fit, conf_level, transform = identity) {
  tr <- transform
  out <- data.frame(
    term = c("estimate", "se", fit$stat_label, "p_value",
             "lower_limit", "upper_limit",
             "prediction_lower", "prediction_upper",
             "tau2", "tau2_lower", "tau2_upper", "tau",
             "I2", "I2_lower", "I2_upper", "H2",
             "Q", "Q_df", "Q_p", "k"),
    value = c(tr(fit$estimate), fit$se, fit$stat, fit$p,
              tr(fit$lower), tr(fit$upper),
              tr(fit$prediction_lower), tr(fit$prediction_upper),
              fit$tau2, fit$tau2_lower, fit$tau2_upper, sqrt(fit$tau2),
              fit$I2, fit$I2_lower, fit$I2_upper, fit$H2,
              fit$Q, fit$Q_df, fit$Q_p, fit$k),
    stringsAsFactors = FALSE
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      p_terms = c("p_value", "Q_p"))
  attr(out, "method") <- fit$method
  attr(out, "hartung_knapp") <- fit$hartung_knapp
  out
}

# Shared input validation for (yi, vi) pairs.
.meta_check_yivi <- function(yi, vi) {
  if (!is.numeric(yi) || length(yi) < 2L || anyNA(yi)) {
    stop("Supply two or more effect sizes with no missing values.",
         call. = FALSE)
  }
  if (!is.numeric(vi) || length(vi) != length(yi) || anyNA(vi) ||
      any(vi <= 0)) {
    stop("Supply a positive sampling variance for each effect size.",
         call. = FALSE)
  }
  invisible(TRUE)
}
