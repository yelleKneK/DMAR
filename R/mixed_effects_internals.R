# Internals shared by ss_aipe_mixed_effects() and
# ss_aipe_mixed_effects_sensitivity(), so that the planner and the function
# that checks it describe one design, one estimator, and one interval.
#
# The model is the balanced two-level random-intercept model
#   y_ki = gamma_0 + beta x_ki + u_k + e_ki,
#   u_k ~ N(0, tau2),  e_ki ~ N(0, sigma2),
# with beta the slope of the outcome on the predictor (true_reg_coef in
# ss_aipe_mixed_effects_sensitivity(), reg_coef in .me_one_interval()),
# for K clusters (Level 2 units) of m individuals (Level 1 units) each. With
# V the total variance of the outcome given the predictor (true_total_var)
# and rho the intraclass correlation (true_icc), tau2 = rho V and
# sigma2 = (1 - rho) V.
#
# The predictor has variance var_X. An individual-level (Level 1) predictor
# is x_ki = b_k + w_ki with Var(b_k) = icc_X var_X between clusters and
# Var(w_ki) = (1 - icc_X) var_X within clusters; a cluster-level (Level 2)
# predictor is x_k, constant within a cluster, with variance var_X between
# clusters.

# The predictor's level and the centering, checked together. Returns the
# normalized pair; refuses the one combination that has no slope (a
# cluster-level predictor centered at its own cluster mean is 0 for every
# individual).
.me_design <- function(predictor_level, centering, true_icc_X,
                       icc_arg = "true_icc_X") {
  if (!is.numeric(predictor_level) || length(predictor_level) != 1L ||
      is.na(predictor_level) || !predictor_level %in% c(1, 2)) {
    .dmar_refuse_argument(c(
      "'predictor_level' must be 1 or 2: 1 for a predictor measured on the ",
      "individuals (i.e., Level 1), which varies within a cluster, and 2 ",
      "for a predictor measured on the clusters (i.e., Level 2), which is ",
      "the same for every individual of a cluster."),
      arg = "predictor_level")
  }
  level <- as.integer(predictor_level)
  if (!is.numeric(true_icc_X) || length(true_icc_X) != 1L ||
      !is.finite(true_icc_X) || true_icc_X < 0 || true_icc_X >= 1) {
    .dmar_refuse_argument(sprintf(paste0(
      "'%s' must be a single number in [0, 1): the share of the ",
      "predictor's variance that lies between clusters. At 1 the ",
      "predictor would not vary within a cluster, which is a ",
      "cluster-level predictor (predictor_level = 2)."), icc_arg),
      arg = icc_arg)
  }
  if (level == 2L && centering == "cluster_mean") {
    .dmar_refuse_argument(c(
      "centering = \"cluster_mean\" cannot be used with predictor_level = 2. ",
      "A cluster-level predictor is the same for every individual of a ",
      "cluster, so centered at its cluster mean it is 0 everywhere and has ",
      "no slope. Use centering = \"grand_mean\" or \"none\"."),
      arg = "centering")
  }
  if (level == 2L && true_icc_X != 0) {
    .dmar_refuse_argument(sprintf(paste0(
      "'%s' describes an individual-level predictor and is not used with ",
      "predictor_level = 2, where all of the predictor's variance is ",
      "between clusters. Leave it at 0."), icc_arg),
      arg = icc_arg)
  }
  list(level = level, centering = centering)
}

# The smallest number of clusters for which the interval exists.
.me_min_clusters <- function(level) if (level == 2L) 3L else 2L

# Degrees of freedom of the t quantile of the interval.
#   Level 2 predictor: the regression of K cluster means on the predictor
#     has K - 2.
#   Level 1 predictor: the within-cluster regression has K (m - 1) - 1, and
#     the same value is used when the predictor is not centered at its
#     cluster mean, where the within-cluster part supplies most of the
#     information.
.me_df <- function(K, m, level) {
  if (level == 2L) K - 2 else K * (m - 1) - 1
}

# The square of the expected standard error of the estimated slope.
#   Level 1, centered at the cluster mean (the within-cluster slope): the
#     standard error is sqrt(MS_within / SS_within(x)), with MS_within on
#     K (m - 1) - 1 degrees of freedom and SS_within(x) distributed as
#     var_w times a chi square on K (m - 1), independent of it for a normal
#     predictor. Its expectation is exactly
#       sqrt(sigma2 / var_w) / sqrt(K (m - 1) - 1).
#   Level 2 (centered at the grand mean or not centered): the regression of
#     the K cluster means, whose standard error has expectation exactly
#       sqrt((tau2 + sigma2 / m) / var_X) / sqrt(K - 2).
#   Level 1, centered at the grand mean or not centered (the two give the
#     same slope; only the intercept moves): the generalized least squares
#     slope pools the within-cluster and the between-cluster information,
#       {K (m - 1) - 1} var_w / sigma2 +
#       (K - 2) (var_b + var_w / m) / (tau2 + sigma2 / m),
#     a first-order expression that reduces to the two exact cases when
#     either part is absent.
.me_slope_variance <- function(K, m, total_var, icc, var_X, icc_X, level,
                               centering) {
  tau2   <- icc * total_var
  sigma2 <- (1 - icc) * total_var
  if (level == 2L) {
    return((tau2 + sigma2 / m) / ((K - 2) * var_X))
  }
  var_b <- icc_X * var_X
  var_w <- (1 - icc_X) * var_X
  info_within <- (K * (m - 1) - 1) * var_w / sigma2
  if (centering == "cluster_mean") return(1 / info_within)
  info_between <- (K - 2) * (var_b + var_w / m) / (tau2 + sigma2 / m)
  1 / (info_within + info_between)
}

# Expected full width of the interval at K clusters.
.me_expected_width <- function(K, m, total_var, icc, var_X, icc_X, level,
                               centering, conf_level) {
  crit <- stats::qt(1 - (1 - conf_level) / 2, df = .me_df(K, m, level))
  2 * crit * sqrt(.me_slope_variance(K, m, total_var, icc, var_X, icc_X,
                                     level, centering))
}

# One simulated data set and its interval. Returns c(estimate, lower,
# upper), or NAs when no interval can be formed.
#
# Two of the three analyses have closed forms in a balanced design, and
# they are computed directly: the slope of a predictor centered at its
# cluster mean is the within-cluster regression slope, and the slope of a
# cluster-level predictor is the slope of the regression of the cluster
# means of the outcome on the predictor. In both, the restricted maximum
# likelihood fit of the random-intercept model gives the same estimate and
# standard error. The third, an individual-level predictor not centered at
# its cluster mean, has no closed form (its weights depend on the estimated
# variance components) and is fit with lme4::lmer().
.me_one_interval <- function(K, m, reg_coef, tau2, sigma2, var_X, icc_X,
                             level, centering, conf_level) {
  cluster <- rep(seq_len(K), each = m)
  u <- rep(stats::rnorm(K, 0, sqrt(tau2)), each = m)
  e <- stats::rnorm(K * m, 0, sqrt(sigma2))
  crit <- stats::qt(1 - (1 - conf_level) / 2, df = .me_df(K, m, level))
  if (level == 2L) {
    x_k <- stats::rnorm(K, 0, sqrt(var_X))
    y <- reg_coef * rep(x_k, each = m) + u + e
    y_bar <- as.numeric(rowsum(y, cluster)) / m
    xc <- x_k - mean(x_k)
    ssx <- sum(xc^2)
    if (!is.finite(ssx) || ssx <= 0) return(rep(NA_real_, 3L))
    b <- sum(xc * y_bar) / ssx
    res <- (y_bar - mean(y_bar)) - b * xc
    se <- sqrt(sum(res^2) / (K - 2) / ssx)
    return(c(b, b - crit * se, b + crit * se))
  }
  x <- rep(stats::rnorm(K, 0, sqrt(icc_X * var_X)), each = m) +
    stats::rnorm(K * m, 0, sqrt((1 - icc_X) * var_X))
  y <- reg_coef * x + u + e
  if (centering == "cluster_mean") {
    xw <- x - rep(as.numeric(rowsum(x, cluster)) / m, each = m)
    yw <- y - rep(as.numeric(rowsum(y, cluster)) / m, each = m)
    ssw <- sum(xw^2)
    if (!is.finite(ssw) || ssw <= 0) return(rep(NA_real_, 3L))
    b <- sum(xw * yw) / ssw
    se <- sqrt(sum((yw - b * xw)^2) / (K * (m - 1) - 1) / ssw)
    return(c(b, b - crit * se, b + crit * se))
  }
  # Centered at the grand mean or not centered: the same slope.
  if (centering == "grand_mean") x <- x - mean(x)
  cl <- factor(cluster)
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(y ~ x + (1 | cl), REML = TRUE))),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, 3L))
  b <- unname(lme4::fixef(fit)["x"])
  se <- sqrt(as.matrix(stats::vcov(fit))["x", "x"])
  if (!is.finite(b) || !is.finite(se)) return(rep(NA_real_, 3L))
  c(b, b - crit * se, b + crit * se)
}
