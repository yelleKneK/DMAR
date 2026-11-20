# Shared numerics for the two-phase (double) sampling family:
# ci_two_phase(), ss_aipe_two_phase(), ss_aipe_two_phase_cost(),
# ss_aipe_two_phase_sensitivity(), and ss_seq_two_phase().
#
# Notation used throughout. A pool of M units is drawn at random from a
# population; an inexpensive proxy measurement X is observed on every unit
# of the pool, and the target measurement Y on a simple random subsample of
# n of them (the labeled units). The remaining N = M - n units carry the
# proxy alone. The estimand is mu_Y, the population mean of the target.
# The regression estimator of double sampling is
#   mean_Y = ybar_n + b (xbar_M - xbar_n),
# with b the least squares slope of Y on X in the labeled units, and to
# first order
#   Var(mean_Y) = sigma_Y^2 (1 - rho^2) / n + rho^2 sigma_Y^2 / M.
# The estimate is named mean_Y, as the row of ci_two_phase() that reports it
# is.


# Validate and clean the data of one two-phase sample. Labeled units enter
# only when the target, the proxy, and (if given) the cluster are all
# present; unlabeled units enter only when the proxy (and, if given, the
# cluster) is present. Logical vectors are taken as 0/1 measurements, the
# usual coding of a binary label. Returns a list with the cleaned vectors
# and the counts dropped, or stops with a message naming the problem.
.two_phase_prepare <- function(y, x, x_unlabeled, cluster = NULL,
                               cluster_unlabeled = NULL,
                               min_unlabeled = 1L) {
  as_measure <- function(v, name) {
    if (is.logical(v)) v <- as.numeric(v)
    if (!is.numeric(v) || !is.null(dim(v)))
      .dmar_refuse_argument(c(
        "'", name, "' must be a numeric (or logical) vector."), arg = name)
    if (any(is.infinite(v)))
      .dmar_refuse_argument(c(
        "'", name, "' contains infinite values."), arg = name)
    as.numeric(v)
  }
  y <- as_measure(y, "y")
  x <- as_measure(x, "x")
  x_unlabeled <- as_measure(x_unlabeled, "x_unlabeled")
  if (length(y) != length(x))
    .dmar_refuse_argument(c(
      "'y' and 'x' must have the same length: one target and one proxy ",
      "value per labeled unit."), arg = c("y", "x"))
  if (is.null(cluster) != is.null(cluster_unlabeled))
    .dmar_refuse_argument(
      "Supply both 'cluster' and 'cluster_unlabeled', or neither.",
      arg = c("cluster", "cluster_unlabeled"))
  clustered <- !is.null(cluster)
  if (clustered) {
    if (length(cluster) != length(y))
      .dmar_refuse_argument(c(
        "'cluster' must have one entry per labeled unit (the length of ",
        "'y')."), arg = c("cluster", "y"))
    if (length(cluster_unlabeled) != length(x_unlabeled))
      .dmar_refuse_argument(c(
        "'cluster_unlabeled' must have one entry per unlabeled unit (the ",
        "length of 'x_unlabeled')."),
        arg = c("cluster_unlabeled", "x_unlabeled"))
  }

  keep_lab <- !is.na(y) & !is.na(x)
  keep_unl <- !is.na(x_unlabeled)
  if (clustered) {
    keep_lab <- keep_lab & !is.na(cluster)
    keep_unl <- keep_unl & !is.na(cluster_unlabeled)
  }
  out <- list(
    y = y[keep_lab], x = x[keep_lab], x_unlabeled = x_unlabeled[keep_unl],
    cluster = if (clustered) as.character(cluster[keep_lab]),
    cluster_unlabeled = if (clustered) as.character(cluster_unlabeled[keep_unl]),
    dropped_labeled = sum(!keep_lab), dropped_unlabeled = sum(!keep_unl)
  )
  if (length(out$y) < 3L)
    .dmar_refuse_argument(c(
      "At least 3 labeled units with both a target and a proxy value are ",
      "required; ", length(out$y), " remain after units with a missing ",
      "value are set aside."), arg = NULL)
  if (length(out$x_unlabeled) < min_unlabeled)
    .dmar_refuse_argument(c(
      "At least ", min_unlabeled, " unlabeled unit",
      if (min_unlabeled > 1L) "s", " with a proxy value ",
      if (min_unlabeled > 1L) "are" else "is", " required."), arg = NULL)
  out
}


# The regression estimator of double sampling from the labeled pairs and
# three summaries of the unlabeled proxies: their count n_u, mean mean_u,
# and sum of squared deviations ss_u. Taking summaries rather than the
# vector is what lets the sensitivity loop draw the unlabeled part of a
# normal pool exactly from its sampling distribution (the mean is normal
# and the sum of squares a scaled chi square, independently) instead of
# simulating every unit of a pool of tens of thousands. Returns a named
# numeric vector: the estimate, the least squares slope, the residual mean
# square of the target on the proxy (n - 2 degrees of freedom), the proxy's
# variance over the whole pool, and the two sizes. The slope is NA when the
# proxy is constant on the labeled units, where no regression can be fitted.
.two_phase_regression <- function(y, x, n_u, mean_u, ss_u) {
  # Both sizes are doubles. length() returns an integer, and so does n_u
  # when it is a length or a difference of integer sizes; n * n_u below
  # would then be formed as an integer, which overflows to NA once it
  # passes 2^31 - 1 (2,000 labeled units with 1.1 million unlabeled do).
  n      <- as.numeric(length(y))
  n_u    <- as.numeric(n_u)
  M      <- n + n_u
  ybar   <- sum(y) / n
  xbar_n <- sum(x) / n
  dx     <- x - xbar_n
  sxx    <- sum(dx * dx)
  if (!(sxx > 0)) {
    return(c(mean_Y = NA_real_, slope = NA_real_, var_residual = NA_real_,
             var_pool = NA_real_, n_labeled = n, n_pool = M))
  }
  b      <- sum(dx * (y - ybar)) / sxx
  e      <- y - ybar - b * dx
  xbar_M <- (n * xbar_n + n_u * mean_u) / M
  # Pooled sum of squares of the proxy over all M units: the two within
  # sums plus the between term for the difference of the two means.
  ss_M   <- sxx + ss_u + n * n_u / M * (xbar_n - mean_u)^2
  c(mean_Y       = ybar + b * (xbar_M - xbar_n),
    slope        = b,
    var_residual = sum(e * e) / (n - 2),
    var_pool     = ss_M / (M - 1),
    n_labeled    = n,
    n_pool       = M)
}


# Fast path for the unclustered regression interval, the one
# ss_aipe_two_phase() plans for: a bare named numeric vector rather than
# the dmar_tbl that ci_two_phase() returns, which skips the data.frame
# construction in the Monte Carlo loop of ss_aipe_two_phase_sensitivity().
# ci_two_phase() calls it too, so the interval the planner, the sensitivity
# sibling, and the sequential rule work with is the one the user is shown.
# The standard error estimates sigma_e^2 / n + b^2 sigma_X^2 / M, with the
# residual mean square for the first term and the proxy's variance over the
# whole pool for the second, and the reference distribution is t on the
# n - 2 residual degrees of freedom of the first term.
.ci_two_phase_fast <- function(y, x, n_u, mean_u, ss_u, conf_level) {
  r  <- .two_phase_regression(y, x, n_u, mean_u, ss_u)
  n  <- r[["n_labeled"]]
  se <- sqrt(r[["var_residual"]] / n + r[["slope"]]^2 * r[["var_pool"]] / r[["n_pool"]])
  df <- n - 2
  q  <- stats::qt(1 - (1 - conf_level) / 2, df)
  c(mean_Y      = r[["mean_Y"]],
    se          = se,
    df          = df,
    ci_lower = r[["mean_Y"]] - q * se,
    ci_upper = r[["mean_Y"]] + q * se,
    slope       = r[["slope"]])
}


# The search behind ss_aipe_two_phase(): the smallest n >= n_min from which
# a criterion holds at every larger n. meets(n) evaluates the criterion at
# one n; holds_from(n) is TRUE only when the criterion holds at n and at
# every larger n, and holds_on(a, b) only when it holds at every n from a to
# b. Either certificate may decline when the criterion does hold, which
# costs only further evaluation, but neither may accept when it fails
# somewhere. The search doubles from `start` to a point certified by
# holds_from(), then walks down from it in blocks, doubling the block after
# a certified one and halving it after one that is not, until the criterion
# fails at the unit just below the certified stretch; the bottom of that
# stretch is the answer. `start` need not bound the answer from either side,
# so the search is valid for any assurance, including one below one half,
# and for a criterion that holds at a few small n, fails over a stretch of
# larger ones, and then holds for good, as the assurance criterion can with
# a large kurtosis.
.two_phase_search <- function(meets, holds_from, holds_on, start, n_min) {
  u <- max(n_min, ceiling(start))
  while (!holds_from(u)) {
    u <- 2 * u
    if (u > 2^52)
      .dmar_refuse_argument(c(
        "The sample size search found no n below 4.5e15 from which the ",
        "criterion is certain to hold; the width is too narrow for these ",
        "planning values."), arg = NULL)
  }
  b <- u
  k <- 1
  while (b > n_min) {
    k <- min(k, b - n_min)
    if (k == 1) {
      if (!meets(b - 1)) return(b)
      b <- b - 1
      k <- 2
    } else if (holds_on(b - k, b - 1)) {
      b <- b - k
      k <- 2 * k
    } else {
      k <- k %/% 2
    }
  }
  b
}


# Checks shared by the planners and the sensitivity sibling. Each argument
# is checked only when supplied, so a caller passes the ones it has.
.two_phase_check_planning <- function(sd_y = NULL, rho = NULL, width = NULL,
                                      conf_level = NULL, assurance = NULL,
                                      true_kurtosis = NULL,
                                      sd_name = "true_sd_Y",
                                      rho_name = "true_cor") {
  single <- function(v) is.numeric(v) && length(v) == 1L && !is.na(v)
  if (!is.null(sd_y) && !(single(sd_y) && is.finite(sd_y) && sd_y > 0))
    .dmar_refuse_argument(c(
      "'", sd_name, "' must be a single finite positive number."), arg = sd_name)
  if (!is.null(rho) && !(single(rho) && abs(rho) < 1))
    .dmar_refuse_argument(c(
      "'", rho_name, "' must be a single number in (-1, 1)."), arg = rho_name)
  if (!is.null(width) && !(single(width) && is.finite(width) && width > 0))
    .dmar_refuse_argument(c(
      "'width' must be a single finite positive number (the full width ",
      "of the interval)."), arg = "width")
  if (!is.null(conf_level) && !(single(conf_level) && conf_level > 0 &&
                                conf_level < 1))
    .dmar_refuse_argument(
      "'conf_level' must be a single number strictly between 0 and 1.",
      arg = "conf_level")
  if (!is.null(assurance) && !(single(assurance) && assurance > 0 &&
                               assurance < 1))
    .dmar_refuse_argument(c(
      "'assurance' must be NULL or a single number strictly between 0 ",
      "and 1."), arg = "assurance")
  if (!is.null(true_kurtosis) && !(single(true_kurtosis) && is.finite(true_kurtosis) &&
                              true_kurtosis > -2))
    .dmar_refuse_argument(c(
      "'true_kurtosis' must be a single number greater than -2 (excess ",
      "kurtosis, the scale of kurtosis())."), arg = "true_kurtosis")
  invisible(TRUE)
}
