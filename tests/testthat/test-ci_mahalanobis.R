# Tests for ci_mahalanobis() against the Hotelling T^2 / noncentral F
# derivation (Reiser, 2001; Anderson, 2003, sec. 5.2).

test_that("ci_mahalanobis() two-sample point estimate matches the quadratic-form definition", {
  set.seed(113)
  g1 <- matrix(stats::rnorm(60, mean = 0),   ncol = 3)
  g2 <- matrix(stats::rnorm(60, mean = 1.2), ncol = 3)
  S_p <- ((nrow(g1) - 1) * stats::cov(g1) + (nrow(g2) - 1) * stats::cov(g2)) /
    (nrow(g1) + nrow(g2) - 2)
  d   <- colMeans(g1) - colMeans(g2)
  D2_check <- as.numeric(crossprod(d, solve(S_p, d)))
  res <- ci_mahalanobis(group_1 = g1, group_2 = g2)
  expect_equal(res$D2, D2_check, tolerance = 1e-10)
  expect_equal(res$sample_type, "two-sample")
  expect_equal(res$p, 3)
  expect_equal(res$df_1, 3)
  expect_equal(res$df_2, nrow(g1) + nrow(g2) - 3 - 1)
})

test_that("ci_mahalanobis() one-sample point estimate matches the quadratic-form definition", {
  set.seed(113)
  g <- matrix(stats::rnorm(60, mean = 0.8), ncol = 3)
  null_mean <- c(0, 0, 0)
  S <- stats::cov(g)
  d <- colMeans(g) - null_mean
  D2_check <- as.numeric(crossprod(d, solve(S, d)))
  res <- ci_mahalanobis(group_1 = g, null_mean = null_mean)
  expect_equal(res$D2, D2_check, tolerance = 1e-10)
  expect_equal(res$sample_type, "one-sample")
  expect_true(is.na(res$n_2))
  expect_equal(res$df_1, 3)
  expect_equal(res$df_2, nrow(g) - 3)
})

test_that("ci_mahalanobis() F-value matches Hotelling T^2 transformation (two-sample)", {
  set.seed(113)
  g1 <- matrix(stats::rnorm(80, mean = 0),   ncol = 2)
  g2 <- matrix(stats::rnorm(80, mean = 1.0), ncol = 2)
  res <- ci_mahalanobis(group_1 = g1, group_2 = g2)
  n_1 <- nrow(g1); n_2 <- nrow(g2); p <- 2
  T2 <- (n_1 * n_2 / (n_1 + n_2)) * res$D2
  F_check <- T2 * (n_1 + n_2 - p - 1) / (p * (n_1 + n_2 - 2))
  expect_equal(res$F_value, F_check, tolerance = 1e-10)
})

test_that("ci_mahalanobis() confidence limits bracket the point estimate when F is large", {
  set.seed(113)
  g1 <- matrix(stats::rnorm(80, mean = 0),   ncol = 2)
  g2 <- matrix(stats::rnorm(80, mean = 1.5), ncol = 2)  # clearly separated
  res <- ci_mahalanobis(group_1 = g1, group_2 = g2)
  expect_lt(res$ci_lower, res$D2)
  expect_gt(res$ci_upper, res$D2)
  expect_gte(res$ci_lower, 0)
})

test_that("ci_mahalanobis() clamps the lower limit to 0 when F is below the critical value, with a warning", {
  set.seed(113)
  g1 <- matrix(stats::rnorm(60), ncol = 3)
  g2 <- matrix(stats::rnorm(60), ncol = 3)  # same distribution -> small D^2
  expect_warning(res <- ci_mahalanobis(group_1 = g1, group_2 = g2),
                 "below the alpha_lower critical value")
  expect_equal(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})

test_that("ci_mahalanobis() accepts pre-computed D2 and matches raw-data call", {
  set.seed(113)
  g1 <- matrix(stats::rnorm(60, mean = 0),   ncol = 3)
  g2 <- matrix(stats::rnorm(60, mean = 1.2), ncol = 3)
  raw <- ci_mahalanobis(group_1 = g1, group_2 = g2)
  pre <- ci_mahalanobis(D2 = raw$D2, n_1 = nrow(g1), n_2 = nrow(g2), p = 3)
  expect_equal(pre$ci_lower, raw$ci_lower, tolerance = 1e-8)
  expect_equal(pre$ci_upper, raw$ci_upper, tolerance = 1e-8)
  expect_equal(pre$F_value,     raw$F_value,     tolerance = 1e-10)
})

test_that("ci_mahalanobis() complains about malformed input", {
  expect_error(ci_mahalanobis(D2 = -0.1, n_1 = 30, p = 2),
               "non-negative")
  expect_error(ci_mahalanobis(), "supply raw data")
  expect_error(ci_mahalanobis(D2 = 5, n_1 = 3, p = 5),
               "degrees of freedom")
})

test_that("alpha_lower/alpha_upper require conf_level = NULL, as documented", {
  # Anchors the ?ci_mahalanobis tail-argument contract: the tails pass
  # straight through to ci_nc_F(), which refuses a non-NULL
  # conf_level beside them, so the page directs users to set
  # conf_level = NULL and supply both alphas (the page previously promised
  # a conf_level recomputation, and the documented call errored).
  expect_error(
    ci_mahalanobis(D2 = 1.2, n_1 = 40, n_2 = 45, p = 3,
                   alpha_lower = .025, alpha_upper = .025),
    "cannot mix")
  asym <- ci_mahalanobis(D2 = 1.2, n_1 = 40, n_2 = 45, p = 3,
                         conf_level = NULL,
                         alpha_lower = .025, alpha_upper = .025)
  sym  <- ci_mahalanobis(D2 = 1.2, n_1 = 40, n_2 = 45, p = 3)
  expect_equal(asym$ci_lower, sym$ci_lower)
  expect_equal(asym$ci_upper, sym$ci_upper)
})

test_that("ci_mahalanobis() reports the lower-limit clamp once, in terms of the squared Mahalanobis distance", {
  set.seed(113)
  g1 <- matrix(stats::rnorm(60), ncol = 3)
  g2 <- matrix(stats::rnorm(60), ncol = 3)
  msgs <- capture_warnings(res <- ci_mahalanobis(group_1 = g1, group_2 = g2))
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the squared Mahalanobis distance is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})

test_that("ci_mahalanobis() passes the root search controls on and refuses any other name in ...", {
  # The ... reaches ci_nc_F(), which takes tol and passes maxiter,
  # check.conv, and trace to uniroot(). Those are used; any other name
  # (a misspelling, uniroot()'s own lower or upper, which would move the
  # bracket ci_nc_F() sets, or the DMAR 1.0.x mu_0) is refused by class
  # before anything is computed, with the name in the condition's arg.
  base <- ci_mahalanobis(D2 = 0.608, n_1 = 145, n_2 = 156, p = 4)
  tight <- ci_mahalanobis(D2 = 0.608, n_1 = 145, n_2 = 156, p = 4,
                          tol = 1e-12, maxiter = 5000)
  expect_equal(tight$ci_lower, base$ci_lower, tolerance = 1e-6)
  expect_equal(tight$ci_upper, base$ci_upper, tolerance = 1e-6)
  # The limits solve the noncentral F equations: P(F <= F_obs; lambda) is
  # .975 at the lower limit and .025 at the upper, with lambda = D2 n_1 n_2
  # / (n_1 + n_2).
  lam <- function(D2) D2 * 145 * 156 / (145 + 156)
  expect_equal(pf(tight$F_value, 4, 296, ncp = lam(tight$ci_lower)), .975,
               tolerance = 1e-9)
  expect_equal(pf(tight$F_value, 4, 296, ncp = lam(tight$ci_upper)), .025,
               tolerance = 1e-9)

  for (nm in c("tolerance", "lower", "f.upper", "extendInt", "verbose",
               "conf.level")) {
    err <- expect_error(
      do.call(ci_mahalanobis, c(list(D2 = 0.608, n_1 = 145, n_2 = 156,
                                     p = 4), stats::setNames(list(1), nm))),
      class = "dmar_argument_invalid", info = nm)
    expect_identical(err$arg, nm, info = nm)
    expect_match(conditionMessage(err), nm, fixed = TRUE, info = nm)
  }
  err <- expect_error(
    ci_mahalanobis(D2 = 0.608, n_1 = 145, n_2 = 156, p = 4, tolerance = 1,
                   lower = 0),
    class = "dmar_argument_invalid")
  expect_identical(err$arg, c("tolerance", "lower"))
  # An argument without a name beyond the formals is refused as well.
  expect_error(ci_mahalanobis(0.608, NULL, NULL, NULL, 145, 156, 4, .95,
                              NULL, NULL, 1e-10),
               class = "dmar_argument_invalid")
  # The DMAR 1.0.x name keeps its own message.
  err <- expect_error(
    ci_mahalanobis(D2 = 0.608, n_1 = 145, p = 4, mu_0 = rep(0, 4)),
    "null_mean", class = "dmar_argument_invalid")
  expect_identical(err$arg, "mu_0")
})
