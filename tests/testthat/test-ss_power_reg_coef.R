# An oracle for the necessary N, coded apart from the package: the power of
# the t test of one coefficient with p predictors, from stats::pt() with
# noncentrality sqrt(N) * f on N - p - 1 df, and the smallest N from p + 2
# up at which it reaches the target, found by brute force.
.rc_oracle_power <- function(N, f, p, alpha = .05, directional = FALSE) {
  df <- N - p - 1
  if (directional) {
    return(stats::pt(stats::qt(1 - alpha, df), df, sqrt(N) * abs(f),
                     lower.tail = FALSE))
  }
  cv <- stats::qt(1 - alpha / 2, df)
  stats::pt(-cv, df, sqrt(N) * abs(f)) +
    stats::pt(cv, df, sqrt(N) * abs(f), lower.tail = FALSE)
}
.rc_oracle_N <- function(f, p, target, alpha = .05, directional = FALSE) {
  N <- p + 2
  while (.rc_oracle_power(N, f, p, alpha, directional) < target) N <- N + 1
  N
}
.rc_value <- function(res, term) res$value[res$term == term]

test_that("a whole p gives a whole necessary_N that matches the oracle (DR39)", {
  for (p in c(1, 2, 3L, 5)) {
    res <- ss_power_reg_coef(true_f2 = .2, p = p, desired_power = .85)
    N <- .rc_value(res, "necessary_N")
    expect_identical(N, round(N))
    expect_identical(N, .rc_oracle_N(sqrt(.2), p, .85))
    expect_equal(.rc_value(res, "actual_power"),
                 .rc_oracle_power(N, sqrt(.2), p), tolerance = 1e-12)
  }
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_power_reg_coef() returns a tidy data.frame with the documented columns", {
  res <- ss_power_reg_coef(true_R2_Y_X = 0.50, true_R2_Y_X_without_j = 0.40, p = 5)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_N", "actual_power") %in% res$term))
})

test_that("ss_power_reg_coef() actual_power at returned N reaches desired_power", {
  res <- ss_power_reg_coef(true_R2_Y_X = 0.50, true_R2_Y_X_without_j = 0.40, p = 5,
                           desired_power = 0.85)
  expect_gte(res$value[res$term == "actual_power"], 0.85 - 1e-6)
})

test_that("ss_power_reg_coef() effect size is sqrt((R2 - R2_without_j)/(1 - R2))", {
  rho2 <- 0.50; rho2_j <- 0.40
  res  <- ss_power_reg_coef(true_R2_Y_X = rho2, true_R2_Y_X_without_j = rho2_j, p = 5)
  f2 <- (rho2 - rho2_j) / (1 - rho2)
  expect_equal(res$value[res$term == "true_f"], sqrt(f2),
               tolerance = 1e-9)
})

test_that("ss_power_reg_coef() rc and reg_coef agree when the parameterizations align", {
  rc_res    <- ss_power_rc(true_R2_Y_X = 0.50, true_R2_Y_X_without_j = 0.40, p = 5)
  coef_res  <- ss_power_reg_coef(true_R2_Y_X = 0.50, true_R2_Y_X_without_j = 0.40, p = 5)
  # Both functions plan for the same noncentrality, so the n's should match.
  expect_equal(rc_res$value[rc_res$term == "necessary_N"],
               coef_res$value[coef_res$term == "necessary_N"])
})

test_that("ss_power_reg_coef() Cohen's f^2 path matches the rho2 path", {
  rho2 <- 0.50; rho2_j <- 0.40
  f2 <- (rho2 - rho2_j) / (1 - rho2)
  via_cohen <- ss_power_reg_coef(true_f2 = f2, p = 5)
  via_rho   <- ss_power_reg_coef(true_R2_Y_X = rho2, true_R2_Y_X_without_j = rho2_j, p = 5)
  expect_equal(via_cohen$value[via_cohen$term == "necessary_N"],
               via_rho$value[via_rho$term == "necessary_N"])
})

test_that("ss_power_reg_coef() rejects inconsistent rho specifications", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(
    ss_power_reg_coef(true_R2_Y_X = 0.5, true_R2_Y_X_without_j = 0.4, p = 5,
                      true_cor_YX = 1:5, true_cor_XX = diag(5)),
    "do not specify"
  )
})

test_that("ss_power_reg_coef two-sided power counts both tails (HIGH-02)", {
  # A nondirectional level-alpha test must reject with probability alpha at a
  # null effect, and match the base-R noncentral-t two-tail expression at any
  # effect. The bug summed only the upper tail, returning alpha/2 at the null.
  two_tail <- function(ncp, df, a = 0.05) {
    cv <- stats::qt(1 - a / 2, df)
    stats::pt(-cv, df, ncp) + stats::pt(cv, df, ncp, lower.tail = FALSE)
  }
  for (f2 in c(0, 0.001, 0.05, 0.2)) {
    r <- ss_power_reg_coef(true_f2 = f2, p = 2, specified_N = 20)
    got <- r$value[r$term == "actual_power"]
    expect_equal(got, two_tail(sqrt(20) * sqrt(f2), 20 - 2 - 1), tolerance = 1e-10)
  }
  # At the null the power is exactly the nominal alpha, not alpha/2.
  r0 <- ss_power_reg_coef(true_f2 = 0, p = 2, specified_N = 20)
  expect_equal(r0$value[r0$term == "actual_power"], 0.05, tolerance = 1e-10)
})

test_that("ss_power_reg_coef returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_reg_coef(true_f2 = 0.6, p = 5, specified_N = n)
    r$value[r$term == "actual_power"] }
  target <- pw(7) - 1e-4                        # minimum admissible is p + 2 = 7 (df >= 1)
  r <- ss_power_reg_coef(true_f2 = 0.6, p = 5, desired_power = target)
  expect_equal(r$value[r$term == "necessary_N"], 7)
})

test_that("ss_power_reg_coef rejects invalid fixed sample sizes (MEDIUM-03)", {
  expect_error(ss_power_reg_coef(true_f2 = 0.2, p = 5, specified_N = 5), "whole number")
  expect_error(ss_power_reg_coef(true_f2 = 0.2, p = 2, specified_N = 12.5), "whole number")
})

test_that("p must be a whole number of predictors on every route (DR39)", {
  skip_on_cran()
  # Through the first builds of 1.1.0 p = 2.5 returned necessary_N = 47.5
  # (true_f2 = .2), p = 0 and p = -1 returned 47, and p = Inf or a vector
  # stopped with R's own message. Every refusal comes before any power is
  # computed: the power function is replaced by one that fails if called.
  testthat::local_mocked_bindings(
    .power_noncentral_t = function(...) stop("the search ran"))
  bad <- list(2.5, Inf, -Inf, c(2, 3), 0, -1, NA_real_)
  for (p_bad in bad) {
    calls <- list(
      quote(ss_power_reg_coef(true_f2 = .2, p = p_bad)),
      quote(ss_power_reg_coef(true_f2 = .2, p = p_bad, specified_N = 30)),
      quote(ss_power_reg_coef(true_R2_Y_X = .5, true_R2_Y_X_without_j = .4,
                              p = p_bad)),
      quote(ss_power_reg_coef(true_R2_Y_X = .5, true_R2_j_X_without_j = .3,
                              true_reg_coef_j = .3, true_sd_X = 1, true_sd_Y = 1,
                              p = p_bad)),
      quote(ss_power_rc(true_f2 = .2, p = p_bad)))
    for (cl in calls) {
      err <- tryCatch(eval(cl), error = function(e) e)
      expect_s3_class(err, "dmar_argument_invalid")
      expect_identical(err$arg, "p")
    }
  }
  # The route through true_reg_coef_j needs p as well (it stopped with
  # "argument is of length zero").
  err <- tryCatch(
    ss_power_reg_coef(true_R2_Y_X = .5, true_R2_j_X_without_j = .3,
                      true_reg_coef_j = .3, true_sd_X = 1, true_sd_Y = 1),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "p")
})

test_that("one predictor through true_cor_XX matches the direct route (DR40)", {
  # With p = 1 the restricted model has no predictors (R2 = 0); this route
  # stopped with "'a' is 0-diml" through the first builds of 1.1.0.
  f2 <- .09 / .91
  via_matrix <- ss_power_reg_coef(true_cor_XX = matrix(1, 1, 1),
                                  true_cor_YX = .3, which_predictor = 1)
  via_f2 <- ss_power_reg_coef(true_f2 = f2, p = 1)
  via_R2 <- ss_power_reg_coef(true_R2_Y_X = .09, true_R2_Y_X_without_j = 0,
                              p = 1)
  via_alias <- ss_power_rc(true_cor_XX = matrix(1, 1, 1), true_cor_YX = .3,
                           which_predictor = 1)
  expect_identical(.rc_value(via_matrix, "necessary_N"), .rc_oracle_N(sqrt(f2), 1, .85))
  expect_equal(via_matrix$value, via_f2$value, tolerance = 1e-12)
  expect_equal(via_matrix$value, via_R2$value, tolerance = 1e-12)
  expect_equal(via_alias$value, via_matrix$value)
  at_40 <- ss_power_reg_coef(true_cor_XX = matrix(1, 1, 1), true_cor_YX = .3,
                             which_predictor = 1, specified_N = 40)
  expect_equal(.rc_value(at_40, "actual_power"),
               .rc_oracle_power(40, sqrt(f2), 1), tolerance = 1e-12)
})

test_that("two or more predictors through true_cor_XX are unchanged (DR40)", {
  skip_on_cran()
  # The direct algebra: R2 = r' R^-1 r for the full and the restricted sets.
  R <- matrix(c(1, .4, .2, .4, 1, .3, .2, .3, 1), 3)
  r <- c(.5, .3, .2)
  R2_full <- drop(t(r) %*% solve(R) %*% r)
  for (j in 1:3) {
    R2_red <- drop(t(r[-j]) %*% solve(R[-j, -j]) %*% r[-j])
    f2 <- (R2_full - R2_red) / (1 - R2_full)
    res <- ss_power_reg_coef(true_cor_XX = R, true_cor_YX = r, which_predictor = j)
    expect_identical(.rc_value(res, "necessary_N"), .rc_oracle_N(sqrt(f2), 3, .85))
    expect_equal(.rc_value(res, "true_f"), sqrt(f2), tolerance = 1e-12)
  }
  # Two predictors, where the restricted set is a single one.
  res2 <- ss_power_reg_coef(true_cor_XX = R[1:2, 1:2], true_cor_YX = r[1:2],
                            which_predictor = 2)
  R2_2 <- drop(t(r[1:2]) %*% solve(R[1:2, 1:2]) %*% r[1:2])
  expect_equal(.rc_value(res2, "true_f"), sqrt((R2_2 - .25) / (1 - R2_2)),
               tolerance = 1e-12)
})

test_that("which_predictor is a whole number from 1 to p (DR40)", {
  # 3 of 2 predictors removed none and searched to N = 10,000,000 with an
  # effect of 0; 1.5 was cut to 1 by indexing.
  for (w in list(3, 0, 1.5, NA_real_, c(1, 2))) {
    err <- tryCatch(
      ss_power_reg_coef(true_cor_XX = diag(2), true_cor_YX = c(.3, .2),
                        which_predictor = w),
      error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "which_predictor")
  }
})

test_that("a zero effect with desired_power above alpha is refused at once (DR49)", {
  skip_on_cran()
  # The search ran to N = 10,000,000 (about 40 seconds) and then stopped.
  # The count of power evaluations shows the refusal comes before any.
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  real_power <- .power_noncentral_t
  testthat::local_mocked_bindings(.power_noncentral_t = function(...) {
    state$n <- state$n + 1L
    real_power(...)
  })
  for (directional in c(FALSE, TRUE)) {
    calls <- list(
      quote(ss_power_reg_coef(true_f2 = 0, p = 2, desired_power = .85,
                              directional = directional)),
      quote(ss_power_reg_coef(true_R2_Y_X = .3, true_R2_Y_X_without_j = .3,
                              p = 4, directional = directional)),
      quote(ss_power_reg_coef(true_R2_Y_X = .3, true_R2_j_X_without_j = .2,
                              true_reg_coef_j = 0, true_sd_X = 1, true_sd_Y = 1, p = 3,
                              directional = directional)),
      # A predictor uncorrelated with the outcome and the other predictor.
      quote(ss_power_reg_coef(true_cor_XX = diag(2), true_cor_YX = c(.3, 0),
                              which_predictor = 2, directional = directional)),
      quote(ss_power_rc(true_f2 = 0, p = 2, desired_power = .06,
                        directional = directional)))
    for (cl in calls) {
      state$n <- 0L
      err <- tryCatch(eval(cl), error = function(e) e)
      expect_s3_class(err, "dmar_effect_zero")
      expect_identical(state$n, 0L)
    }
  }
})

test_that("a zero effect with desired_power at or below alpha returns p + 2 (DR49)", {
  # The null rejection probability is alpha for the two-sided test and
  # for the directional one, by the central t distribution of the
  # statistic; the oracle checks that at the returned size.
  for (directional in c(FALSE, TRUE)) {
    for (target in c(.01, .05)) {
      res <- ss_power_reg_coef(true_f2 = 0, p = 3, desired_power = target,
                               directional = directional)
      expect_identical(.rc_value(res, "necessary_N"), 5)
      expect_identical(.rc_value(res, "actual_power"), .05)
      expect_equal(.rc_oracle_power(5, 0, 3, directional = directional), .05,
                   tolerance = 1e-12)
    }
  }
  # With a sample size given, the power at a zero effect is alpha.
  for (directional in c(FALSE, TRUE)) {
    r0 <- ss_power_reg_coef(true_f2 = 0, p = 2, specified_N = 30,
                            directional = directional)
    expect_equal(.rc_value(r0, "actual_power"), .05, tolerance = 1e-10)
  }
})

test_that("an effect too small for any N up to 10,000,000 is refused at once (DR49)", {
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  real_power <- .power_noncentral_t
  testthat::local_mocked_bindings(.power_noncentral_t = function(...) {
    state$n <- state$n + 1L
    real_power(...)
  })
  err <- tryCatch(ss_power_reg_coef(true_f2 = 1e-12, p = 2, desired_power = .85),
                  error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_lte(state$n, 2L)
})

test_that("the power at the planned N for one predictor holds in simulated data (DR40)", {
  skip_on_cran()
  # Monte Carlo check from raw data: a fixed predictor with mean 0 and
  # variance 1 (divisor N), a coefficient giving f2 = .09 / .91, and the
  # two-sided t test of the slope at the N the planner returns. 2000
  # replications; the tolerance is four Monte Carlo standard errors.
  res <- ss_power_reg_coef(true_cor_XX = matrix(1, 1, 1), true_cor_YX = .3,
                           which_predictor = 1, desired_power = .85)
  N <- .rc_value(res, "necessary_N")
  power <- .rc_value(res, "actual_power")
  set.seed(4021)
  x <- stats::rnorm(N)
  x <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  X <- cbind(1, x)
  beta <- sqrt(.09 / .91)
  G <- 2000L
  crit <- stats::qt(.975, N - 2)
  XtX_inv <- solve(crossprod(X))
  reject <- vapply(seq_len(G), function(g) {
    y <- beta * x + stats::rnorm(N)
    fit <- stats::lm.fit(X, y)
    s2 <- sum(fit$residuals^2) / (N - 2)
    abs(fit$coefficients[2] / sqrt(s2 * XtX_inv[2, 2])) > crit
  }, logical(1))
  mcse <- sqrt(power * (1 - power) / G)
  expect_lt(abs(mean(reject) - power), 4 * mcse)
})

test_that("a p too large for the search is refused with a class, not R's error (review)", {
  # With p above about 10,000,000 the power at N = 10,000,000 had no
  # degrees of freedom: the first round of the DR49 repair stopped with
  # "missing value where TRUE/FALSE needed" and a NaN warning. The search
  # now ends at the larger of 10,000,000 and p + 2.
  err <- tryCatch(ss_power_reg_coef(true_f2 = 1e-6, p = 2e7),
                  error = function(e) e, warning = function(w) w)
  expect_s3_class(err, "dmar_computation_failed")
})

test_that("true_cor_XX must be a correlation matrix (review)", {
  # An off-diagonal 1.2 can be inverted but is not positive definite; it
  # gave R's "NaNs produced" warning and then the general refusal of an
  # effect size specification. A matrix that is not symmetric or has a
  # diagonal other than 1 is not a correlation matrix either.
  bad <- list(matrix(c(1, 1.2, 1.2, 1), 2), matrix(c(1, .3, .2, 1), 2),
              matrix(c(2, .3, .3, 1), 2))
  for (R in bad) {
    err <- tryCatch(
      ss_power_reg_coef(true_cor_XX = R, true_cor_YX = c(.3, .2),
                        which_predictor = 1),
      error = function(e) e, warning = function(w) w)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "true_cor_XX")
  }
  # Correlations with the outcome that cannot occur with these predictors
  # (squared multiple correlation 1.25 for two uncorrelated predictors).
  err <- tryCatch(
    ss_power_reg_coef(true_cor_XX = diag(2), true_cor_YX = c(.9, .66),
                      which_predictor = 1),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, c("true_cor_XX", "true_cor_YX"))
})

test_that("ss_power_reg_coef() power is exact past a noncentrality of 37.62", {
  # N = 5, p = 2 (2 residual df), a predictor that adds .799 to R2 = .999:
  # ncp 63.2. pt() (a normal approximation there) gave .9876481868 in DMAR
  # through 3338679; the power comes from .power_noncentral_t(), shared
  # with ss_power_pcm().
  r <- ss_power_reg_coef(true_R2_Y_X = .999, true_R2_Y_X_without_j = .2, p = 2,
                         specified_N = 5, alpha_level = .001)
  ncp <- r$value[r$term == "noncentral_t_parm"]
  crit <- qt(1 - .001 / 2, 2)
  expect_equal(r$value[r$term == "actual_power"], 0.9815741922, tolerance = 1e-9)
  expect_equal(r$value[r$term == "actual_power"],
               .oracle_nct(crit, 2, ncp, lower_tail = FALSE) + .oracle_nct(-crit, 2, ncp),
               tolerance = 1e-12)
})
