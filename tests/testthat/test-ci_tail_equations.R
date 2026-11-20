# The limits of an interval formed by inverting a noncentral distribution
# satisfy two tail equations. For an observed statistic t and an interval
# [L, U] at level 1 - alpha_lower - alpha_upper, with the statistic's
# distribution indexed by the noncentrality parameter that each limit
# implies,
#   P(T >= t | L) = alpha_lower   and   P(T <= t | U) = alpha_upper
# (Kelley, 2007, Journal of Statistical Software, 20(8)). The tests below
# evaluate both equations at the limits the public functions return, with
# stats::pt(), stats::pchisq(), and stats::pf() only, so that they do not
# depend on the root finder that produced the limits. They are independent
# of the root finder, not of R's noncentral distribution functions, and the
# map from a limit to its noncentrality parameter is written here from the
# cited source. Each grid is fixed;
# nothing is drawn at random. A limit that sits on the boundary of the
# parameter space (a lower limit of 0 for a noncentrality parameter that
# cannot be negative) solves no equation and is left out, and each grid
# is required to keep interior limits to check.

tail_tolerance <- 1e-6

# (tail probability - nominal tail area) / nominal tail area.
relative_tail_residual <- function(tail_probability, tail_area) {
  (tail_probability - tail_area) / tail_area
}

long_value <- function(x, term) x$value[match(term, x$term)]

# A lower limit of 0 is announced in a warning; the limit is checked
# against the central distribution where it occurs, so the warning is
# muffled, and no other is.
without_boundary_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("so the lower confidence limit", conditionMessage(w),
              fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("ci_smd() limits solve the noncentral t tail equations", {
  # t = d sqrt(n_1 n_2 / (n_1 + n_2)) on n_1 + n_2 - 2 degrees of freedom,
  # and a limit for delta implies the noncentrality parameter
  # delta sqrt(n_1 n_2 / (n_1 + n_2)).
  grid <- expand.grid(
    smd = c(-1.2, 0, 0.2, 0.5, 1.5),
    n = c("5,5", "12,30", "50,50", "200,150"),
    stringsAsFactors = FALSE)
  grid$conf_level <- rep(c(0.80, 0.90, 0.95, 0.99), length.out = nrow(grid))
  expect_identical(nrow(grid), 20L)
  for (i in seq_len(nrow(grid))) {
    n <- as.numeric(strsplit(grid$n[i], ",")[[1]])
    d <- grid$smd[i]
    conf_level <- grid$conf_level[i]
    tail_area <- (1 - conf_level) / 2
    x <- ci_smd(smd = d, n_1 = n[1], n_2 = n[2], conf_level = conf_level)
    scale <- sqrt(n[1] * n[2] / (n[1] + n[2]))
    df <- n[1] + n[2] - 2
    lower <- stats::pt(d * scale, df, ncp = long_value(x, "ci_lower") * scale,
                       lower.tail = FALSE)
    upper <- stats::pt(d * scale, df, ncp = long_value(x, "ci_upper") * scale)
    case <- sprintf("smd = %g, n = (%s), conf_level = %g", d, grid$n[i],
                    conf_level)
    expect_lt(abs(relative_tail_residual(lower, tail_area)), tail_tolerance,
              label = paste("lower tail residual,", case))
    expect_lt(abs(relative_tail_residual(upper, tail_area)), tail_tolerance,
              label = paste("upper tail residual,", case))
  }
})

test_that("ci_smd() limits solve the tail equations for unequal tail areas", {
  x <- ci_smd(smd = 0.4, n_1 = 25, n_2 = 40, conf_level = NULL,
              alpha_lower = 0.01, alpha_upper = 0.04)
  scale <- sqrt(25 * 40 / 65)
  lower <- stats::pt(0.4 * scale, 63, ncp = long_value(x, "ci_lower") * scale,
                     lower.tail = FALSE)
  upper <- stats::pt(0.4 * scale, 63, ncp = long_value(x, "ci_upper") * scale)
  expect_lt(abs(relative_tail_residual(lower, 0.01)), tail_tolerance)
  expect_lt(abs(relative_tail_residual(upper, 0.04)), tail_tolerance)
})

test_that("ci_rmsea() limits solve the noncentral chi square tail equations", {
  # Under the default convention the statistic an RMSEA implies is
  # T = df + rmsea^2 df (N - 1), and a limit implies the noncentrality
  # parameter limit^2 df (N - 1).
  grid <- expand.grid(
    rmsea = c(0.03, 0.05, 0.08, 0.12),
    df = c(3, 24, 60, 150, 300),
    stringsAsFactors = FALSE)
  grid$N <- rep(c(150, 400, 1000, 2500, 250), length.out = nrow(grid))
  grid$conf_level <- rep(c(0.90, 0.95, 0.80), length.out = nrow(grid))
  expect_identical(nrow(grid), 20L)
  interior_lower <- 0L
  for (i in seq_len(nrow(grid))) {
    rmsea <- grid$rmsea[i]
    df <- grid$df[i]
    multiplier <- df * (grid$N[i] - 1)
    conf_level <- grid$conf_level[i]
    tail_area <- (1 - conf_level) / 2
    # A lower limit of 0 is announced in a message; the equation for it is
    # left out below.
    x <- suppressMessages(
      ci_rmsea(rmsea = rmsea, df = df, N = grid$N[i], conf_level = conf_level))
    statistic <- df + rmsea^2 * multiplier
    ci_lower <- long_value(x, "ci_lower")
    ci_upper <- long_value(x, "ci_upper")
    case <- sprintf("rmsea = %g, df = %g, N = %g, conf_level = %g", rmsea, df,
                    grid$N[i], conf_level)
    if (ci_lower > 0) {
      interior_lower <- interior_lower + 1L
      lower <- stats::pchisq(statistic, df, ncp = ci_lower^2 * multiplier,
                             lower.tail = FALSE)
      expect_lt(abs(relative_tail_residual(lower, tail_area)), tail_tolerance,
                label = paste("lower tail residual,", case))
    } else {
      # The boundary is right only when no noncentrality parameter solves
      # the lower equation: the central tail area is already too large.
      expect_gte(stats::pchisq(statistic, df, lower.tail = FALSE), tail_area,
                 label = paste("central upper tail area,", case))
    }
    upper <- stats::pchisq(statistic, df, ncp = ci_upper^2 * multiplier)
    expect_lt(abs(relative_tail_residual(upper, tail_area)), tail_tolerance,
              label = paste("upper tail residual,", case))
  }
  expect_gte(interior_lower, 15L)
})

test_that("ci_rmsea() limits solve the tail equations under the N convention", {
  x <- ci_rmsea(rmsea = 0.07, df = 40, N = 300, conf_level = 0.90,
                sample_size_convention = "N")
  multiplier <- 40 * 300
  statistic <- 40 + 0.07^2 * multiplier
  lower <- stats::pchisq(statistic, 40,
                         ncp = long_value(x, "ci_lower")^2 * multiplier,
                         lower.tail = FALSE)
  upper <- stats::pchisq(statistic, 40,
                         ncp = long_value(x, "ci_upper")^2 * multiplier)
  expect_lt(abs(relative_tail_residual(lower, 0.05)), tail_tolerance)
  expect_lt(abs(relative_tail_residual(upper, 0.05)), tail_tolerance)
})

test_that("ci_R2() limits for fixed predictors solve the noncentral F tail equations", {
  skip_on_cran()
  # F = (R2 / p) / ((1 - R2) / (N - p - 1)) on p and N - p - 1 degrees of
  # freedom, and a limit for the population squared multiple correlation
  # implies the noncentrality parameter N limit / (1 - limit).
  grid <- expand.grid(
    R2 = c(0.10, 0.25, 0.50, 0.75),
    p = c(1, 3, 6, 12, 20),
    stringsAsFactors = FALSE)
  grid$N <- rep(c(60, 120, 300, 45, 800), length.out = nrow(grid))
  grid$conf_level <- rep(c(0.95, 0.90, 0.99), length.out = nrow(grid))
  expect_identical(nrow(grid), 20L)
  interior_lower <- 0L
  for (i in seq_len(nrow(grid))) {
    R2 <- grid$R2[i]
    p <- grid$p[i]
    N <- grid$N[i]
    conf_level <- grid$conf_level[i]
    tail_area <- (1 - conf_level) / 2
    df_2 <- N - p - 1
    x <- without_boundary_warning(
      ci_R2(R2 = R2, N = N, p = p, conf_level = conf_level,
            random_predictors = FALSE))
    F_value <- (R2 / p) / ((1 - R2) / df_2)
    ci_lower <- long_value(x, "ci_lower")
    ci_upper <- long_value(x, "ci_upper")
    case <- sprintf("R2 = %g, p = %g, N = %g, conf_level = %g", R2, p, N,
                    conf_level)
    if (ci_lower > 0) {
      interior_lower <- interior_lower + 1L
      lower <- stats::pf(F_value, p, df_2,
                         ncp = N * ci_lower / (1 - ci_lower),
                         lower.tail = FALSE)
      expect_lt(abs(relative_tail_residual(lower, tail_area)), tail_tolerance,
                label = paste("lower tail residual,", case))
    } else {
      expect_gte(stats::pf(F_value, p, df_2, lower.tail = FALSE), tail_area,
                 label = paste("central upper tail area,", case))
    }
    upper <- stats::pf(F_value, p, df_2,
                       ncp = N * ci_upper / (1 - ci_upper))
    expect_lt(abs(relative_tail_residual(upper, tail_area)), tail_tolerance,
              label = paste("upper tail residual,", case))
  }
  expect_gte(interior_lower, 15L)
})

test_that("ci_omega_squared() limits solve the noncentral F tail equations", {
  skip_on_cran()
  # A limit for omega squared implies the noncentrality parameter
  # N limit / (1 - limit) (Steiger, 2004, Equation 16, inverted).
  grid <- expand.grid(
    F_value = c(2.5, 4, 8, 20),
    df_effect = c(1, 2, 4, 6, 9),
    stringsAsFactors = FALSE)
  grid$df_error <- rep(c(30, 57, 120, 18, 400), length.out = nrow(grid))
  grid$conf_level <- rep(c(0.95, 0.90, 0.80), length.out = nrow(grid))
  expect_identical(nrow(grid), 20L)
  interior_lower <- 0L
  for (i in seq_len(nrow(grid))) {
    F_value <- grid$F_value[i]
    df_effect <- grid$df_effect[i]
    df_error <- grid$df_error[i]
    N <- df_effect + df_error + 1
    conf_level <- grid$conf_level[i]
    tail_area <- (1 - conf_level) / 2
    x <- without_boundary_warning(
      ci_omega_squared(F_value = F_value, df_effect = df_effect,
                       df_error = df_error, N = N, estimand = "partial",
                       conf_level = conf_level))
    lo <- x$omega_squared_partial_ci_lower
    up <- x$omega_squared_partial_ci_upper
    case <- sprintf("F = %g, df = (%g, %g), conf_level = %g", F_value,
                    df_effect, df_error, conf_level)
    if (lo > 0) {
      interior_lower <- interior_lower + 1L
      lower <- stats::pf(F_value, df_effect, df_error,
                         ncp = N * lo / (1 - lo),
                         lower.tail = FALSE)
      expect_lt(abs(relative_tail_residual(lower, tail_area)), tail_tolerance,
                label = paste("lower tail residual,", case))
    } else {
      expect_gte(stats::pf(F_value, df_effect, df_error, lower.tail = FALSE),
                 tail_area, label = paste("central upper tail area,", case))
    }
    upper <- stats::pf(F_value, df_effect, df_error,
                       ncp = N * up / (1 - up))
    expect_lt(abs(relative_tail_residual(upper, tail_area)), tail_tolerance,
              label = paste("upper tail residual,", case))
  }
  expect_gte(interior_lower, 15L)

  # A design in which N is not df_effect + df_error + 1 (a factorial
  # design, where other effects take degrees of freedom): the noncentrality
  # parameter is N limit / (1 - limit) with the N given, which the grid
  # above cannot tell from the sum of the degrees of freedom plus 1.
  x <- ci_omega_squared(F_value = 6, df_effect = 2, df_error = 54, N = 60,
                        estimand = "partial", conf_level = 0.95)
  lo <- x$omega_squared_partial_ci_lower
  up <- x$omega_squared_partial_ci_upper
  expect_gt(lo, 0)
  lower <- stats::pf(6, 2, 54, ncp = 60 * lo / (1 - lo),
                     lower.tail = FALSE)
  upper <- stats::pf(6, 2, 54, ncp = 60 * up / (1 - up))
  expect_lt(abs(relative_tail_residual(lower, 0.025)), tail_tolerance)
  expect_lt(abs(relative_tail_residual(upper, 0.025)), tail_tolerance)
})
