test_that("power_equivalence_md() reproduces Phillips (1990) Table 1, row 5, col 5", {
  res <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -.2, upper_bound = .2, true_effect = .05,
                              true_sd = .20, n_per_treatment = 24, df_error = 22)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "power")
  # Expected to ~7 sig figs from Phillips (1990).
  expect_equal(res$value, 0.8029678, tolerance = 1e-5)
})

test_that("power_equivalence_md() reproduces Diletti (1991) on the log scale", {
  res <- power_equivalence_md(alpha_level = .05, logscale = TRUE,
                              lower_bound = .8, upper_bound = 1.25, true_effect = 1.05,
                              true_sd = .20, n_per_treatment = 18, df_error = 16)
  expect_equal(res$value, 0.7922796, tolerance = 1e-5)
})

test_that("power_equivalence_md() power increases with sample size at the midpoint", {
  diffs_zero <- 0
  pw <- vapply(c(12, 24, 60), function(n_i) {
    nu_i <- 2 * n_i - 2  # parallel two-group design df
    power_equivalence_md(alpha_level = .05, logscale = FALSE,
                         lower_bound = -.2, upper_bound = .2, true_effect = diffs_zero,
                         true_sd = .20, n_per_treatment = n_i, df_error = nu_i)$value
  }, numeric(1))
  expect_true(all(diff(pw) > 0))
})

test_that("power_equivalence_md() power approaches alpha as |diff| approaches the boundary", {
  # At the equivalence boundary, the test is operating at its size (alpha).
  res <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -.2, upper_bound = .2, true_effect = .2,
                              true_sd = .20, n_per_treatment = 24, df_error = 22)
  # Power at the boundary should be very small (well below 0.5; in fact ~ alpha).
  expect_lt(res$value, 0.10)
})

test_that("power_equivalence_md() validates its arguments", {
  expect_error(power_equivalence_md(alpha_level = 0, logscale = FALSE,
                                    lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                    true_sd = .2, n_per_treatment = 24, df_error = 22),
               "alpha")
  expect_error(power_equivalence_md(alpha_level = .05, logscale = "yes",
                                    lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                    true_sd = .2, n_per_treatment = 24, df_error = 22),
               "logical")
  expect_error(power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = .2, upper_bound = -.2, true_effect = 0,
                                    true_sd = .2, n_per_treatment = 24, df_error = 22),
               "upper_bound.*greater")
  expect_error(power_equivalence_md(alpha_level = .05, logscale = TRUE,
                                    lower_bound = -1, upper_bound = 1, true_effect = 0.5,
                                    true_sd = .2, n_per_treatment = 24, df_error = 22),
               "log scale")
})

test_that("power_density_equivalence_md() returns a tidy two-column data.frame", {
  res <- power_density_equivalence_md(
    sd_values = c(0.05, 0.10, 0.20),
    alpha_level = .05, lower_bound = -.2, upper_bound = .2, true_effect = .05,
    true_sd = .20, n_per_treatment = 24, df_error = 22
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("sd_values", "power_density"))
  expect_equal(nrow(res), 3L)
  expect_equal(res$sd_values, c(0.05, 0.10, 0.20))
  expect_true(all(res$power_density >= 0))
})

test_that("power_density_equivalence_md() integrates back to the same power as power_equivalence_md()", {
  # Direct numerical check: integrate the public density function over the
  # support and compare against power_equivalence_md().
  alpha  <- 0.05
  theta1 <- -.2; theta2 <- .2; diff <- .05
  sigma  <- .20; n <- 24; nu <- 22

  power_t <- qt(1 - alpha, df = nu)
  upper   <- (theta2 - theta1) / (2 * power_t * sqrt(2 / n)) - 1e-6

  integ <- integrate(
    function(s) power_density_equivalence_md(s, alpha_level = alpha,
                                             lower_bound = theta1, upper_bound = theta2,
                                             true_effect = diff, true_sd = sigma,
                                             n_per_treatment = n, df_error = nu)$power_density,
    lower = 0, upper = upper
  )$value

  ref <- power_equivalence_md(alpha_level = alpha, logscale = FALSE,
                              lower_bound = theta1, upper_bound = theta2, true_effect = diff,
                              true_sd = sigma, n_per_treatment = n, df_error = nu)$value

  expect_equal(integ, ref, tolerance = 1e-6)
})

test_that("power_equivalence_md_plot() returns a ggplot with the expected attribute", {
  skip_if_not_installed("ggplot2")
  n  <- c(12, 24, 40)
  nu <- c(10, 22, 38)
  p <- power_equivalence_md_plot(
    alpha_level = .05, logscale = FALSE,
    lower_bound = -.2, upper_bound = .2, true_sd = .20,
    n_per_treatment = n, df_error = nu, subtitle = "test"
  )
  expect_s3_class(p, "ggplot")
  grid <- attr(p, "power_grid")
  expect_true(is.matrix(grid))
  expect_equal(nrow(grid), 201L)
  expect_equal(ncol(grid), length(n) + 1L)
  expect_equal(colnames(grid), c("true_effect", paste0("n=", n)))
  # First column should span [lower_bound, upper_bound].
  expect_equal(range(grid[, "true_effect"]), c(-.2, .2))
})

test_that("power_equivalence_md_plot() rejects mismatched n and nu lengths", {
  skip_if_not_installed("ggplot2")
  expect_error(
    power_equivalence_md_plot(.05, FALSE, -.2, .2, .20,
                              n_per_treatment = c(12, 24), df_error = c(10)),
    "same length"
  )
})

test_that("power_equivalence_md() is invariant to the measurement scale", {
  # CRITICAL-01 regression, mean-difference path. On the raw response scale,
  # rescaling all inputs by a common factor leaves the equivalence power
  # unchanged. The legacy absolute '- 1e-6' endpoint nudge pushed the
  # integration limit negative on a small scale, producing log() of a
  # negative value; a small raw scale must now return the same finite power.
  power_at_scale <- function(s) {
    power_equivalence_md(alpha_level = .05, logscale = FALSE,
                         lower_bound = -0.2 * s, upper_bound = 0.2 * s, true_effect = 0.05 * s,
                         true_sd = 0.20 * s, n_per_treatment = 24, df_error = 22)$value
  }
  base <- power_at_scale(1)
  for (s in c(1e-6, 1e-3, 1e3, 1e6, 1e9)) {
    expect_equal(power_at_scale(s), base, tolerance = 1e-10)
  }
  expect_equal(base, 0.8029678, tolerance = 1e-6)
})

# An independent calculation of the TOST power, integrating over the
# estimated difference rather than over the error standard deviation (the
# order power_equivalence_md() uses). With D ~ N(effect, (sd k)^2), k =
# sqrt(2 / n), and S^2 / sd^2 ~ chi square on df / df, both one-sided tests
# reject when D - t S k > lower and D + t S k < upper, that is when
# S < min(D - lower, upper - D) / (t k), whose probability given D is a chi
# square probability. Everything is on the analysis scale.
tost_power_by_difference <- function(alpha, lower, upper, effect, sd, n, df) {
  k  <- sqrt(2 / n)
  tq <- stats::qt(1 - alpha, df)
  f <- function(d) {
    c_max <- pmin(d - lower, upper - d) / (tq * k)
    stats::dnorm(d, effect, sd * k) * stats::pchisq(df * (c_max / sd)^2, df)
  }
  stats::integrate(f, lower, upper, rel.tol = 1e-12, abs.tol = 0)$value
}

test_that("power_equivalence_md() takes differences on the original scale when logscale = FALSE", {
  for (case in list(c(-.2, .2, .05, .20, 24, 22), c(-5, 5, 1.5, 15.67, 61, 120),
                    c(-.3, .1, -.05, .25, 40, 78))) {
    got <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                lower_bound = case[1], upper_bound = case[2],
                                true_effect = case[3], true_sd = case[4],
                                n_per_treatment = case[5],
                                df_error = case[6])$value
    ref <- tost_power_by_difference(.05, case[1], case[2], case[3], case[4],
                                    case[5], case[6])
    expect_equal(got, ref, tolerance = 1e-8, info = paste(case, collapse = " "))
  }
})

test_that("power_equivalence_md() takes ratios on the original scale when logscale = TRUE and logs them", {
  for (case in list(c(.8, 1.25, 1.05, .20, 18, 16), c(.8, 1.25, 1, .30, 24, 22),
                    c(.9, 1.11, .95, .15, 30, 28))) {
    got <- power_equivalence_md(alpha_level = .05, logscale = TRUE,
                                lower_bound = case[1], upper_bound = case[2],
                                true_effect = case[3], true_sd = case[4],
                                n_per_treatment = case[5],
                                df_error = case[6])$value
    ref <- tost_power_by_difference(.05, log(case[1]), log(case[2]),
                                    log(case[3]), case[4], case[5], case[6])
    expect_equal(got, ref, tolerance = 1e-8, info = paste(case, collapse = " "))
    # The same design given on the analysis scale with logscale = FALSE.
    on_logs <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = log(case[1]),
                                    upper_bound = log(case[2]),
                                    true_effect = log(case[3]),
                                    true_sd = case[4],
                                    n_per_treatment = case[5],
                                    df_error = case[6])$value
    expect_equal(got, on_logs, tolerance = 1e-12)
  }
  # Diletti et al. (1991), Table 1, as in the example.
  expect_equal(power_equivalence_md(alpha_level = .05, logscale = TRUE,
                                    lower_bound = .8, upper_bound = 1.25,
                                    true_effect = 1.05, true_sd = .20,
                                    n_per_treatment = 18,
                                    df_error = 16)$value,
               0.7922796, tolerance = 1e-6)
})

test_that("power_density_equivalence_md() works on the analysis scale and transforms nothing", {
  # Integrated over the error standard deviation, the density given the
  # logarithms of the ratios is the power power_equivalence_md() gives for
  # the ratios with logscale = TRUE.
  alpha <- .05; sd <- .20; n <- 18; df <- 16
  lower <- log(.8); upper <- log(1.25); effect <- log(1.05)
  s_max <- (upper - lower) / (2 * stats::qt(1 - alpha, df) * sqrt(2 / n))
  integ <- stats::integrate(function(s) {
    power_density_equivalence_md(s, alpha_level = alpha, lower_bound = lower,
                                 upper_bound = upper, true_effect = effect,
                                 true_sd = sd, n_per_treatment = n,
                                 df_error = df)$power_density
  }, lower = 0, upper = s_max, rel.tol = 1e-10)$value
  ratio_power <- power_equivalence_md(alpha_level = alpha, logscale = TRUE,
                                      lower_bound = .8, upper_bound = 1.25,
                                      true_effect = 1.05, true_sd = sd,
                                      n_per_treatment = n,
                                      df_error = df)$value
  expect_equal(integ, ratio_power, tolerance = 1e-7)
  expect_equal(integ, tost_power_by_difference(alpha, lower, upper, effect,
                                               sd, n, df), tolerance = 1e-7)
  # Given the ratios themselves, it does not log them: the region (0.8,
  # 1.25) as differences is another design, with another density.
  on_ratios <- power_density_equivalence_md(.1, alpha_level = alpha,
                                            lower_bound = .8,
                                            upper_bound = 1.25,
                                            true_effect = 1.05, true_sd = sd,
                                            n_per_treatment = n,
                                            df_error = df)$power_density
  on_logs <- power_density_equivalence_md(.1, alpha_level = alpha,
                                          lower_bound = lower,
                                          upper_bound = upper,
                                          true_effect = effect, true_sd = sd,
                                          n_per_treatment = n,
                                          df_error = df)$power_density
  expect_false(isTRUE(all.equal(on_ratios, on_logs)))
  expect_equal(on_ratios,
               power_density_equivalence_md(.1, alpha_level = alpha,
                                            lower_bound = .8 - 1.05,
                                            upper_bound = 1.25 - 1.05,
                                            true_effect = 0, true_sd = sd,
                                            n_per_treatment = n,
                                            df_error = df)$power_density,
               tolerance = 1e-12)
})

test_that("power_equivalence_md_plot() takes the bounds on the original scale and grids the true effect on it", {
  skip_if_not_installed("ggplot2")
  p <- power_equivalence_md_plot(alpha_level = .05, logscale = TRUE,
                                 lower_bound = .8, upper_bound = 1.25,
                                 true_sd = .20, n_per_treatment = c(18, 24),
                                 df_error = c(16, 22))
  grid <- attr(p, "power_grid")
  expect_identical(colnames(grid), c("true_effect", "n=18", "n=24"))
  expect_equal(range(grid[, "true_effect"]), c(.8, 1.25))
  j <- which.min(abs(grid[, "true_effect"] - 1.05))
  expect_equal(unname(grid[j, "n=18"]),
               tost_power_by_difference(.05, log(.8), log(1.25),
                                        log(grid[j, "true_effect"]), .20, 18,
                                        16), tolerance = 1e-8)
})
