test_that("power_equivalence_md() reproduces Phillips (1990) Table 1, row 5, col 5", {
  res <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              ltheta1 = -.2, ltheta2 = .2, ldiff = .05,
                              sigma = .20, n = 24, nu = 22)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "power")
  # Expected to ~7 sig figs from Phillips (1990).
  expect_equal(res$value, 0.8029678, tolerance = 1e-5)
})

test_that("power_equivalence_md() reproduces Diletti (1991) on the log scale", {
  res <- power_equivalence_md(alpha_level = .05, logscale = TRUE,
                              ltheta1 = .8, ltheta2 = 1.25, ldiff = 1.05,
                              sigma = .20, n = 18, nu = 16)
  expect_equal(res$value, 0.7922796, tolerance = 1e-5)
})

test_that("power_equivalence_md() power increases with sample size at the midpoint", {
  diffs_zero <- 0
  pw <- vapply(c(12, 24, 60), function(n_i) {
    nu_i <- 2 * n_i - 2  # parallel two-group design df
    power_equivalence_md(alpha_level = .05, logscale = FALSE,
                         ltheta1 = -.2, ltheta2 = .2, ldiff = diffs_zero,
                         sigma = .20, n = n_i, nu = nu_i)$value
  }, numeric(1))
  expect_true(all(diff(pw) > 0))
})

test_that("power_equivalence_md() power approaches alpha as |diff| approaches the boundary", {
  # At the equivalence boundary, the test is operating at its size (alpha).
  res <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              ltheta1 = -.2, ltheta2 = .2, ldiff = .2,
                              sigma = .20, n = 24, nu = 22)
  # Power at the boundary should be very small (well below 0.5; in fact ~ alpha).
  expect_lt(res$value, 0.10)
})

test_that("power_equivalence_md() validates its arguments", {
  expect_error(power_equivalence_md(alpha_level = 0, logscale = FALSE,
                                    ltheta1 = -.2, ltheta2 = .2, ldiff = 0,
                                    sigma = .2, n = 24, nu = 22),
               "alpha")
  expect_error(power_equivalence_md(alpha_level = .05, logscale = "yes",
                                    ltheta1 = -.2, ltheta2 = .2, ldiff = 0,
                                    sigma = .2, n = 24, nu = 22),
               "logical")
  expect_error(power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    ltheta1 = .2, ltheta2 = -.2, ldiff = 0,
                                    sigma = .2, n = 24, nu = 22),
               "ltheta2.*greater")
  expect_error(power_equivalence_md(alpha_level = .05, logscale = TRUE,
                                    ltheta1 = -1, ltheta2 = 1, ldiff = 0.5,
                                    sigma = .2, n = 24, nu = 22),
               "log scale")
})

test_that("power_density_equivalence_md() returns a tidy two-column data.frame", {
  res <- power_density_equivalence_md(
    power_sigma = c(0.05, 0.10, 0.20),
    alpha_level = .05, theta1 = -.2, theta2 = .2, diff = .05,
    sigma = .20, n = 24, nu = 22
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("power_sigma", "power_density"))
  expect_equal(nrow(res), 3L)
  expect_equal(res$power_sigma, c(0.05, 0.10, 0.20))
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
                                             theta1 = theta1, theta2 = theta2,
                                             diff = diff, sigma = sigma,
                                             n = n, nu = nu)$power_density,
    lower = 0, upper = upper
  )$value

  ref <- power_equivalence_md(alpha_level = alpha, logscale = FALSE,
                              ltheta1 = theta1, ltheta2 = theta2, ldiff = diff,
                              sigma = sigma, n = n, nu = nu)$value

  expect_equal(integ, ref, tolerance = 1e-6)
})

test_that("power_equivalence_md_plot() returns a ggplot with the expected attribute", {
  skip_if_not_installed("ggplot2")
  n  <- c(12, 24, 40)
  nu <- c(10, 22, 38)
  p <- power_equivalence_md_plot(
    alpha_level = .05, logscale = FALSE,
    theta1 = -.2, theta2 = .2, sigma = .20,
    n = n, nu = nu, subtitle = "test"
  )
  expect_s3_class(p, "ggplot")
  grid <- attr(p, "power_grid")
  expect_true(is.matrix(grid))
  expect_equal(nrow(grid), 201L)
  expect_equal(ncol(grid), length(n) + 1L)
  expect_equal(colnames(grid), c("diff", paste0("n=", n)))
  # First column should span [theta1, theta2].
  expect_equal(range(grid[, "diff"]), c(-.2, .2))
})

test_that("power_equivalence_md_plot() rejects mismatched n and nu lengths", {
  skip_if_not_installed("ggplot2")
  expect_error(
    power_equivalence_md_plot(.05, FALSE, -.2, .2, .20,
                              n = c(12, 24), nu = c(10)),
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
                         ltheta1 = -0.2 * s, ltheta2 = 0.2 * s, ldiff = 0.05 * s,
                         sigma = 0.20 * s, n = 24, nu = 22)$value
  }
  base <- power_at_scale(1)
  for (s in c(1e-6, 1e-3, 1e3, 1e6, 1e9)) {
    expect_equal(power_at_scale(s), base, tolerance = 1e-10)
  }
  expect_equal(base, 0.8029678, tolerance = 1e-6)
})
