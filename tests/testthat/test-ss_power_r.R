test_that("ss_power_r() returns a tidy data.frame with necessary_N and actual_power", {
  res <- ss_power_r(rho = 0.3)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_N", "actual_power", "rho") %in% res$term))
  n <- res$value[res$term == "necessary_N"]
  expect_true(n > 0 && n == round(n))
})

test_that("ss_power_r() power at returned N is at or above desired_power", {
  for (rho in c(0.2, 0.3, 0.5)) {
    res <- ss_power_r(rho = rho, desired_power = 0.80)
    expect_gte(res$value[res$term == "actual_power"], 0.80 - 1e-6)
  }
})

test_that("ss_power_r() necessary_N shrinks as rho grows", {
  small <- ss_power_r(rho = 0.20)$value[1]
  large <- ss_power_r(rho = 0.50)$value[1]
  expect_lt(large, small)
})

test_that("ss_power_r() necessary_N grows with higher desired_power", {
  low  <- ss_power_r(rho = 0.30, desired_power = 0.80)$value[1]
  high <- ss_power_r(rho = 0.30, desired_power = 0.99)$value[1]
  expect_gt(high, low)
})

test_that("ss_power_r() rho = 0.3, alpha = .05, power = .80 matches the textbook value (~85)", {
  # Standard reference value from Cohen (1988) Table 3.4.1
  res <- ss_power_r(rho = 0.30, desired_power = 0.80, alpha_level = 0.05)
  n <- res$value[res$term == "necessary_N"]
  expect_lte(abs(n - 85), 2)  # within 2 of the Cohen (1988) tabled value of 85
})

# An independent base-R statement of the Fisher's Z normal-approximation power,
# written from the definition rather than from the function under test. Under
# the alternative, atanh(r) is approximately normal with mean atanh(rho) and
# standard deviation 1 / sqrt(N - 3), so the two-sided rejection probability is
# the sum of the two tail probabilities. Vectorized over N so the minimal
# sample size can be found by an exhaustive scan.
oracle_power_r <- function(rho, rho_0, N, alpha_level, directional) {
  z_shift <- abs(atanh(rho) - atanh(rho_0)) * sqrt(N - 3)
  if (directional) {
    stats::pnorm(z_shift - stats::qnorm(1 - alpha_level))
  } else {
    crit <- stats::qnorm(1 - alpha_level / 2)
    stats::pnorm(z_shift - crit) + stats::pnorm(-z_shift - crit)
  }
}

test_that("ss_power_r() fails fast on a degenerate effect instead of hanging", {
  # With rho = 1e-12 the closed-form starting value lands near 7.8e24, where
  # N - 1 is not a distinct double, so an unguarded walk-down step can never
  # advance and the search does not terminate at all. The time limit is here so
  # that a regression fails the suite rather than hanging it.
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  setTimeLimit(elapsed = 30, transient = TRUE)
  elapsed <- system.time(
    expect_error(ss_power_r(rho = 1e-12, desired_power = 0.80),
                 "Could not reach|did not converge")
  )[["elapsed"]]
  setTimeLimit(elapsed = Inf, transient = TRUE)
  expect_lt(elapsed, 5)

  # A merely tiny effect is unreachable for the same reason and errors too.
  expect_error(ss_power_r(rho = 1e-6, desired_power = 0.80),
               "Could not reach|did not converge")
  # A near-null contrast against a nonzero null behaves the same way.
  expect_error(ss_power_r(rho = 0.30, rho_0 = 0.30 - 1e-11, desired_power = 0.80),
               "Could not reach|did not converge")
})

test_that("ss_power_r() returns the minimal N against an independent oracle", {
  # necessary_N must meet the target while necessary_N - 1 does not, and it
  # must equal the smallest N an exhaustive base-R scan accepts.
  grid <- list(
    list(rho = 0.30, rho_0 = 0,    power = 0.80, alpha = 0.05, dir = FALSE, N = 85),
    list(rho = 0.30, rho_0 = 0,    power = 0.80, alpha = 0.05, dir = TRUE,  N = 68),
    list(rho = 0.40, rho_0 = 0.20, power = 0.80, alpha = 0.05, dir = TRUE,  N = 130),
    list(rho = -0.15, rho_0 = 0,   power = 0.99, alpha = 0.01, dir = FALSE, N = 1056),
    list(rho = 0.20, rho_0 = 0,    power = 0.85, alpha = 0.05, dir = FALSE, N = 222),
    list(rho = 0.50, rho_0 = 0,    power = 0.95, alpha = 0.05, dir = FALSE, N = 47),
    list(rho = 0.05, rho_0 = 0,    power = 0.90, alpha = 0.05, dir = FALSE, N = 4199)
  )
  for (g in grid) {
    res <- ss_power_r(rho = g$rho, rho_0 = g$rho_0, desired_power = g$power,
                      alpha_level = g$alpha, directional = g$dir)
    n <- res$value[res$term == "necessary_N"]

    candidates <- 4:(2 * g$N + 50)
    accepted <- oracle_power_r(g$rho, g$rho_0, candidates, g$alpha, g$dir) >= g$power
    expect_equal(n, candidates[which(accepted)[1]])
    expect_equal(n, g$N)

    # Minimal: N meets the target, N - 1 does not.
    expect_gte(oracle_power_r(g$rho, g$rho_0, n, g$alpha, g$dir), g$power)
    expect_lt(oracle_power_r(g$rho, g$rho_0, n - 1, g$alpha, g$dir), g$power)

    # The reported power is the oracle's power at that N, to full precision.
    expect_equal(res$value[res$term == "actual_power"],
                 oracle_power_r(g$rho, g$rho_0, n, g$alpha, g$dir),
                 tolerance = 1e-12)
  }
})

test_that("ss_power_r() realized power for a specified N matches the oracle", {
  res <- ss_power_r(rho = 0.30, N = 100)
  expect_equal(res$value[res$term == "specified_N"], 100)
  expect_equal(res$value[res$term == "actual_power"],
               oracle_power_r(0.30, 0, 100, 0.05, FALSE), tolerance = 1e-12)
  # Value carried by the oracle: pnorm(atanh(.3) * sqrt(97) - qnorm(.975)) +
  # pnorm(-atanh(.3) * sqrt(97) - qnorm(.975)).
  expect_equal(res$value[res$term == "actual_power"], 0.8618021542,
               tolerance = 1e-9)
})

test_that("ss_power_r() analytic power tracks a Monte Carlo rejection rate", {
  skip_on_cran()
  # External check on the approximation itself: draw bivariate normal samples
  # at the planned N and record how often the Fisher's Z test rejects. With
  # G = 20000 the standard error of the estimate is about 0.003. The plain
  # Fisher's Z used here omits the rho / (2 (N - 1)) mean correction, so the
  # analytic 0.8003 sits a few thousandths below the true rejection rate; the
  # tolerance covers that known bias of the approximation.
  set.seed(113)
  G <- 20000
  N <- 85
  rho <- 0.30
  L <- chol(matrix(c(1, rho, rho, 1), 2))
  crit <- stats::qnorm(0.975)
  hits <- vapply(seq_len(G), function(i) {
    X <- matrix(stats::rnorm(2 * N), ncol = 2) %*% L
    abs(atanh(stats::cor(X[, 1], X[, 2])) * sqrt(N - 3)) > crit
  }, logical(1))
  res <- ss_power_r(rho = rho, desired_power = 0.80)
  expect_equal(res$value[res$term == "actual_power"], mean(hits), tolerance = 0.02)
})
