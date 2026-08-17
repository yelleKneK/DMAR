test_that("delta = 1 reproduces the logistic exactly", {
  tt <- seq(0, 12, by = 0.25)
  d_r <- simulate_longitudinal_richards(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9,
                         delta = 1, zeta = 10),
    error_variance = 0
  )
  d_l <- simulate_longitudinal_logistic(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9, zeta = 10),
    error_variance = 0
  )
  expect_equal(d_r$true_score, d_l$true_score)
})

test_that("delta near 0 approaches the Gompertz limit", {
  tt <- seq(0, 12, by = 0.5)
  d_r <- simulate_longitudinal_richards(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9,
                         delta = 1e-6, zeta = 10),
    error_variance = 0
  )
  d_g <- simulate_longitudinal_gompertz(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9, zeta = 10),
    error_variance = 0
  )
  expect_equal(d_r$true_score, d_g$true_score, tolerance = 1e-4)
})

test_that("the inflection ordinate follows alpha (1 + delta)^(-1/delta) + zeta", {
  for (dl in c(0.25, 1, 3)) {
    d <- simulate_longitudinal_richards(
      n = 1, target_times = 6,
      fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9,
                           delta = dl, zeta = 10),
      error_variance = 0
    )
    expect_equal(d$true_score, 80 * (1 + dl)^(-1 / dl) + 10)
  }
})

test_that("nonpositive delta is refused, fixed or drawn", {
  expect_error(
    simulate_longitudinal_richards(
      n = 2, target_times = 0:3,
      fixed_parameters = c(alpha = 80, beta = 2, gamma = 1,
                           delta = 0, zeta = 0),
      error_variance = 1),
    "delta <= 0|must be positive")
  set.seed(113)
  expect_error(
    simulate_longitudinal_richards(
      n = 50, target_times = 0:3,
      fixed_parameters = c(alpha = 80, beta = 2, gamma = 1,
                           delta = 0.1, zeta = 0),
      random_variances = c(0, 0, 0, 1, 0), error_variance = 1),
    "Shrink the delta")
})

test_that("a seed reproduces the draw", {
  set.seed(113)
  d1 <- simulate_longitudinal_richards(
    n = 5, target_times = 0:5,
    fixed_parameters = c(alpha = 80, beta = 3, gamma = 1,
                         delta = 2, zeta = 5),
    random_variances = 0.05, error_variance = 1)
  set.seed(113)
  d2 <- simulate_longitudinal_richards(
    n = 5, target_times = 0:5,
    fixed_parameters = c(alpha = 80, beta = 3, gamma = 1,
                         delta = 2, zeta = 5),
    random_variances = 0.05, error_variance = 1)
  expect_identical(d1, d2)
})
