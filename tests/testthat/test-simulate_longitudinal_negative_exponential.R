test_that("the deterministic negative exponential equals its formula", {
  d <- simulate_longitudinal_negative_exponential(
    n = 2, target_times = 0:6,
    fixed_parameters = c(alpha = 100, zeta = -80, gamma = 0.5),
    error_variance = 0
  )
  expect_equal(d$true_score, rep(100 - 80 * exp(-0.5 * (0:6)), 2))
  expect_equal(d$y, d$true_score)
  # The intercept is phi = alpha + zeta.
  expect_equal(d$true_score[d$target_time == 0], c(20, 20))
})

test_that("the long-format schema and sizes are correct", {
  set.seed(113)
  d <- simulate_longitudinal_negative_exponential(
    n = 7, target_times = c(0, 1, 3, 6),
    fixed_parameters = c(alpha = 10, zeta = -8, gamma = 0.4),
    random_variances = c(1, 1, 0.01), error_variance = 1
  )
  expect_named(d, c("id", "population", "occasion", "target_time", "time",
                    "true_score", "y"))
  expect_equal(nrow(d), 7 * 4)
  expect_s3_class(d$id, "factor")
  expect_equal(nlevels(d$id), 7)
  expect_equal(attr(d, "model"), "negative_exponential")
  expect_equal(dim(attr(d, "random_covariance")), c(3, 3))
})

test_that("timing jitter evaluates the true score at the actual time", {
  set.seed(113)
  d <- simulate_longitudinal_negative_exponential(
    n = 3, target_times = 0:4,
    fixed_parameters = c(alpha = 10, zeta = -8, gamma = 0.4),
    error_variance = 0, timing_sd = 0.2
  )
  expect_false(all(d$time == d$target_time))
  expect_equal(d$true_score, 10 - 8 * exp(-0.4 * d$time))
})

test_that("reliability and error_variance are mutually exclusive and solved", {
  expect_error(
    simulate_longitudinal_negative_exponential(
      n = 5, target_times = 0:3,
      fixed_parameters = c(10, -8, 0.4)),
    "exactly one")
  expect_error(
    simulate_longitudinal_negative_exponential(
      n = 5, target_times = 0:3,
      fixed_parameters = c(10, -8, 0.4), reliability = 0.8),
    "between-unit variation")
  set.seed(113)
  d <- simulate_longitudinal_negative_exponential(
    n = 5, target_times = 0:3,
    fixed_parameters = c(alpha = 10, zeta = -8, gamma = 0.4),
    random_variances = c(4, 4, 0.01), reliability = 0.8
  )
  expect_equal(mean(attr(d, "reliability_by_occasion")), 0.8,
               tolerance = 1e-6)
})

test_that("misnamed or wrong-length parameters are rejected", {
  expect_error(
    simulate_longitudinal_negative_exponential(
      n = 5, target_times = 0:3,
      fixed_parameters = c(a = 10, z = -8, g = 0.4), error_variance = 1),
    "Parameter names")
  expect_error(
    simulate_longitudinal_negative_exponential(
      n = 5, target_times = 0:3,
      fixed_parameters = c(10, -8), error_variance = 1),
    "length 3")
})

test_that("a named random_variances entry varies only that parameter", {
  set.seed(113)
  d <- simulate_longitudinal_negative_exponential(
    n = 40, target_times = c(0, 2, 6),
    fixed_parameters = c(alpha = 100, zeta = -80, gamma = 0.5),
    random_variances = c(gamma = 0.02), error_variance = 0
  )
  # With only gamma varying, every subject shares the intercept
  # alpha + zeta exactly, and trajectories fan out later.
  expect_equal(stats::sd(d$true_score[d$target_time == 0]), 0)
  expect_gt(stats::sd(d$true_score[d$target_time == 2]), 0.5)
  expect_error(
    simulate_longitudinal_negative_exponential(
      n = 5, target_times = 0:2,
      fixed_parameters = c(alpha = 100, zeta = -80, gamma = 0.5),
      random_variances = c(delta = 1), error_variance = 1),
    "Unknown parameter")
})
