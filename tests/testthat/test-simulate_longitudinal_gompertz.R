test_that("the deterministic Gompertz equals its formula and inflection", {
  tt <- seq(0, 8, by = 0.5)
  d <- simulate_longitudinal_gompertz(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 75, beta = 3, gamma = 0.55, zeta = 10),
    error_variance = 0
  )
  expect_equal(d$true_score, 75 * exp(-exp(-0.55 * (tt - 3))) + 10)
  # At t = beta the curve passes through alpha / e + zeta, about 36.8%
  # of the total change.
  expect_equal(d$true_score[d$target_time == 3], 75 / exp(1) + 10)
})

test_that("the Kelley (2005) dissertation curve is reproduced exactly", {
  # The dissertation's polynomial-comparison figure uses the Gompertz
  # y = 5 exp(-exp(7 - 1.75 t)), which in this parameterization is
  # alpha = 5, gamma = 1.75, beta = 7 / 1.75 = 4, zeta = 0.
  tt <- 0:8
  d <- simulate_longitudinal_gompertz(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 5, beta = 4, gamma = 1.75, zeta = 0),
    error_variance = 0
  )
  expect_equal(d$true_score, 5 * exp(-exp(7 - 1.75 * tt)))
})

test_that("random coefficients produce between-subject spread", {
  set.seed(113)
  d <- simulate_longitudinal_gompertz(
    n = 40, target_times = 0:8,
    fixed_parameters = c(alpha = 75, beta = 3, gamma = 0.55, zeta = 10),
    random_variances = c(alpha = 25, beta = 0.5, gamma = 0, zeta = 4),
    error_variance = 0
  )
  final <- d$true_score[d$target_time == 8]
  expect_gt(stats::sd(final), 3)
  expect_lt(abs(mean(final) - (75 * exp(-exp(-0.55 * 5)) + 10)), 3)
})
