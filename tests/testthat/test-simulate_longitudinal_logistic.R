test_that("the deterministic logistic equals its formula and inflection", {
  tt <- seq(0, 12, by = 0.5)
  d <- simulate_longitudinal_logistic(
    n = 1, target_times = tt,
    fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9, zeta = 10),
    error_variance = 0
  )
  expect_equal(d$true_score, 80 / (1 + exp(-0.9 * (tt - 6))) + 10)
  # At t = beta the curve passes through half the total change.
  expect_equal(d$true_score[d$target_time == 6], 80 / 2 + 10)
})

test_that("groups differ in trajectory and sizes follow n", {
  set.seed(113)
  d <- simulate_longitudinal_logistic(
    n = c(4, 6), target_times = 0:8,
    fixed_parameters = list(
      c(alpha = 80, beta = 4, gamma = 1, zeta = 10),
      c(alpha = 60, beta = 5, gamma = 1, zeta = 10)),
    error_variance = 1
  )
  expect_equal(as.integer(table(d$population)) / 9L, c(4L, 6L))
  m <- tapply(d$true_score[d$target_time == 8], d$population[d$target_time == 8],
              mean)
  expect_gt(m[["1"]], m[["2"]])
})

test_that("the ar1 error structure shapes the error covariance", {
  set.seed(113)
  d <- simulate_longitudinal_logistic(
    n = 3, target_times = 0:3,
    fixed_parameters = c(alpha = 80, beta = 2, gamma = 1, zeta = 10),
    error_variance = 4, error_structure = "ar1", error_correlation = 0.5
  )
  Sig <- attr(d, "error_covariance")
  expect_equal(Sig[1, 2], 4 * 0.5)
  expect_equal(Sig[1, 4], 4 * 0.5^3)
  expect_equal(diag(Sig), rep(4, 4))
})

test_that("unit-specific schedules draw times inside the bounds", {
  set.seed(113)
  d <- simulate_longitudinal_logistic(
    n = 15, time_range = c(40, 90), occasions = c(5, 9),
    fixed_parameters = c(alpha = 80, beta = 65, gamma = 0.15, zeta = 10),
    random_variances = c(beta = 16), error_variance = 4
  )
  expect_true(all(d$time >= 40 & d$time <= 90))
  expect_identical(d$time, d$target_time)
  counts <- tapply(d$occasion, d$id, max)
  expect_true(all(counts >= 5 & counts <= 9))
  expect_gt(length(unique(counts)), 1)
  # Times are sorted within unit and differ across units.
  one <- d[d$id == levels(d$id)[1], ]
  expect_identical(one$time, sort(one$time))
  two <- d[d$id == levels(d$id)[2], ]
  expect_false(any(one$time %in% two$time))
  expect_identical(attr(d, "schedule"), "unit_specific")
  # A fixed occasions scalar gives every unit that many times.
  set.seed(113)
  d5 <- simulate_longitudinal_logistic(
    n = 4, time_range = c(0, 10), occasions = 6,
    fixed_parameters = c(alpha = 80, beta = 5, gamma = 0.9, zeta = 10),
    error_variance = 1
  )
  expect_true(all(tapply(d5$occasion, d5$id, max) == 6))
})

test_that("unit-specific schedules refuse the shared-schedule machinery", {
  fp <- c(alpha = 80, beta = 5, gamma = 0.9, zeta = 10)
  expect_error(
    simulate_longitudinal_logistic(
      n = 4, target_times = 0:5, time_range = c(0, 10), occasions = 6,
      fixed_parameters = fp, error_variance = 1),
    "exactly one of 'target_times'")
  expect_error(
    simulate_longitudinal_logistic(
      n = 4, time_range = c(0, 10), fixed_parameters = fp,
      error_variance = 1),
    "give 'occasions'")
  expect_error(
    simulate_longitudinal_logistic(
      n = 4, time_range = c(0, 10), occasions = 6,
      fixed_parameters = fp, random_variances = c(beta = 4),
      reliability = 0.8),
    "needs a shared measurement schedule")
  expect_error(
    simulate_longitudinal_logistic(
      n = 4, time_range = c(0, 10), occasions = 6,
      fixed_parameters = fp, error_variance = 1,
      error_structure = "ar1", error_correlation = 0.5),
    "must be \"independent\"")
  expect_error(
    simulate_longitudinal_logistic(
      n = 4, time_range = c(0, 10), occasions = 6,
      fixed_parameters = fp, error_variance = 1, timing_sd = 0.5),
    "does not combine with 'time_range'")
  expect_error(
    simulate_longitudinal_logistic(
      n = 4, time_range = c(10, 0), occasions = 6,
      fixed_parameters = fp, error_variance = 1),
    "lower < upper")
})
