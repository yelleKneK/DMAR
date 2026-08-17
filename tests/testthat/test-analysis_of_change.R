test_that("analysis_of_change() recovers the generating Gompertz population", {
  set.seed(113)
  d <- simulate_longitudinal_gompertz(
    n = 50, target_times = 0:10,
    fixed_parameters = c(alpha = 75, beta = 3, gamma = 0.55, zeta = 10),
    random_variances = c(alpha = 25, beta = 0.4, gamma = 0.005, zeta = 4),
    error_variance = 4
  )
  fit <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                      model = "gompertz")
  expect_named(fit, c("term", "estimate", "se", "sd_units",
                      "var_units"))
  expect_equal(fit$term, c("alpha", "beta", "gamma", "zeta"))
  v <- stats::setNames(fit$estimate, fit$term)
  expect_lt(abs(v["alpha"] - 75), 4)
  expect_lt(abs(v["beta"] - 3), 0.4)
  expect_lt(abs(v["gamma"] - 0.55), 0.1)
  expect_lt(abs(v["zeta"] - 10), 3)
  # Between-unit spread is reported and reflects the generating SDs
  # (inflated by estimation noise, so bounded rather than pinned).
  s <- stats::setNames(fit$sd_units, fit$term)
  expect_gt(s["alpha"], 2)
  expect_gt(s["beta"], 0.3)
  expect_equal(fit$var_units, fit$sd_units^2)
  expect_equal(attr(fit, "n_used"), 50)
  expect_equal(dim(attr(fit, "per_unit_estimates")), c(50, 4))
})

test_that("analysis_of_change() with one trajectory is a single nls fit", {
  set.seed(113)
  d <- simulate_longitudinal_negative_exponential(
    n = 1, target_times = seq(0, 8, by = 0.5),
    fixed_parameters = c(alpha = 100, zeta = -80, gamma = 0.5),
    error_variance = 0.01
  )
  fit <- analysis_of_change(d, id = NULL, time = "time", outcome = "y",
                      model = "negative_exponential")
  v <- stats::setNames(fit$estimate, fit$term)
  expect_lt(abs(v["alpha"] - 100), 0.5)
  expect_lt(abs(v["zeta"] + 80), 0.5)
  expect_lt(abs(v["gamma"] - 0.5), 0.02)
  expect_true(all(is.na(fit$sd_units)))
  expect_true(all(fit$se > 0))
  expect_equal(attr(fit, "n_used"), 1)
})

test_that("analysis_of_change() fits the logistic and Richards families", {
  set.seed(113)
  d_l <- simulate_longitudinal_logistic(
    n = 40, target_times = 0:12,
    fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9, zeta = 10),
    random_variances = c(alpha = 16, beta = 0.5, gamma = 0.005, zeta = 4),
    error_variance = 4
  )
  f_l <- analysis_of_change(d_l, id = "id", time = "time", outcome = "y",
                      model = "logistic")
  v <- stats::setNames(f_l$estimate, f_l$term)
  expect_lt(abs(v["beta"] - 6), 0.4)
  expect_lt(abs(v["alpha"] - 80), 4)

  # Richards on logistic data recovers the inflection; delta is weakly
  # identified unit by unit, so a few fits may fail to converge and
  # the function warns as designed.
  f_r <- suppressWarnings(
    analysis_of_change(d_l, id = "id", time = "time", outcome = "y",
                 model = "richards"))
  vr <- stats::setNames(f_r$estimate, f_r$term)
  expect_equal(f_r$term, c("alpha", "beta", "gamma", "delta", "zeta"))
  expect_lt(abs(vr["beta"] - 6), 0.6)
})

test_that("analysis_of_change() validates input and reports dropped fits", {
  d <- data.frame(t = 0:5, y = c(1, 2, 3, 3.5, 3.8, 4))
  expect_error(analysis_of_change(d, id = NULL, time = "nope", outcome = "y",
                            model = "gompertz"), "not in 'data'")
  expect_error(analysis_of_change(d, id = NULL, time = "t", outcome = "y",
                            model = "logistic",
                            start = c(a = 1)), "named numeric vector")
  # Two units, one of them a two-point series that cannot support a
  # four parameter fit: it is dropped with a single counted warning.
  # Person a carries a little noise because nls() cannot fit
  # zero-residual data.
  set.seed(113)
  d2 <- rbind(
    data.frame(id = "a", t = 0:8,
               y = 10 + 70 * exp(-exp(-0.9 * (0:8 - 3))) +
                 rnorm(9, 0, 0.2)),
    data.frame(id = "b", t = 0:1, y = c(5, 6))
  )
  expect_warning(
    fit <- analysis_of_change(d2, id = "id", time = "t", outcome = "y",
                        model = "gompertz"),
    "1 of 2 unit fits"
  )
  expect_equal(attr(fit, "n_used"), 1)
})

test_that("the polynomial model recovers a simulated linear population", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 50, target_times = 0:6, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1
  )
  fit <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "polynomial", order = 1)
  expect_equal(fit$term, c("b0", "b1"))
  v <- stats::setNames(fit$estimate, fit$term)
  expect_lt(abs(v["b0"] - 10), 1)
  expect_lt(abs(v["b1"] - 1.5), 0.3)
  expect_true(all(fit$se > 0))
  expect_equal(attr(fit, "n_used"), 50)
  expect_equal(attr(fit, "method"), "two_stage")
})

test_that("the mixed method purges estimation noise from sd_units (lmer)", {
  skip_if_not_installed("lme4")
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 60, target_times = 0:5, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 4
  )
  two <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "polynomial", order = 1)
  mix <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "polynomial", order = 1,
                            method = "mixed")
  expect_equal(mix$term, two$term)
  expect_equal(attr(mix, "method"), "mixed")
  # Fixed effects agree closely across methods.
  expect_equal(mix$estimate, two$estimate, tolerance = 0.05)
  # The two-stage slope SD carries estimation noise on top of the true
  # 0.5; the mixed variance component sits at or below it, closer to
  # the generating value.
  sd_two <- two$sd_units[two$term == "b1"]
  sd_mix <- mix$sd_units[mix$term == "b1"]
  expect_lt(sd_mix, sd_two)
  expect_lt(abs(sd_mix - 0.5), abs(sd_two - 0.5) + 1e-8)
  expect_gt(attr(mix, "sigma"), 1)
  expect_equal(dim(attr(mix, "per_unit_estimates")), c(60, 2))
})

test_that("the mixed method fits a nonlinear curve through nlme", {
  skip_if_not_installed("nlme")
  skip_on_cran()
  set.seed(113)
  d <- simulate_longitudinal_gompertz(
    n = 40, target_times = 0:10,
    fixed_parameters = c(alpha = 75, beta = 3, gamma = 0.55, zeta = 10),
    random_variances = c(alpha = 25, beta = 0.4, gamma = 0, zeta = 0),
    error_variance = 4
  )
  mix <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "gompertz", method = "mixed")
  v <- stats::setNames(mix$estimate, mix$term)
  expect_lt(abs(v["alpha"] - 75), 4)
  expect_lt(abs(v["beta"] - 3), 0.4)
  expect_lt(abs(v["gamma"] - 0.55), 0.1)
  expect_lt(abs(v["zeta"] - 10), 3)
  s <- stats::setNames(mix$sd_units, mix$term)
  expect_gt(s["alpha"], 1)
  expect_equal(attr(mix, "n_used"), 40)
})

test_that("analysis_of_change() validates the new arguments", {
  d <- data.frame(id = rep(c("a", "b"), each = 4), t = rep(0:3, 2),
                  y = rnorm(8))
  expect_error(
    analysis_of_change(d, id = NULL, time = "t", outcome = "y",
                       model = "polynomial", method = "mixed"),
    "needs several units")
  expect_error(
    analysis_of_change(d, id = "id", time = "t", outcome = "y",
                       model = "polynomial", order = 1.5),
    "non-negative integer")
  expect_error(
    analysis_of_change(d, id = "id", time = "t", outcome = "y",
                       model = "polynomial", start = c(b0 = 1, b1 = 1)),
    "closed form")
})

test_that("an intercept-only polynomial (order = 0) fits under both methods", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 12, target_times = 0:4,
    fixed_coefficients = c(b0 = 5),
    random_variances = c(b0 = 1), error_variance = 0.25
  )
  two <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "polynomial", order = 0)
  expect_equal(two$term, "b0")
  expect_lt(abs(two$estimate - 5), 1)
  skip_if_not_installed("lme4")
  mix <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "polynomial", order = 0,
                            method = "mixed")
  expect_equal(mix$term, "b0")
  expect_lt(abs(mix$estimate - 5), 1)
})

test_that("the polynomial no-fit message counts occasions, not 'start'", {
  d <- data.frame(id = rep(c("a", "b"), each = 3), t = rep(0:2, 2),
                  y = rnorm(6))
  expect_error(
    analysis_of_change(d, id = "id", time = "t", outcome = "y",
                       model = "polynomial", order = 4),
    "at least 5 occasions per unit")
})
