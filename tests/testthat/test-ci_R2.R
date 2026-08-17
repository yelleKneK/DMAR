test_that("ci_R2() returns a data frame with lower_limit and upper_limit", {
  result <- ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_true("lower_limit" %in% result$term)
  expect_true("upper_limit" %in% result$term)
})

test_that("ci_R2() lower_limit < R2 < upper_limit for random predictors", {
  R2 <- 0.25
  result <- ci_R2(R2 = R2, N = 100, p = 5, conf_level = 0.95, random_predictors = TRUE)
  ll <- result[result$term == "lower_limit", "value"]
  ul <- result[result$term == "upper_limit", "value"]
  expect_lt(ll, R2)
  expect_gt(ul, R2)
})

test_that("ci_R2() lower_limit < R2 < upper_limit for fixed predictors", {
  R2 <- 0.25
  result <- ci_R2(R2 = R2, N = 100, p = 5, conf_level = 0.95, random_predictors = FALSE)
  ll <- result[result$term == "lower_limit", "value"]
  ul <- result[result$term == "upper_limit", "value"]
  expect_lt(ll, R2)
  expect_gt(ul, R2)
})

test_that("ci_R2() gives same result from F_value as from R2", {
  r_from_r2 <- ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 0.95)
  r_from_f  <- ci_R2(F_value = 6.266667, N = 100, p = 5, conf_level = 0.95)
  expect_equal(r_from_r2$value, r_from_f$value, tolerance = 1e-4)
})

test_that("ci_R2() wider CI with lower confidence level", {
  r95 <- ci_R2(R2 = 0.20, N = 80, p = 3, conf_level = 0.95)
  r80 <- ci_R2(R2 = 0.20, N = 80, p = 3, conf_level = 0.80)
  width95 <- r95[r95$term == "upper_limit", "value"] - r95[r95$term == "lower_limit", "value"]
  width80 <- r80[r80$term == "upper_limit", "value"] - r80[r80$term == "lower_limit", "value"]
  expect_gt(width95, width80)
})

test_that("ci_R2() rejects a conf_level outside (0, 1)", {
  # A percentage-scale conf_level (95 instead of .95) once slipped through a
  # commented-out guard and silently returned a [0, 1] interval.
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 95),
               "'conf_level' must be a single number in (0, 1).", fixed = TRUE)
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 0),
               "'conf_level' must be a single number in (0, 1).", fixed = TRUE)
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 1),
               "'conf_level' must be a single number in (0, 1).", fixed = TRUE)
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = "0.95"),
               "'conf_level' must be a single number in (0, 1).", fixed = TRUE)
})

test_that("ci_R2() reproduces the Kelley (2008) worked example, p. 549", {
  # Kelley (2008, Multivariate Behavioral Research, 43, 524-555, p. 549)
  # forms the interval for the Bunce and West (1995) study (R2 = .39, N = 77,
  # eight random regressors) as ci.R2(R2 = .39, N = 77, K = 8,
  # conf.level = .95), "which yields .140 and .503 for the lower and upper
  # 95% confidence limits"; the prose also reports the width as .363. The
  # limits are printed to three decimals, so the anchors are exact at that
  # precision.
  result <- ci_R2(R2 = 0.39, N = 77, p = 8, conf_level = 0.95,
                  random_predictors = TRUE)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_equal(round(ll, 3), 0.140)
  expect_equal(round(ul, 3), 0.503)
  expect_equal(round(ul - ll, 3), 0.363)
})

test_that("ci_R2() reproduces the Kelley (2008) Figure 1 crossing values, p. 533", {
  # Kelley (2008, p. 533) reads two equal-width pairs off Figure 1 (N = 100,
  # K = 5, 95% CIs; width defined by Equation 2, p. 528): the interval at
  # R2 = .1676 is as wide as the interval at R2 = .60, and the interval at
  # R2 = .5527 is as wide as the interval at R2 = .20. The crossing values
  # are printed to four decimals, which supports agreement to about 1e-4.
  width <- function(R2) {
    x <- ci_R2(R2 = R2, N = 100, p = 5, conf_level = 0.95,
               random_predictors = TRUE)
    x$value[x$term == "upper_limit"] - x$value[x$term == "lower_limit"]
  }
  expect_equal(width(0.1676), width(0.60), tolerance = 1e-4)
  expect_equal(width(0.5527), width(0.20), tolerance = 1e-4)
})

test_that("ci_R2() with fixed predictors reports the lower-limit clamp once, in terms of the squared multiple correlation, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_R2(R2 = 0.02, N = 50, p = 5, random_predictors = FALSE)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the population squared multiple correlation coefficient is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$value[result$term == "lower_limit"], 0)
  expect_gt(result$value[result$term == "upper_limit"], 0)
})
