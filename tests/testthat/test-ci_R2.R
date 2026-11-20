test_that("ci_R2() returns a data frame with ci_lower and ci_upper", {
  result <- ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_true("ci_lower" %in% result$term)
  expect_true("ci_upper" %in% result$term)
})

test_that("ci_R2() ci_lower < R2 < ci_upper for random predictors", {
  R2 <- 0.25
  result <- ci_R2(R2 = R2, N = 100, p = 5, conf_level = 0.95, random_predictors = TRUE)
  ll <- result[result$term == "ci_lower", "value"]
  ul <- result[result$term == "ci_upper", "value"]
  expect_lt(ll, R2)
  expect_gt(ul, R2)
})

test_that("ci_R2() ci_lower < R2 < ci_upper for fixed predictors", {
  R2 <- 0.25
  result <- ci_R2(R2 = R2, N = 100, p = 5, conf_level = 0.95, random_predictors = FALSE)
  ll <- result[result$term == "ci_lower", "value"]
  ul <- result[result$term == "ci_upper", "value"]
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
  width95 <- r95[r95$term == "ci_upper", "value"] - r95[r95$term == "ci_lower", "value"]
  width80 <- r80[r80$term == "ci_upper", "value"] - r80[r80$term == "ci_lower", "value"]
  expect_gt(width95, width80)
})

test_that("ci_R2() rejects a conf_level outside (0, 1)", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # A percentage-scale conf_level (95 instead of .95) once slipped through a
  # commented-out guard and silently returned a [0, 1] interval.
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 95),
               "'conf_level' must be a single number in (0, 1)", fixed = TRUE)
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 0),
               "'conf_level' must be a single number in (0, 1)", fixed = TRUE)
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = 1),
               "'conf_level' must be a single number in (0, 1)", fixed = TRUE)
  expect_error(ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = "0.95"),
               "'conf_level' must be a single number in (0, 1)", fixed = TRUE)
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
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
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
    x$value[x$term == "ci_upper"] - x$value[x$term == "ci_lower"]
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
  expect_equal(result$value[result$term == "ci_lower"], 0)
  expect_gt(result$value[result$term == "ci_upper"], 0)
})

test_that("ci_R2() with random predictors reports an upper limit that does not exist as NA, with a warning", {
  # R2 = .003 with N = 50 and p = 3 is below the .025 quantile of the
  # distribution of R2 at rho^2 = 0, Beta(3 / 2, 46 / 2): no population
  # value leaves an area of .025 at or below it.
  q <- qbeta(.025, 3 / 2, (50 - 3 - 1) / 2)
  expect_equal(q, 0.004630471, tolerance = 1e-6)
  expect_warning(r <- ci_R2(R2 = .003, N = 50, p = 3),
                 "the upper limit is undefined and reported as NA",
                 class = "dmar_result_qualified")
  up <- r$term == "ci_upper"
  expect_true(is.na(r$value[up]))
  expect_true(is.na(r$prob_less[up]))
  expect_true(is.na(r$prob_greater[up]))
  expect_identical(r$value[r$term == "ci_lower"], 0)
  # The same as with fixed predictors, where ci_nc_F() reports it.
  f <- suppressWarnings(ci_R2(R2 = .003, N = 50, p = 3,
                              random_predictors = FALSE))
  expect_true(is.na(f$value[f$term == "ci_upper"]))
  # A one-sided upper bound and ci_R() report it the same way.
  expect_warning(o <- ci_R2(R2 = .003, N = 50, p = 3, conf_level = NULL,
                            alpha_lower = 0, alpha_upper = .05),
                 class = "dmar_result_qualified")
  expect_true(is.na(o$value[o$term == "ci_upper"]))
  expect_warning(rr <- ci_R(R = sqrt(.003), N = 50, p = 3),
                 class = "dmar_result_qualified")
  expect_true(is.na(rr$value[rr$term == "ci_upper"]))
  # Just above the quantile the upper limit exists, with no warning.
  expect_no_warning(a <- ci_R2(R2 = q + 1e-4, N = 50, p = 3))
  expect_gt(a$value[a$term == "ci_upper"], 0)
  expect_equal(a$prob_greater[a$term == "ci_upper"], .025)
  # Inside ss_aipe_R2() the procedure takes that interval as running from 0
  # to 0, without the warning.
  expect_no_warning(pl <- DMAR:::.ci_R2_planning(R2 = .003, N = 50, p = 3))
  expect_identical(pl$value[pl$term == "ci_upper"], 0)
  expect_equal(pl$prob_greater[pl$term == "ci_upper"], .025)
})
