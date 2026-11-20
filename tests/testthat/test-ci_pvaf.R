test_that("ci_pvaf() returns the expected 4-row tidy data frame", {
  result <- ci_pvaf(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term,
               c("ci_lower", "pvaf", "ci_upper", "actual_coverage"))
})

test_that("ci_pvaf() proportion-of-variance limits lie in [0, 1]", {
  result <- ci_pvaf(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  expect_gte(ll, 0)
  expect_lte(ul, 1)
  expect_lt(ll, ul)
})

test_that("ci_pvaf() reports the sample point estimate between its limits", {
  result <- ci_pvaf(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  est <- result$value[result$term == "pvaf"]
  # Closed form: df_1 * F / (df_1 * F + df_2), the same value eta squared
  # reports for these degrees of freedom.
  expect_equal(est, 5 * 5 / (5 * 5 + 100), tolerance = 1e-12)
  expect_gte(est, result$value[result$term == "ci_lower"])
  expect_lte(est, result$value[result$term == "ci_upper"])
  # The tail-error columns belong to the limits alone.
  no_tail <- result$term %in% c("pvaf", "actual_coverage")
  expect_true(all(is.na(result$prob_less[no_tail])))
  expect_true(all(is.na(result$prob_greater[no_tail])))
})

test_that("ci_pvaf() reports actual_coverage close to the requested conf_level", {
  result <- ci_pvaf(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  cov <- result$value[result$term == "actual_coverage"]
  expect_equal(cov, 0.95, tolerance = 1e-6)
})

test_that("ci_pvaf() reports the stated coverage when the lower limit is 0", {
  # F = 1.5 on 4 and 50 degrees of freedom is not significant (p = .216), so
  # the lower limit is 0. The interval from 0 to the upper limit has the
  # stated coverage: the lower limit can miss only by lying above the
  # population value, and 0 never does. The coverage row is therefore
  # 1 - alpha_lower - alpha_upper, and not 1 - p - alpha_upper = .759.
  result <- suppressWarnings(
    ci_pvaf(F_value = 1.5, df_1 = 4, df_2 = 50, N = 55))
  v <- stats::setNames(result$value, result$term)
  expect_identical(v[["ci_lower"]], 0)
  expect_equal(v[["actual_coverage"]], 0.95, tolerance = 1e-8)
  # The upper limit, from the noncentral F with base R alone: the
  # noncentrality that leaves .025 at or below the observed F.
  lambda_upper <- stats::uniroot(
    function(l) stats::pf(1.5, 4, 50, ncp = l) - 0.025,
    c(0, 100), tol = 1e-12)$root
  expect_equal(v[["ci_upper"]], lambda_upper / (lambda_upper + 55),
               tolerance = 1e-7)
  # The tail columns of the lower limit are the tail areas of the central
  # F at the observed F; they are not the coverage of the interval.
  p_central <- stats::pf(1.5, 4, 50, lower.tail = FALSE)
  expect_equal(result$prob_less[result$term == "ci_lower"], p_central,
               tolerance = 1e-12)
  expect_equal(result$prob_greater[result$term == "ci_lower"], 1 - p_central,
               tolerance = 1e-12)
  expect_equal(result$prob_greater[result$term == "ci_upper"], 0.025,
               tolerance = 1e-8)

  # Unequal tail areas: the coverage row is 1 - .01 - .04.
  unequal <- suppressWarnings(
    ci_pvaf(F_value = 1.5, df_1 = 4, df_2 = 50, N = 55, conf_level = NULL,
            alpha_lower = 0.01, alpha_upper = 0.04))
  expect_identical(unequal$value[unequal$term == "ci_lower"], 0)
  expect_equal(unequal$value[unequal$term == "actual_coverage"], 0.95,
               tolerance = 1e-8)

  # An F so small that the upper limit is undefined has no interval, and
  # the coverage row is NA with it.
  none <- suppressWarnings(
    ci_pvaf(F_value = 0.1, df_1 = 4, df_2 = 50, N = 55))
  expect_true(is.na(none$value[none$term == "ci_upper"]))
  expect_true(is.na(none$value[none$term == "actual_coverage"]))
  # The tail columns of the undefined upper limit are NA as well.
  expect_true(is.na(none$prob_less[none$term == "ci_upper"]))
  expect_true(is.na(none$prob_greater[none$term == "ci_upper"]))
})

test_that("ci_pvaf() tail columns hold the requested 0 and its complement on the side of a one-sided bound", {
  # With a tail area of 0 no limit is computed on that side, and the row of
  # that limit holds the requested area and its complement, whether the F
  # is significant (11.221) or not (1.5); they are not the tail areas of
  # the central F at the observed F.
  for (F_obs in c(11.221, 1.5)) {
    upper_bound <- ci_pvaf(F_value = F_obs, df_1 = 4, df_2 = 50, N = 55,
                           alpha_lower = 0, alpha_upper = 0.05)
    low <- upper_bound[upper_bound$term == "ci_lower", ]
    expect_identical(low$value, 0)
    expect_identical(low$prob_less, 0)
    expect_identical(low$prob_greater, 1)
    expect_equal(upper_bound$prob_greater[upper_bound$term == "ci_upper"],
                 0.05, tolerance = 1e-8)
    lower_bound <- suppressWarnings(
      ci_pvaf(F_value = F_obs, df_1 = 4, df_2 = 50, N = 55,
              alpha_lower = 0.05, alpha_upper = 0))
    up <- lower_bound[lower_bound$term == "ci_upper", ]
    expect_identical(up$value, 1)
    expect_identical(up$prob_less, 1)
    expect_identical(up$prob_greater, 0)
  }
  # The columns are the reverse of ci_nc_F()'s of the same names: prob_less
  # here is the area at or above the observed F, which ci_nc_F() reports as
  # prob_greater.
  both <- ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
  nc <- ci_nc_F(F_value = 11.221, df_1 = 4, df_2 = 50, verbose = TRUE)
  for (limit in c("ci_lower", "ci_upper")) {
    expect_equal(both$prob_less[both$term == limit],
                 nc$prob_greater[nc$term == limit], tolerance = 1e-12)
    expect_equal(both$prob_greater[both$term == limit],
                 nc$prob_less[nc$term == limit], tolerance = 1e-12)
  }
  expect_equal(both$prob_less[both$term == "ci_lower"], 0.025,
               tolerance = 1e-8)
})

test_that("ci_pvaf() reports the lower-limit clamp once, in terms of the proportion of variance, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_pvaf(F_value = 1.2, df_1 = 4, df_2 = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the proportion of variance accounted for is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$value[result$term == "ci_lower"], 0)
  expect_gt(result$value[result$term == "ci_upper"], 0)
})

test_that("ci_pvaf() names itself and the inputs when the noncentral root search fails", {
  expect_error(
    ci_pvaf(F_value = Inf, df_1 = 4, df_2 = 50, N = 55),
    "In ci_pvaf\\(\\).*F-statistic Inf"
  )
})
