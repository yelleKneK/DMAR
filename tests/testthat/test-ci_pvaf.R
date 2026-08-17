test_that("ci_pvaf() returns the expected 4-row tidy data frame", {
  result <- ci_pvaf(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term,
               c("lower_limit", "pvaf", "upper_limit", "actual_coverage"))
})

test_that("ci_pvaf() proportion-of-variance limits lie in [0, 1]", {
  result <- ci_pvaf(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
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
  expect_gte(est, result$value[result$term == "lower_limit"])
  expect_lte(est, result$value[result$term == "upper_limit"])
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

test_that("ci_pvaf() reports the lower-limit clamp once, in terms of the proportion of variance, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_pvaf(F_value = 1.2, df_1 = 4, df_2 = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the proportion of variance accounted for is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$value[result$term == "lower_limit"], 0)
  expect_gt(result$value[result$term == "upper_limit"], 0)
})

test_that("ci_pvaf() names itself and the inputs when the noncentral root search fails", {
  expect_error(
    ci_pvaf(F_value = Inf, df_1 = 4, df_2 = 50, N = 55),
    "In ci_pvaf\\(\\).*F-statistic Inf"
  )
})
