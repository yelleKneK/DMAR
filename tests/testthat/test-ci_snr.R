test_that("ci_snr() returns the expected 2-row tidy data frame", {
  result <- ci_snr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("lower_limit", "upper_limit"))
})

test_that("ci_snr() lower_limit <= upper_limit and both non-negative", {
  result <- ci_snr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_gte(ll, 0)
  expect_lt(ll, ul)
})

test_that("ci_snr() reports the lower-limit clamp once, in terms of the signal-to-noise ratio, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_snr(F_value = 1.2, df_1 = 4, df_2 = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the signal-to-noise ratio is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$value[result$term == "lower_limit"], 0)
  expect_gt(result$value[result$term == "upper_limit"], 0)
})

test_that("ci_snr() names itself and the inputs when the noncentral root search fails", {
  expect_error(
    ci_snr(F_value = Inf, df_1 = 4, df_2 = 50, N = 55),
    "In ci_snr\\(\\).*F-statistic Inf"
  )
})
