test_that("ss_power_sem() returns a tidy data.frame with necessary_N term", {
  res <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.10, df = 20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  n <- res$value[res$term == "necessary_N"]
  expect_true(n > 0 && n == round(n))
})

test_that("ss_power_sem() necessary_N grows as the gap between RMSEAs shrinks", {
  near <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.06, df = 20)$value[1]
  far  <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.15, df = 20)$value[1]
  expect_gt(near, far)
})

test_that("ss_power_sem() necessary_N grows with higher desired power", {
  low_pwr  <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.10, df = 20,
                           desired_power = 0.80)$value[1]
  high_pwr <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.10, df = 20,
                           desired_power = 0.99)$value[1]
  expect_gt(high_pwr, low_pwr)
})

test_that("ss_power_sem() RMSEA_true above RMSEA_null gives finite N for 'test of close fit'", {
  res <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.08, df = 30)
  n <- res$value[res$term == "necessary_N"]
  expect_true(is.finite(n))
  expect_gt(n, 0)
})

test_that("ss_power_sem() smaller df requires more N at fixed RMSEA differential", {
  small_df <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.08, df =  5)$value[1]
  large_df <- ss_power_sem(RMSEA_null = 0.05, RMSEA_true = 0.08, df = 50)$value[1]
  expect_gt(small_df, large_df)
})

test_that("ss_power_sem validates probabilities and reports realized power (MEDIUM-05)", {
  expect_error(ss_power_sem(RMSEA_null = 0, RMSEA_true = .05, df = 20, desired_power = 0),
               "strictly between 0 and 1")
  expect_error(ss_power_sem(RMSEA_null = 0, RMSEA_true = .05, df = 20, alpha_level = 1),
               "strictly between 0 and 1")
  r <- ss_power_sem(RMSEA_null = 0, RMSEA_true = .05, df = 20)
  expect_true("actual_power" %in% r$term)
  expect_gte(r$value[r$term == "actual_power"], 0.85)     # realized power meets the target
  expect_false(is.na(generics::tidy(r)$power))            # broom now finds the power
})
