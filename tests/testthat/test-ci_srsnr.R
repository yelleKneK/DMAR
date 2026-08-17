test_that("ci_srsnr() returns the expected 2-row tidy data frame", {
  result <- ci_srsnr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("lower_limit", "upper_limit"))
})

test_that("ci_srsnr() limits are the square roots of the corresponding ci_snr() limits", {
  snr_result    <- ci_snr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  srsnr_result  <- ci_srsnr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_equal(srsnr_result$value, sqrt(pmax(0, snr_result$value)), tolerance = 1e-8)
})

test_that("ci_srsnr() design-stage call matches the implied F-value path", {
  means <- c(94, 91, 92, 83)
  sig2  <- 67.375
  n     <- 6
  J     <- length(means)
  N     <- J * n
  df1   <- J - 1
  df2   <- N - J
  mu_bar      <- mean(means)
  SS_between  <- n * sum((means - mu_bar)^2)
  F_implied   <- (SS_between / df1) / sig2

  # The small implied F legitimately clamps the lower limit to 0 on both paths.
  expect_warning(
    r_design <- ci_srsnr(means = means, sigma_squared = sig2, n_per_group = n,
                         conf_level = 0.95),
    "below the alpha_lower critical value")
  expect_warning(
    r_F <- ci_srsnr(F_value = F_implied, df_1 = df1, df_2 = df2, N = N,
                    conf_level = 0.95),
    "below the alpha_lower critical value")
  expect_equal(r_design$value, r_F$value, tolerance = 1e-12)
})

test_that("ci_srsnr() rejects mixing F_value and means", {
  expect_error(
    ci_srsnr(F_value = 5, means = c(1, 2, 3),
             sigma_squared = 1, n_per_group = 5),
    "either 'F_value'"
  )
})

test_that("ci_srsnr() reports the lower-limit clamp once, in terms of the root signal-to-noise ratio, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_srsnr(F_value = 1.2, df_1 = 4, df_2 = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the square root of the signal-to-noise ratio is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$value[result$term == "lower_limit"], 0)
  expect_gt(result$value[result$term == "upper_limit"], 0)
})

test_that("ci_srsnr() names itself and the inputs when the noncentral root search fails", {
  expect_error(
    ci_srsnr(F_value = Inf, df_1 = 4, df_2 = 50, N = 55),
    "In ci_srsnr\\(\\).*F-statistic Inf"
  )
})
