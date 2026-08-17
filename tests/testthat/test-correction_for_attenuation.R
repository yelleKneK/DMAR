v <- function(tab, t) tab$value[tab$term == t]

test_that("correction_for_attenuation() applies the Spearman formula exactly", {
  res <- correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 0.70)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(v(res, "correlation_observed"), 0.30)
  expect_equal(v(res, "correlation_corrected"), 0.30 / sqrt(0.80 * 0.70))
  expect_identical(res$term, c("correlation_observed", "correlation_corrected",
                               "reliability_x", "reliability_y"))
})

test_that("correction_for_attenuation() with N disattenuates the Fisher's Z interval endpoints", {
  r <- 0.30; rxx <- 0.80; ryy <- 0.70; N <- 120
  res <- correction_for_attenuation(r, rxx, ryy, N = N)
  half <- qnorm(0.975) / sqrt(N - 3)
  expect_equal(v(res, "lower_limit"), tanh(atanh(r) - half) / sqrt(rxx * ryy))
  expect_equal(v(res, "upper_limit"), tanh(atanh(r) + half) / sqrt(rxx * ryy))
  expect_identical(attr(res, "conf_level"), 0.95)
  # The corrected estimate sits inside its corrected interval.
  expect_gt(v(res, "correlation_corrected"), v(res, "lower_limit"))
  expect_lt(v(res, "correlation_corrected"), v(res, "upper_limit"))
})

test_that("correction_for_attenuation() identities and limits", {
  # Perfect reliabilities change nothing.
  expect_equal(v(correction_for_attenuation(0.42, 1, 1), "correlation_corrected"), 0.42)
  # One-sided correction (error-free criterion).
  expect_equal(v(correction_for_attenuation(0.30, 0.64, 1), "correlation_corrected"),
               0.30 / 0.8)
  # The correction preserves sign.
  expect_equal(v(correction_for_attenuation(-0.30, 0.8, 0.7), "correlation_corrected"),
               -v(correction_for_attenuation(0.30, 0.8, 0.7), "correlation_corrected"))
})

test_that("correction_for_attenuation() warns, not truncates, when the correction exceeds 1", {
  expect_warning(res <- correction_for_attenuation(r = 0.9, reliability_x = 0.6,
                                            reliability_y = 0.6),
                 "exceeds 1")
  expect_equal(v(res, "correlation_corrected"), 0.9 / 0.6)  # 1.5, reported as is
})

test_that("correction_for_attenuation() validates its arguments", {
  expect_error(correction_for_attenuation(1.2, 0.8, 0.8), "\\[-1, 1\\]")
  expect_error(correction_for_attenuation(0.3, 0, 0.8), "\\(0, 1\\]")
  expect_error(correction_for_attenuation(0.3, 0.8, 1.2), "\\(0, 1\\]")
  expect_error(correction_for_attenuation(0.3, 0.8, 0.8, N = 2), "at least 4")
  expect_error(correction_for_attenuation(0.3, 0.8, 0.8, N = 100, conf_level = 1.1),
               "\\(0, 1\\)")
})
