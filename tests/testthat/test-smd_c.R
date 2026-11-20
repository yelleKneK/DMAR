test_that("smd_c() returns a 1-row tidy data frame with term smd_c", {
  result <- smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 30)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "smd_c")
})

test_that("smd_c() returns the textbook value (mean_T - mean_C) / s_C", {
  result <- smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 30)
  expect_equal(result$value, (5 - 4) / 1.2, tolerance = 1e-12)
})

test_that("smd_c() with unbiased = TRUE applies the Hedges correction (smaller in magnitude)", {
  biased   <- smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 30, unbiased = FALSE)$value
  unbiased <- smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 30, unbiased = TRUE)$value
  expect_lt(abs(unbiased), abs(biased))
})

test_that("smd_c() refuses the Hedges correction on 1 degree of freedom", {
  # J(1) = 0, so 1.0.x returned 0 whatever the data.
  err <- expect_error(
    smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 2, unbiased = TRUE),
    class = "dmar_hedges_correction_undefined")
  expect_s3_class(err, "dmar_refusal")
  expect_identical(err$df, 1)
  expect_match(conditionMessage(err), "smd_c(unbiased = TRUE)", fixed = TRUE)
  expect_match(conditionMessage(err), "n_C = 2", fixed = TRUE)
  expect_match(conditionMessage(err), "unbiased = FALSE", fixed = TRUE)
  expect_error(smd_c(group_T = c(1, 4, 6), group_C = c(2, 3), unbiased = TRUE),
               class = "dmar_hedges_correction_undefined")
  # The uncorrected estimate is returned, and 2 degrees of freedom are
  # enough for the corrected one: J(2) = Gamma(1) / Gamma(1/2).
  expect_equal(smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 2)$value,
               1 / 1.2)
  expect_equal(
    smd_c(mean_T = 5, mean_C = 4, s_C = 1.2, n_C = 3, unbiased = TRUE)$value,
    (1 / 1.2) * gamma(1) / gamma(1 / 2), tolerance = 1e-12)
})
