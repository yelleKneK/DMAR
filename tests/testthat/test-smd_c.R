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
