test_that("sd_unbiased() returns the expected 1-row tidy data frame", {
  result <- sd_unbiased(s = 2, N = 10)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "sd")
})

test_that("sd_unbiased() value is greater than the sample s for finite N", {
  result <- sd_unbiased(s = 2, N = 10)
  expect_gt(result$value, 2)
})

test_that("sd_unbiased() approaches s as N grows large", {
  small <- sd_unbiased(s = 2, N = 5)$value
  large <- sd_unbiased(s = 2, N = 1000)$value
  expect_gt(small - 2, large - 2)
  expect_equal(large, 2, tolerance = 1e-3)
})

test_that("sd_unbiased() accepts a raw vector X in place of s + N", {
  x <- c(1, 2, 3, 4, 5)
  r1 <- sd_unbiased(s = sd(x), N = length(x))
  r2 <- sd_unbiased(X = x)
  expect_equal(r1$value, r2$value, tolerance = 1e-10)
})
