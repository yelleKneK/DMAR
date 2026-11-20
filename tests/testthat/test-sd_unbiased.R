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

test_that("sd_unbiased() corrects upward by 1 / c4(N) at a very large N", {
  # c_N = 1 / c4(N), about 1 + 1/(4 (N - 1)). Through 1.0.x it was a
  # difference of lgamma() values that returned exactly s at N = 1e8 and
  # .998 s at N = 1e12. The pins are 1 / c4(N) from Rmpfr at 320 bits.
  N <- c(1e6, 1e8, 1e12)
  pins <- c(1.0000002500002811, 1.0000000025, 1.00000000000025)
  for (i in seq_along(N)) {
    value <- sd_unbiased(s = 1, N = N[i])$value
    expect_gt(value, 1)
    expect_equal(value, pins[i], tolerance = 2e-15, info = paste("N =", N[i]))
  }
  expect_equal(sd_unbiased(s = 2.5, N = 1e8)$value, 2.5 * pins[2],
               tolerance = 2e-15)
  # The correction never falls below 1 and never grows with N, out to and
  # past the N at which c4(N) rounds to 1 (about 4.5e15). Formed through
  # lbeta(), c4(N) came out above 1 at some N from about 1.7e14 on, and the
  # correction below 1.
  N <- c(2:200, 10^seq(2.5, 20, by = 0.05))
  value <- vapply(N, function(n) sd_unbiased(s = 1, N = n)$value, numeric(1))
  expect_true(all(value >= 1))
  expect_true(all(diff(value) <= 0))
  expect_identical(sd_unbiased(s = 1, N = 1e16)$value, 1)
})
