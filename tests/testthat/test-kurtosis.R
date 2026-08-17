test_that("kurtosis() returns a single numeric value", {
  result <- kurtosis(c(1, 2, 3, 4, 5, 5, 5))
  expect_type(result, "double")
  expect_length(result, 1L)
})

test_that("kurtosis() of a normal sample is approximately 0 (excess kurtosis)", {
  set.seed(113)
  x <- rnorm(10000)
  expect_equal(kurtosis(x), 0, tolerance = 0.1)
})

test_that("kurtosis() of a heavy-tailed t(3) sample is positive", {
  set.seed(113)
  x <- rt(10000, df = 3)
  expect_gt(kurtosis(x), 1)
})

test_that("kurtosis() equals the hand-computed Type 2 value on a fixed vector", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  n <- length(x)
  z <- (x - mean(x)) / sd(x)
  term1      <- (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))
  correction <- 3 * (n - 1)^2 / ((n - 2) * (n - 3))
  hand <- term1 * sum(z^4) - correction  # 0.940625
  expect_equal(kurtosis(x), hand, tolerance = 1e-12)
  expect_equal(kurtosis(x), 0.940625, tolerance = 1e-9)
})

test_that("kurtosis() handles NAs via na_rm", {
  x <- c(1, 2, 3, NA, 5)
  expect_true(is.na(kurtosis(x, na_rm = FALSE)))
  expect_false(is.na(kurtosis(x, na_rm = TRUE)))
})
