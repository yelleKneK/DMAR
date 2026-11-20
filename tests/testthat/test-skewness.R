test_that("skewness() returns a single numeric value", {
  result <- skewness(c(1, 2, 3, 4, 5))
  expect_type(result, "double")
  expect_length(result, 1L)
})

test_that("skewness() of a symmetric sample is approximately 0", {
  set.seed(113)
  x <- rnorm(10000)
  expect_equal(skewness(x), 0, tolerance = 0.1)
})

test_that("skewness() detects right-skewness", {
  set.seed(113)
  x <- rexp(1000, rate = 1)
  expect_gt(skewness(x), 1)
})

test_that("skewness() detects left-skewness", {
  set.seed(113)
  x <- -rexp(1000, rate = 1)
  expect_lt(skewness(x), -1)
})

test_that("skewness() handles NAs via na_rm", {
  x <- c(1, 2, 3, NA, 5)
  expect_true(is.na(skewness(x, na_rm = FALSE)))
  expect_false(is.na(skewness(x, na_rm = TRUE)))
})
