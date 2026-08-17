test_that("cv() returns a data frame with term='cv'", {
  result <- cv(mean = 100, sd = 15)
  expect_s3_class(result, "data.frame")
  expect_equal(result$term, "cv")
})

test_that("cv() basic calculation: sd/mean", {
  result <- cv(mean = 100, sd = 15)
  expect_equal(result$value, 0.15)
})

test_that("cv() unbiased estimate is larger than biased for small N", {
  biased   <- cv(mean = 100, sd = 15, N = 10)
  unbiased <- cv(mean = 100, sd = 15, N = 10, unbiased = TRUE)
  expect_gt(unbiased$value, biased$value)
})

test_that("cv() errors when mean and sd given alongside cv", {
  expect_error(cv(mean = 100, sd = 15, cv = 0.15), "do not specify")
})

test_that("cv() from cv= argument gives same result as from mean/sd", {
  r_from_raw <- cv(mean = 50, sd = 10)
  r_from_cv  <- cv(cv = 0.2)
  expect_equal(r_from_raw$value, r_from_cv$value)
})
