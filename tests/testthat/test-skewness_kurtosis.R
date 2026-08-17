test_that("skewness() is zero for symmetric data", {
  expect_equal(skewness(1:5), 0)
  expect_equal(skewness(c(-3, -2, -1, 0, 1, 2, 3)), 0)
})

test_that("kurtosis() of 1:5 equals -1.2 (SAS/SPSS Type 2)", {
  expect_equal(kurtosis(1:5), -1.2)
})

test_that("skewness() detects asymmetry: positive for right-tail data", {
  set.seed(113)
  expect_gt(skewness(rexp(2000, rate = 1)), 1.0)
})

test_that("kurtosis() detects heavy tails: positive for t(4)", {
  set.seed(113)
  expect_gt(kurtosis(rt(2000, df = 4)), 1.0)
})

test_that("skewness() and kurtosis() return NA for too-small samples", {
  expect_true(is.na(skewness(c(1, 2))))
  expect_true(is.na(kurtosis(c(1, 2, 3))))
})

test_that("skewness() and kurtosis() return NA when sd = 0", {
  expect_true(is.na(skewness(rep(3, 10))))
  expect_true(is.na(kurtosis(rep(3, 10))))
})

test_that("na_rm = TRUE strips NAs; na_rm = FALSE returns NA", {
  x <- c(1, 2, 3, 4, 5, NA)
  expect_equal(skewness(x, na_rm = TRUE), 0)
  expect_true(is.na(skewness(x, na_rm = FALSE)))

  expect_equal(kurtosis(x, na_rm = TRUE), -1.2)
  expect_true(is.na(kurtosis(x, na_rm = FALSE)))
})

test_that("skewness() and kurtosis() reject non-numeric input", {
  expect_error(skewness("a"),  "numeric")
  expect_error(kurtosis(TRUE), "numeric")
})

test_that("descriptives() still works after the promotion", {
  res <- descriptives(data.frame(v = c(1, 2, 3, 4, 5)))
  expect_equal(res$descriptives$skewness, 0)
  expect_equal(res$descriptives$kurtosis, -1.2)
})
