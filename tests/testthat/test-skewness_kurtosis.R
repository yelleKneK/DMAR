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

test_that("skewness() is G1 and differs from g1 and b1 by the factors on its help page", {
  # Oracle: the sample moments m_r = sum((x - mean)^r) / n computed here,
  # g1 = m3 / m2^(3/2), G1 = g1 sqrt(n(n - 1)) / (n - 2),
  # b1 = g1 ((n - 1) / n)^(3/2) = m3 / s^3.
  x <- c(2.1, 3.4, 3.9, 5.0, 7.7, 12.6, 4.4, 3.1)
  n <- length(x)
  m2 <- mean((x - mean(x))^2); m3 <- mean((x - mean(x))^3)
  g1 <- m3 / m2^1.5
  expect_equal(skewness(x), g1 * sqrt(n * (n - 1)) / (n - 2),
               tolerance = 1e-12)
  expect_equal(g1 * ((n - 1) / n)^1.5, m3 / sd(x)^3, tolerance = 1e-12)
  expect_equal(descriptives(data.frame(v = x))$descriptives$skewness,
               skewness(x), tolerance = 1e-12)
})
