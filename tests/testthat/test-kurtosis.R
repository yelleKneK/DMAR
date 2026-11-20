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

test_that("kurtosis() is G2 and relates to g2 and b2 as its help page says", {
  # Oracle: the sample moments m_r = sum((x - mean)^r) / n computed here,
  # g2 = m4 / m2^2 - 3, G2 = ((n + 1) g2 + 6)(n - 1) / ((n - 2)(n - 3)),
  # b2 = (g2 + 3)(1 - 1/n)^2 - 3 = m4 / s^4 - 3.
  x <- c(2.1, 3.4, 3.9, 5.0, 7.7, 12.6, 4.4, 3.1)
  n <- length(x)
  m2 <- mean((x - mean(x))^2); m4 <- mean((x - mean(x))^4)
  g2 <- m4 / m2^2 - 3
  expect_equal(kurtosis(x), ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3)),
               tolerance = 1e-12)
  expect_equal((g2 + 3) * (1 - 1 / n)^2 - 3, m4 / sd(x)^4 - 3, tolerance = 1e-12)
  expect_equal(descriptives(data.frame(v = x))$descriptives$kurtosis,
               kurtosis(x), tolerance = 1e-12)
})

test_that("kurtosis() has expected value 0 for normal samples, as its help page says", {
  # The help page says that subtracting 3(n - 1)^2 / ((n - 2)(n - 3)) makes
  # the expected value 0 under normality (it said "asymptotic correction"
  # before 1.1.0; the subtraction is exact in finite samples). Monte Carlo
  # at n = 10: 20,000 samples (all usable), the mean of G2 within 4 Monte
  # Carlo standard errors of 0, and the mean of g2 near its exact
  # expectation -6 / (n + 1), which the subtraction removes.
  set.seed(2026)
  n <- 10L; R <- 20000L
  X <- matrix(rnorm(n * R), n, R)
  G2 <- apply(X, 2L, kurtosis)
  expect_identical(sum(!is.finite(G2)), 0L)
  expect_lt(abs(mean(G2)), 4 * sd(G2) / sqrt(R))
  g2 <- ((G2 * (n - 2) * (n - 3)) / (n - 1) - 6) / (n + 1)
  expect_lt(abs(mean(g2) + 6 / (n + 1)), 4 * sd(g2) / sqrt(R))
})
