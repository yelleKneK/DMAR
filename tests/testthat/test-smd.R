test_that("smd() returns a data frame with term and value columns", {
  result <- smd(mean_1 = 1, mean_2 = 0, s = 1)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "smd")
})

test_that("smd() from summary data matches raw data calculation", {
  set.seed(113)
  g1 <- rnorm(25, mean = 0.5, sd = 1)
  g2 <- rnorm(25, mean = 0,   sd = 1)
  result_raw     <- smd(group_1 = g1, group_2 = g2)
  result_summary <- smd(
    mean_1 = mean(g1), mean_2 = mean(g2),
    s_1 = sd(g1), s_2 = sd(g2), n_1 = 25, n_2 = 25
  )
  expect_equal(result_raw$value, result_summary$value, tolerance = 1e-10)
})

test_that("smd() biased estimate is correct for known inputs", {
  # (1 - 0) / 1 = 1
  result <- smd(mean_1 = 1, mean_2 = 0, s = 1)
  expect_equal(result$value, 1)
})

test_that("smd() unbiased estimate is smaller in magnitude than biased for small samples", {
  result_biased   <- smd(mean_1 = 1, mean_2 = 0, s_1 = 1, s_2 = 1, n_1 = 10, n_2 = 10)
  result_unbiased <- smd(mean_1 = 1, mean_2 = 0, s_1 = 1, s_2 = 1, n_1 = 10, n_2 = 10, unbiased = TRUE)
  expect_lt(abs(result_unbiased$value), abs(result_biased$value))
})

test_that("smd() errors when both raw data and summary values are provided", {
  expect_error(smd(group_1 = c(1, 2, 3), group_2 = c(0, 1, 2), mean_1 = 2), "raw data")
})

test_that("smd() errors when n is needed for unbiased but not given", {
  expect_error(
    smd(mean_1 = 1, mean_2 = 0, s = 1, unbiased = TRUE),
    "sample sizes"
  )
})

test_that("smd() unbiased correction does not overflow at large samples (gamma -> lgamma)", {
  # gamma() overflows above df/2 ~ 171.6; the log-scale factor must stay finite.
  J <- function(df) exp(lgamma(df / 2) - lgamma((df - 1) / 2)) / sqrt(df / 2)
  for (n in c(173, 200, 500)) {
    v <- smd(mean_1 = 1, mean_2 = 0, s_1 = 1, s_2 = 1, n_1 = n, n_2 = n,
             unbiased = TRUE)$value
    expect_true(is.finite(v))
    expect_equal(v, J(2 * n - 2), tolerance = 1e-12)
  }
  # Small df is unchanged relative to the exact factor.
  expect_equal(
    smd(mean_1 = 1, mean_2 = 0, s_1 = 1, s_2 = 1, n_1 = 10, n_2 = 10,
        unbiased = TRUE)$value,
    J(18), tolerance = 1e-12)
})

test_that("smd() errors on unresolvable input rather than returning NULL", {
  expect_error(smd(mean_1 = 1, s = 2))
  expect_error(smd(group_1 = c(1, 2, 3)))
  expect_error(smd(mean_1 = 1, mean_2 = 0))
})
