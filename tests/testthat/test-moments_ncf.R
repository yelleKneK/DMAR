v <- function(tab, t) tab$value[tab$term == t]

test_that("moments_ncf() returns a tidy dmar_tbl with the documented terms", {
  m <- moments_ncf(df_1 = 3, df_2 = 40, ncp = 8)
  expect_s3_class(m, "dmar_tbl")
  expect_identical(m$term, c("mean", "variance", "sd", "skewness",
                             "excess_kurtosis", "df_1", "df_2", "ncp"))
  expect_type(m$value, "double")
  expect_equal(v(m, "df_1"), 3)
  expect_equal(v(m, "df_2"), 40)
  expect_equal(v(m, "ncp"), 8)
  expect_equal(v(m, "sd"), sqrt(v(m, "variance")))
})

test_that("moments_ncf() reduces to the central F at ncp = 0", {
  d1 <- 3; d2 <- 40
  m  <- moments_ncf(df_1 = d1, df_2 = d2)
  expect_equal(v(m, "mean"), d2 / (d2 - 2))
  # Central F variance: 2 d2^2 (d1 + d2 - 2) / [d1 (d2 - 2)^2 (d2 - 4)].
  expect_equal(v(m, "variance"),
               2 * d2^2 * (d1 + d2 - 2) / (d1 * (d2 - 2)^2 * (d2 - 4)))
})

test_that("moments_ncf() mean and variance match the closed forms", {
  d1 <- 3; d2 <- 40; lam <- 8
  mean_cf <- d2 * (d1 + lam) / (d1 * (d2 - 2))
  var_cf  <- 2 * (d2 / d1)^2 *
    ((d1 + lam)^2 + (d1 + 2 * lam) * (d2 - 2)) /
    ((d2 - 2)^2 * (d2 - 4))
  m <- moments_ncf(d1, d2, lam)
  expect_equal(v(m, "mean"), mean_cf)
  expect_equal(v(m, "variance"), var_cf)
})

test_that("moments_ncf() returns NA for moments that do not exist", {
  # df_2 = 4: the variance needs more than 4 denominator df; the mean exists.
  m4 <- moments_ncf(df_1 = 2, df_2 = 4, ncp = 5)
  expect_false(is.na(v(m4, "mean")))
  expect_equal(v(m4, "mean"), 4 * (2 + 5) / (2 * (4 - 2)))   # = 7
  expect_true(is.na(v(m4, "variance")))
  # df_2 = 5: variance exists, but skewness (needs > 6) and kurtosis do not.
  m5 <- moments_ncf(df_1 = 3, df_2 = 5, ncp = 2)
  expect_false(is.na(v(m5, "variance")))
  expect_true(is.na(v(m5, "skewness")))
  expect_true(is.na(v(m5, "excess_kurtosis")))
})

test_that("moments_ncf() validates its arguments", {
  expect_error(moments_ncf(df_1 = 0, df_2 = 10), "positive number")
  expect_error(moments_ncf(df_1 = 3, df_2 = -1), "positive number")
  expect_error(moments_ncf(df_1 = 3, df_2 = 10, ncp = -1), "non-negative")
  expect_error(moments_ncf(df_1 = c(2, 3), df_2 = 10), "positive number")
})

test_that("moments_ncf() matches a Monte Carlo simulation", {
  skip_on_cran()
  set.seed(113)
  d1 <- 3; d2 <- 40; lam <- 8; G <- 1e6
  y <- rf(G, df1 = d1, df2 = d2, ncp = lam)
  m <- moments_ncf(d1, d2, lam)
  z <- (y - mean(y)) / sd(y)
  expect_lt(abs(v(m, "mean")     - mean(y)),    0.02)
  expect_lt(abs(v(m, "variance") - var(y)),     0.10)
  expect_lt(abs(v(m, "skewness") - mean(z^3)),  0.05)
})
