v <- function(tab, t) tab$value[tab$term == t]

test_that("moments_nc_F() returns a tidy dmar_tbl with the documented terms", {
  m <- moments_nc_F(df_1 = 3, df_2 = 40, ncp = 8)
  expect_s3_class(m, "dmar_tbl")
  expect_identical(m$term, c("mean", "variance", "sd", "skewness",
                             "excess_kurtosis", "df_1", "df_2", "ncp"))
  expect_type(m$value, "double")
  expect_equal(v(m, "df_1"), 3)
  expect_equal(v(m, "df_2"), 40)
  expect_equal(v(m, "ncp"), 8)
  expect_equal(v(m, "sd"), sqrt(v(m, "variance")))
})

test_that("moments_nc_F() reduces to the central F at ncp = 0", {
  d1 <- 3; d2 <- 40
  m  <- moments_nc_F(df_1 = d1, df_2 = d2)
  expect_equal(v(m, "mean"), d2 / (d2 - 2))
  # Central F variance: 2 d2^2 (d1 + d2 - 2) / [d1 (d2 - 2)^2 (d2 - 4)].
  expect_equal(v(m, "variance"),
               2 * d2^2 * (d1 + d2 - 2) / (d1 * (d2 - 2)^2 * (d2 - 4)))
})

test_that("moments_nc_F() mean and variance match the closed forms", {
  d1 <- 3; d2 <- 40; lam <- 8
  mean_cf <- d2 * (d1 + lam) / (d1 * (d2 - 2))
  var_cf  <- 2 * (d2 / d1)^2 *
    ((d1 + lam)^2 + (d1 + 2 * lam) * (d2 - 2)) /
    ((d2 - 2)^2 * (d2 - 4))
  m <- moments_nc_F(d1, d2, lam)
  expect_equal(v(m, "mean"), mean_cf)
  expect_equal(v(m, "variance"), var_cf)
})

test_that("moments_nc_F() returns NA for moments that do not exist", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # df_2 = 4: the variance needs more than 4 denominator df; the mean exists.
  m4 <- moments_nc_F(df_1 = 2, df_2 = 4, ncp = 5)
  expect_false(is.na(v(m4, "mean")))
  expect_equal(v(m4, "mean"), 4 * (2 + 5) / (2 * (4 - 2)))   # = 7
  expect_true(is.na(v(m4, "variance")))
  # df_2 = 5: variance exists, but skewness (needs > 6) and kurtosis do not.
  m5 <- moments_nc_F(df_1 = 3, df_2 = 5, ncp = 2)
  expect_false(is.na(v(m5, "variance")))
  expect_true(is.na(v(m5, "skewness")))
  expect_true(is.na(v(m5, "excess_kurtosis")))
})

test_that("moments_nc_F() validates its arguments", {
  expect_error(moments_nc_F(df_1 = 0, df_2 = 10), "positive number")
  expect_error(moments_nc_F(df_1 = 3, df_2 = -1), "positive number")
  expect_error(moments_nc_F(df_1 = 3, df_2 = 10, ncp = -1), "non-negative")
  expect_error(moments_nc_F(df_1 = c(2, 3), df_2 = 10), "positive number")
})

test_that("moments_nc_F() matches a Monte Carlo simulation", {
  skip_on_cran()
  set.seed(113)
  d1 <- 3; d2 <- 40; lam <- 8; G <- 1e6
  y <- rf(G, df1 = d1, df2 = d2, ncp = lam)
  m <- moments_nc_F(d1, d2, lam)
  z <- (y - mean(y)) / sd(y)
  expect_lt(abs(v(m, "mean")     - mean(y)),    0.02)
  expect_lt(abs(v(m, "variance") - var(y)),     0.10)
  expect_lt(abs(v(m, "skewness") - mean(z^3)),  0.05)
})

test_that("moments_nc_F() matches exact rational arithmetic, large df included", {
  # Each row: df_1, df_2, ncp, then the mean, variance, skewness, and
  # excess kurtosis from exact rational arithmetic on the raw moments
  # (central moments by exact subtraction), rounded to 20 digits.
  ref <- rbind(
    c(1, 9, 0, 1.2857142857142857143, 5.2897959183673469388,
      6.7082039324993690892, 214.5),
    c(3, 12, 2, 2, 3.8, 3.1679331573155635228, 25.337950138504155125),
    c(5, 20, 7.5, 2.7777777777777777778, 3.1867283950617283951,
      1.8130544855365979731, 6.6898686162198289255),
    c(10, 1000, 0, 1.0020040080160320641, 0.2032217125174365439,
      0.91055417553121950508, 1.2582675198657382081),
    c(1e6, 1e6, 0, 1.000002000004000008, 4.000028000144000656e-6,
      0.0060000230001282507474, 0.000066000618005562048186),
    c(3, 1e6, 1e6, 333335.00000333334, 666673.11115000019289,
      0.0059876636390060254182, 0.000065333932894327602611),
    c(2, 1e9, 5, 3.500000007000000014, 6.000000060500000364,
      1.1567035031793417224, 1.8333333966944451011),
    c(1e7, 1000, 1e5, 1.0120240480961923848, 0.002056817024766710748,
      0.17960493868962514719, 0.060581229431996648115),
    c(1, 1e12, 0, 1.000000000002, 2.000000000014, 2.8284271247589180197,
      12.000000000144),
    c(50, 9, 300, 9, 33.603428571428571429, 4.2164846577817891602,
      92.004325751883978781),
    c(2.5, 8.5, 0.25, 1.4384615384615384615, 3.2909138724523339908,
      5.792403147796932049, 275.30280709801641696),
    c(4, 60, 500, 130.34482758620689655, 745.87863088160353321,
      0.78435601019801162851, 1.1890601192229498535),
    c(7, 33, 10000, 1521.8018433179723502, 160705.12063949337315,
      1.1282627026606071425, 2.5422222964112758136),
    c(1000, 1e6, 1, 1.001002002004004008, 0.002006016042088216417,
      0.089576962798088043746, 0.012048060659923853366))
  terms <- c("mean", "variance", "skewness", "excess_kurtosis")
  for (i in seq_len(nrow(ref))) {
    m <- moments_nc_F(ref[i, 1], ref[i, 2], ref[i, 3])
    got <- m$value[match(terms, m$term)]
    expect_lt(max(abs(got - ref[i, 4:7]) / abs(ref[i, 4:7])), 2e-15,
              label = paste("relative error, row", i))
  }
  # The case the differences of raw moments lost through 1.0.x (they gave
  # -5.29e-05).
  expect_equal(v(moments_nc_F(1e6, 1e6, 0), "excess_kurtosis"),
               6.60006e-05, tolerance = 1e-6)
})

test_that("moments_nc_F() agrees with the differences of raw moments at ordinary df", {
  skip_on_cran()
  # The form used through 1.0.x, which is accurate while the degrees of
  # freedom are moderate.
  raw_form <- function(df_1, df_2, lam) {
    k1 <- df_1 + lam; k2 <- 2 * (df_1 + 2 * lam)
    k3 <- 8 * (df_1 + 3 * lam); k4 <- 48 * (df_1 + 4 * lam)
    ex <- c(k1, k2 + k1^2, k3 + 3 * k2 * k1 + k1^3,
            k4 + 4 * k3 * k1 + 3 * k2^2 + 6 * k2 * k1^2 + k1^4)
    raw <- vapply(1:4, function(k) if (df_2 > 2 * k)
      (df_2 / df_1)^k * ex[k] * prod(1 / (df_2 - 2 * seq_len(k))) else NA_real_,
      0)
    va <- raw[2] - raw[1]^2
    c(raw[1], va, sqrt(va),
      (raw[3] - 3 * raw[1] * raw[2] + 2 * raw[1]^3) / va^1.5,
      (raw[4] - 4 * raw[1] * raw[3] + 6 * raw[1]^2 * raw[2] -
         3 * raw[1]^4) / va^2 - 3)
  }
  g <- expand.grid(df_1 = c(1, 2.5, 3, 10, 60, 200),
                   df_2 = c(1, 2, 3, 4, 4.5, 6, 7, 8, 8.5, 9, 12, 40, 200),
                   ncp = c(0, 0.5, 8, 50))
  worst <- 0
  for (i in seq_len(nrow(g))) {
    m <- moments_nc_F(g$df_1[i], g$df_2[i], g$ncp[i])$value[1:5]
    o <- raw_form(g$df_1[i], g$df_2[i], g$ncp[i])
    expect_identical(is.na(m), is.na(o), info = paste("row", i))
    ok <- !is.na(m)
    if (any(ok)) worst <- max(worst, abs(m[ok] - o[ok]) / abs(o[ok]))
  }
  expect_lt(worst, 1e-10)
})
