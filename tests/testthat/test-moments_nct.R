v <- function(tab, t) tab$value[tab$term == t]

test_that("moments_nct() returns a tidy dmar_tbl with the documented terms", {
  m <- moments_nct(df = 20, ncp = 2.5)
  expect_s3_class(m, "dmar_tbl")
  expect_s3_class(m, "data.frame")
  expect_identical(m$term, c("mean", "variance", "sd", "skewness",
                             "excess_kurtosis", "df", "ncp"))
  expect_type(m$value, "double")
  expect_equal(v(m, "df"), 20)
  expect_equal(v(m, "ncp"), 2.5)
  expect_equal(v(m, "sd"), sqrt(v(m, "variance")))
})

test_that("moments_nct() reduces to the central t at ncp = 0", {
  df <- 10
  m  <- moments_nct(df = df)
  expect_equal(v(m, "mean"), 0)
  expect_equal(v(m, "variance"), df / (df - 2))            # 1.25
  expect_equal(v(m, "skewness"), 0)
  expect_equal(v(m, "excess_kurtosis"), 6 / (df - 4))      # 1
})

test_that("moments_nct() mean matches the Owen (1968) closed form", {
  df <- 18; ncp <- 1.2
  closed <- ncp * sqrt(df / 2) * exp(lgamma((df - 1) / 2) - lgamma(df / 2))
  expect_equal(v(moments_nct(df, ncp), "mean"), closed)
})

test_that("moments_nct() returns NA for moments that do not exist", {
  # df = 2: the variance (and everything above it) is undefined; the mean is not.
  m2 <- moments_nct(df = 2, ncp = 1)
  expect_false(is.na(v(m2, "mean")))
  expect_true(is.na(v(m2, "variance")))
  expect_true(is.na(v(m2, "sd")))
  expect_true(is.na(v(m2, "skewness")))
  expect_true(is.na(v(m2, "excess_kurtosis")))
  # df = 4: the excess kurtosis needs more than 4 df, but the skewness exists.
  m4 <- moments_nct(df = 4, ncp = 1)
  expect_false(is.na(v(m4, "skewness")))
  expect_true(is.na(v(m4, "excess_kurtosis")))
  # A fractional df below 1 leaves even the mean undefined.
  expect_true(is.na(v(moments_nct(df = 0.7, ncp = 1), "mean")))
})

test_that("moments_nct() validates its arguments", {
  expect_error(moments_nct(df = 0, ncp = 1), "positive number")
  expect_error(moments_nct(df = -2, ncp = 1), "positive number")
  expect_error(moments_nct(df = c(5, 6), ncp = 1), "positive number")
  expect_error(moments_nct(df = 10, ncp = c(1, 2)), "single number")
  expect_error(moments_nct(df = 10, ncp = NA), "single number")
})

test_that("moments_nct() matches a Monte Carlo simulation", {
  skip_on_cran()
  set.seed(113)
  df <- 20; ncp <- 2.5; G <- 1e6
  x <- rt(G, df = df, ncp = ncp)
  m <- moments_nct(df, ncp)
  z <- (x - mean(x)) / sd(x)
  expect_lt(abs(v(m, "mean")     - mean(x)),     0.01)
  expect_lt(abs(v(m, "variance") - var(x)),      0.03)
  expect_lt(abs(v(m, "skewness") - mean(z^3)),   0.03)
})
