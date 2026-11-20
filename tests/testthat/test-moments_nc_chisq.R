v <- function(tab, t) tab$value[tab$term == t]

test_that("moments_nc_chisq() returns a tidy dmar_tbl with the documented terms", {
  m <- moments_nc_chisq(df = 5, ncp = 3)
  expect_s3_class(m, "dmar_tbl")
  expect_identical(m$term, c("mean", "variance", "sd", "skewness",
                             "excess_kurtosis", "df", "ncp"))
  expect_type(m$value, "double")
  expect_equal(v(m, "df"), 5)
  expect_equal(v(m, "ncp"), 3)
  expect_equal(v(m, "sd"), sqrt(v(m, "variance")))
})

test_that("moments_nc_chisq() matches the cumulant closed forms", {
  df <- 5; lam <- 3
  m <- moments_nc_chisq(df, lam)
  expect_equal(v(m, "mean"), df + lam)
  expect_equal(v(m, "variance"), 2 * (df + 2 * lam))
  expect_equal(v(m, "skewness"),
               sqrt(8) * (df + 3 * lam) / (df + 2 * lam)^(3 / 2))
  expect_equal(v(m, "excess_kurtosis"),
               12 * (df + 4 * lam) / (df + 2 * lam)^2)
})

test_that("moments_nc_chisq() reduces to the central chi square at ncp = 0", {
  df <- 5
  m  <- moments_nc_chisq(df = df)
  expect_equal(v(m, "mean"), df)
  expect_equal(v(m, "variance"), 2 * df)
  expect_equal(v(m, "skewness"), sqrt(8 / df))
  expect_equal(v(m, "excess_kurtosis"), 12 / df)
})

test_that("moments_nc_chisq() moments always exist (never NA)", {
  # Unlike the t and F, the chi square has finite moments of every order.
  expect_false(anyNA(moments_nc_chisq(df = 1, ncp = 10)$value))
  expect_false(anyNA(moments_nc_chisq(df = 0.5, ncp = 0)$value))
})

test_that("moments_nc_chisq() validates its arguments", {
  expect_error(moments_nc_chisq(df = 0, ncp = 1), "positive number")
  expect_error(moments_nc_chisq(df = -3, ncp = 1), "positive number")
  expect_error(moments_nc_chisq(df = 5, ncp = -1), "non-negative")
  expect_error(moments_nc_chisq(df = c(5, 6), ncp = 1), "positive number")
})

test_that("moments_nc_chisq() matches a Monte Carlo simulation", {
  skip_on_cran()
  set.seed(113)
  df <- 5; lam <- 3; G <- 1e6
  x <- rchisq(G, df = df, ncp = lam)
  m <- moments_nc_chisq(df, lam)
  z <- (x - mean(x)) / sd(x)
  expect_lt(abs(v(m, "mean")     - mean(x)),    0.02)
  expect_lt(abs(v(m, "variance") - var(x)),     0.15)
  expect_lt(abs(v(m, "skewness") - mean(z^3)),  0.03)
})
