v <- function(tab, t) tab$value[tab$term == t]

test_that("moments_nc_t() returns a tidy dmar_tbl with the documented terms", {
  m <- moments_nc_t(df = 20, ncp = 2.5)
  expect_s3_class(m, "dmar_tbl")
  expect_s3_class(m, "data.frame")
  expect_identical(m$term, c("mean", "variance", "sd", "skewness",
                             "excess_kurtosis", "df", "ncp"))
  expect_type(m$value, "double")
  expect_equal(v(m, "df"), 20)
  expect_equal(v(m, "ncp"), 2.5)
  expect_equal(v(m, "sd"), sqrt(v(m, "variance")))
})

test_that("moments_nc_t() reduces to the central t at ncp = 0", {
  df <- 10
  m  <- moments_nc_t(df = df)
  expect_equal(v(m, "mean"), 0)
  expect_equal(v(m, "variance"), df / (df - 2))            # 1.25
  expect_equal(v(m, "skewness"), 0)
  expect_equal(v(m, "excess_kurtosis"), 6 / (df - 4))      # 1
})

test_that("moments_nc_t() mean matches the Owen (1968) closed form", {
  df <- 18; ncp <- 1.2
  closed <- ncp * sqrt(df / 2) * exp(lgamma((df - 1) / 2) - lgamma(df / 2))
  expect_equal(v(moments_nc_t(df, ncp), "mean"), closed)
})

test_that("moments_nc_t() returns NA for moments that do not exist", {
  # df = 2: the variance (and everything above it) is undefined; the mean is not.
  m2 <- moments_nc_t(df = 2, ncp = 1)
  expect_false(is.na(v(m2, "mean")))
  expect_true(is.na(v(m2, "variance")))
  expect_true(is.na(v(m2, "sd")))
  expect_true(is.na(v(m2, "skewness")))
  expect_true(is.na(v(m2, "excess_kurtosis")))
  # df = 4: the excess kurtosis needs more than 4 df, but the skewness exists.
  m4 <- moments_nc_t(df = 4, ncp = 1)
  expect_false(is.na(v(m4, "skewness")))
  expect_true(is.na(v(m4, "excess_kurtosis")))
  # A fractional df below 1 leaves even the mean undefined.
  expect_true(is.na(v(moments_nc_t(df = 0.7, ncp = 1), "mean")))
})

test_that("moments_nc_t() validates its arguments", {
  expect_error(moments_nc_t(df = 0, ncp = 1), "positive number")
  expect_error(moments_nc_t(df = -2, ncp = 1), "positive number")
  expect_error(moments_nc_t(df = c(5, 6), ncp = 1), "positive number")
  expect_error(moments_nc_t(df = 10, ncp = c(1, 2)), "single number")
  expect_error(moments_nc_t(df = 10, ncp = NA), "single number")
})

test_that("moments_nc_t() matches a Monte Carlo simulation", {
  skip_on_cran()
  set.seed(113)
  df <- 20; ncp <- 2.5; G <- 1e6
  x <- rt(G, df = df, ncp = ncp)
  m <- moments_nc_t(df, ncp)
  z <- (x - mean(x)) / sd(x)
  expect_lt(abs(v(m, "mean")     - mean(x)),     0.01)
  expect_lt(abs(v(m, "variance") - var(x)),      0.03)
  expect_lt(abs(v(m, "skewness") - mean(z^3)),   0.03)
})

# The reference values in moments_nc_t_reference.csv were computed with
# mpmath at 4000 digits from the raw moments, E(T^k) = E{(Z + ncp)^k}
# (df/2)^(k/2) Gamma((df - k)/2) / Gamma(df/2), with the central moments
# formed by subtraction at that precision. df, ncp, and the values are
# stored as hexadecimal doubles, which R reads exactly; its decimal parser
# is not exact at extreme magnitudes.
test_that("moments_nc_t() matches high-precision references across df and ncp", {
  skip_on_cran()
  ref <- utils::read.csv(test_path("moments_nc_t_reference.csv"),
                         colClasses = "character", na.strings = character(0))
  ref$df    <- as.numeric(ref$df)
  ref$ncp   <- as.numeric(ref$ncp)
  ref$exist <- ref$value != "NA"
  ref$value <- suppressWarnings(as.numeric(ref$value))
  expect_identical(nrow(ref), 355L)
  cells <- unique(ref[c("df", "ncp")])
  expect_identical(nrow(cells), 71L)
  for (i in seq_len(nrow(cells))) {
    m <- moments_nc_t(cells$df[i], cells$ncp[i])
    r <- ref[ref$df == cells$df[i] & ref$ncp == cells$ncp[i], ]
    got <- m$value[match(r$term, m$term)]
    label <- sprintf("df = %.17g, ncp = %g", cells$df[i], cells$ncp[i])
    # NA exactly where the moment does not exist.
    expect_identical(is.na(got), !r$exist, info = label)
    zero <- r$exist & r$value == 0
    expect_identical(got[zero], r$value[zero], info = label)
    rest <- r$exist & r$value != 0
    if (any(rest)) {
      expect_lt(max(abs(got[rest] - r$value[rest]) / abs(r$value[rest])),
                1e-11, label = paste("largest relative error at", label))
    }
  }
})

test_that("moments_nc_t() at df = Inf is the normal distribution", {
  expect_identical(moments_nc_t(Inf, 2)$value[1:5], c(2, 1, 1, 0, 0))
  expect_identical(moments_nc_t(Inf)$value[1:5], c(0, 1, 1, 0, 0))
})

test_that("moments_nc_t() keeps the sign of the excess kurtosis at a large df", {
  # The difference of raw moments gave -0.0022675 here.
  k <- v(moments_nc_t(1e4, 100), "excess_kurtosis")
  expect_gt(k, 0)
  expect_lt(abs(k - 0.00133400558) / 0.00133400558, 1e-9)
})

test_that("moments_nc_t() has a positive variance and a finite sd at any df", {
  skip_on_cran()
  for (ncp in c(0, 1, 1e3)) {
    for (df in 10^seq(0.31, 300, length.out = 200)) {
      m <- moments_nc_t(df, ncp)
      expect_gt(v(m, "variance"), 0)
      expect_true(is.finite(v(m, "sd")))
    }
  }
})

test_that("moments_nc_t() mirrors the distribution with the sign of ncp", {
  for (df in c(7.3, 1e6)) {
    p <- moments_nc_t(df, 2.5)
    n <- moments_nc_t(df, -2.5)
    expect_identical(v(n, "mean"), -v(p, "mean"))
    expect_identical(v(n, "skewness"), -v(p, "skewness"))
    expect_gt(v(p, "mean"), 0)
    expect_gt(v(p, "skewness"), 0)
    for (term in c("variance", "sd", "excess_kurtosis")) {
      expect_identical(v(n, term), v(p, term))
    }
  }
})

test_that("moments_nc_t() stays finite at a very large ncp and refuses an infinite one", {
  expect_true(all(is.finite(moments_nc_t(20, 1e80)$value[1:5])))
  expect_error(moments_nc_t(20, Inf), "finite")
  expect_error(moments_nc_t(20, -Inf), "finite")
})

test_that("the reference grid holds the series branch between df 100 and 1e4", {
  # The expansions in 1 / df are used from df = 100 on. With anchors only at
  # 100, 100.0001, and 1e4, a branch point moved to 1000 or 1e4 passed; the
  # closed forms lose digits in that range (2.5e-10 in the excess kurtosis).
  ref <- utils::read.csv(test_path("moments_nc_t_reference.csv"),
                         colClasses = "character", na.strings = character(0))
  df <- unique(as.numeric(ref$df))
  expect_true(all(c(150, 1000.5, 5000) %in% df))
  expect_gte(sum(df > 100.001 & df < 1e4), 3L)
})

test_that("moments_nc_t() reports an overflowing standard deviation as Inf", {
  # abs(ncp) times the standard deviation factor is above the largest
  # double; the variance and the standard deviation were NaN (Inf / Inf).
  val <- function(df, ncp) {
    m <- moments_nc_t(df, ncp)
    stats::setNames(m$value, m$term)
  }
  for (cell in list(c(2.2, 9e307), c(2.5, 1.7e308), c(2.01, 2e307))) {
    m <- val(cell[1], cell[2])
    expect_identical(unname(m[c("variance", "sd")]), c(Inf, Inf))
  }
  # The skewness has a finite limit as ncp grows, and keeps it where the
  # standard deviation overflows (which needs df below about 3.1, where the
  # standard deviation exceeds abs(ncp)).
  far  <- val(3.05, 1.79e308)
  near <- val(3.05, 1e150)
  expect_identical(unname(far["sd"]), Inf)
  expect_true(is.finite(near[["sd"]]))
  expect_equal(far[["skewness"]], near[["skewness"]], tolerance = 1e-12)
  expect_equal(val(3.05, -1.79e308)[["skewness"]], -near[["skewness"]],
               tolerance = 1e-12)
})
