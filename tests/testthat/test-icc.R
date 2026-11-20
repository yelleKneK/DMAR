# Shrout & Fleiss (1979) Table 1 ratings
sf_data <- matrix(
  c(9, 2, 5, 8,
    6, 1, 3, 2,
    8, 4, 6, 8,
    7, 1, 2, 6,
    10, 5, 6, 9,
    6, 2, 4, 7),
  nrow = 6, byrow = TRUE
)

test_that("icc() reproduces Shrout & Fleiss (1979) Table 2 point estimates", {
  res <- icc(sf_data, type = "all")
  # S&F report values to 2 decimal places.
  v <- function(t) round(res$value[res$type == t], 2)
  expect_equal(v("ICC(1,1)"), 0.17)
  expect_equal(v("ICC(2,1)"), 0.29)
  expect_equal(v("ICC(3,1)"), 0.71)
  expect_equal(v("ICC(1,k)"), 0.44)
  expect_equal(v("ICC(2,k)"), 0.62)
  expect_equal(v("ICC(3,k)"), 0.91)
})

test_that("icc() reproduces Shrout & Fleiss (1979) Table 3 confidence intervals", {
  res <- icc(sf_data, type = "all")
  # S&F report intervals to 2 decimal places. Round and compare exactly.
  rounded <- function(r) c(round(r$ci_lower, 2), round(r$ci_upper, 2))

  expect_equal(rounded(res[res$type == "ICC(1,1)", ]), c(-0.13, 0.72))
  expect_equal(rounded(res[res$type == "ICC(2,1)", ]), c( 0.02, 0.76))
  expect_equal(rounded(res[res$type == "ICC(3,1)", ]), c( 0.34, 0.95))
})

# On CRAN the file stops here, after the anchor above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

test_that("icc() returns the documented columns", {
  res <- icc(sf_data, type = "ICC(2,1)")
  expect_named(res, c("type", "value", "ci_lower", "ci_upper",
                      "F_value", "df_1", "df_2", "p_value"))
  expect_equal(res$df_1, 5)
  expect_equal(res$df_2, 15)
})

test_that("icc() type aliases (1, 2, 3, 1k, 2k, 3k) all resolve correctly", {
  v_alias <- icc(sf_data, type = c("1", "2", "3", "1k", "2k", "3k"))$value
  v_full  <- icc(sf_data, type = "all")$value
  # Same set of values, possibly in a different order, but order matches our spec.
  expect_equal(v_alias, v_full)
})

test_that("icc() average-of-k versions are functions of single-rater versions", {
  res <- icc(sf_data, type = "all")
  k <- ncol(sf_data)
  spear_brown <- function(rho, k) (k * rho) / (1 + (k - 1) * rho)
  for (pair in list(c("ICC(1,1)", "ICC(1,k)"),
                    c("ICC(2,1)", "ICC(2,k)"),
                    c("ICC(3,1)", "ICC(3,k)"))) {
    single  <- res$value[res$type == pair[1]]
    average <- res$value[res$type == pair[2]]
    expect_equal(average, spear_brown(single, k), tolerance = 1e-12)
  }
})

test_that("icc() rejects bad inputs", {
  expect_error(icc(sf_data, type = "ICC(7,1)"),  "Unrecognized")
  expect_error(icc(matrix(1:6, nrow = 6)),       "at least 2 rows.*2 columns")
  expect_error(icc(matrix(1:4, nrow = 1)),       "at least 2 rows.*2 columns")
  expect_error(icc(sf_data, conf_level = 1.1),   "conf_level")
  expect_error(icc("not a matrix"),              "matrix")

  bad <- sf_data; bad[1, 1] <- NA
  expect_error(icc(bad), "Missing values")
})

test_that("icc() perfect rater agreement gives ICC(3,1) = 1", {
  perfect <- matrix(rep(1:5, each = 4), nrow = 5, byrow = TRUE)
  # All raters give identical scores -> consistency ICC = 1.
  res <- icc(perfect, type = "ICC(3,1)")
  expect_equal(res$value, 1, tolerance = 1e-12)
})

test_that("icc() with type as a vector returns rows in the requested order", {
  res <- icc(sf_data, type = c("ICC(3,1)", "ICC(1,1)"))
  expect_equal(res$type, c("ICC(3,1)", "ICC(1,1)"))
})

test_that("icc() refuses SPSS user-defined missing values as it refuses NA", {
  # Through 1.0.1 a rating declared missing in SPSS (a haven_labelled_spss
  # column) was analyzed as a rating of 9: ICC(2,1) 0.318 here, against
  # 0.668 on the complete rows.
  set.seed(1)
  n <- 60
  t <- rnorm(n)
  r <- data.frame(a = round(3 + t + rnorm(n, sd = .5)),
                  b = round(3 + t + rnorm(n, sd = .5)),
                  c = round(3 + t + rnorm(n, sd = .5)))
  r$c[1:6] <- NA
  s <- r
  s$c <- structure(as.double(ifelse(is.na(r$c), 9, r$c)),
                   labels = c(Refused = 9), na_values = 9,
                   class = c("haven_labelled_spss", "haven_labelled",
                             "vctrs_vctr", "double"))
  expect_message(
    expect_error(icc(s), "Missing values are not supported by icc()",
                 fixed = TRUE),
    "6 in c", fixed = TRUE)
  expect_equal(icc(r[complete.cases(r), ])$value[1], 0.6684228,
               tolerance = 1e-6)
})
