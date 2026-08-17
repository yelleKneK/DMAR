test_that("convert_r_Z and convert_Z_r are exact inverses", {
  rs <- c(-0.9, -0.5, -0.2, 0, 0.2, 0.5, 0.9)
  for (r in rs) {
    Z <- convert_r_Z(r)$value
    r_back <- convert_Z_r(Z)$value
    expect_equal(r_back, r, tolerance = 1e-12)
  }
})


test_that("convert_Z_r is numerically stable at large |Z|", {
  big <- 400
  out <- convert_Z_r(big)
  expect_true(is.finite(out$value))
  expect_equal(out$value, 1, tolerance = 1e-12)
  out_neg <- convert_Z_r(-big)
  expect_true(is.finite(out_neg$value))
  expect_equal(out_neg$value, -1, tolerance = 1e-12)
})


test_that("convert_Z_r matches the Hays (1994) reference value", {
  # Hays (1994, pp. 649--650): Fisher Z of 0.3654438 -> r of approximately 0.35.
  out <- convert_Z_r(0.3654438)
  expect_equal(out$value, 0.35, tolerance = 1e-4)
})


test_that("convert_r_Z rejects |r| >= 1 with an informative error", {
  expect_error(convert_r_Z(1), "must satisfy")
  expect_error(convert_r_Z(-1), "must satisfy")
  expect_error(convert_r_Z(1.5), "must satisfy")
})


test_that("convert_r_Z rejects non-numeric input", {
  expect_error(convert_r_Z("0.5"), "single numeric value")
})


test_that("convert_Z_r returns a tidy term/value data.frame", {
  out <- convert_Z_r(0.5)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "value"))
  expect_equal(out$term, "r_from_Z")
  expect_true(is.numeric(out$value))
})


test_that("convert_r_Z returns a tidy term/value data.frame", {
  out <- convert_r_Z(0.5)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "value"))
  expect_equal(out$term, "Z_from_r")
  expect_true(is.numeric(out$value))
})


test_that("convert_r_Z / convert_Z_r aliases match the canonical functions", {
  expect_identical(convert_r_Z, convert_r_Z)
  expect_identical(convert_Z_r, convert_Z_r)
  expect_equal(convert_r_Z(0.5)$value, convert_r_Z(0.5)$value, tolerance = 1e-12)
  expect_equal(convert_Z_r(0.5)$value, convert_Z_r(0.5)$value, tolerance = 1e-12)
})


test_that("convert_R2_f and convert_f_R2 are exact inverses", {
  R2s <- c(0.05, 0.2, 0.5, 0.8)
  N <- 100; p <- 3
  df_1 <- p
  df_2 <- N - p - 1
  for (R2 in R2s) {
    F_val <- convert_R2_f(R2 = R2, df_1 = df_1, df_2 = df_2,
                          p = p, N = N)$value
    R2_back <- convert_f_R2(F_value = F_val, df_1 = df_1,
                            df_2 = df_2)$value
    expect_equal(R2_back, R2, tolerance = 1e-10)
  }
})


test_that("convert_lambda_R2 and convert_R2_lambda are exact inverses", {
  R2s <- c(0.1, 0.3, 0.6)
  for (R2 in R2s) {
    lam <- convert_R2_lambda(R2 = R2, N = 100)$value
    R2_back <- convert_lambda_R2(lambda = lam, N = 100)$value
    expect_equal(R2_back, R2, tolerance = 1e-10)
  }
})


test_that("convert_delta_lambda and convert_lambda_delta are exact inverses", {
  d <- 0.4
  lam <- convert_delta_lambda(delta = d, n_1 = 30, n_2 = 30)$value
  d_back <- convert_lambda_delta(lambda = lam, n_1 = 30, n_2 = 30)$value
  expect_equal(d_back, d, tolerance = 1e-12)
})

test_that("every scalar convert function refuses vector and missing input", {
  # The 2026-08-17 family ruling: no recycling, ever. Each function takes
  # a single value; vector users get pointed at the underlying arithmetic.
  expect_error(convert_R2_f(R2 = c(.2, .5), df_1 = 3, df_2 = 96),
               "single number")
  expect_error(convert_f_R2(F_value = c(2, 4), df_1 = 3, df_2 = 96),
               "single number")
  expect_error(convert_lambda_R2(lambda = c(5, 10), N = 100),
               "single number")
  expect_error(convert_R2_lambda(R2 = c(.2, .5), N = 100),
               "single number")
  expect_error(convert_R2_lambda(R2 = NA_real_, N = 100), "single number")
  expect_error(convert_delta_lambda(delta = c(.2, .5), n_1 = 20, n_2 = 20),
               "single number")
  expect_error(convert_lambda_delta(lambda = c(1, 2), n_1 = 20, n_2 = 20),
               "single number")
  expect_error(convert_z_normal(z = c(1, 2)), "single number")
})

test_that("the convert_R2 domain guards stop out-of-range input", {
  expect_error(convert_R2_f(R2 = 1, df_1 = 3, df_2 = 96), "less than 1")
  expect_error(convert_R2_lambda(R2 = 1.2, N = 100), "less than 1")
  expect_error(convert_f_R2(F_value = -1, df_1 = 3, df_2 = 96),
               "at least 0")
  expect_error(convert_delta_lambda(delta = 0.5, n_1 = 0, n_2 = 20),
               "greater than 0")
})
