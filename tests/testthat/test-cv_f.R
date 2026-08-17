test_that("cv_f() reproduces Appendix Table A.2 of MDK (2027)", {
  # Denominator df = 4, alpha = .05, numerator df = 1, 2, ..., 10.
  printed <- c(7.71, 6.94, 6.59, 6.39, 6.26, 6.16, 6.09, 6.04, 6.00, 5.96)
  got <- vapply(1:10, function(k)
    cv_f(alpha_level = .05, df_numerator = k, df_denominator = 4,
         verbose = FALSE)$value[2], numeric(1))
  expect_equal(round(got, 2), printed)
})

test_that("cv_f() is the upper-tail F quantile by default", {
  out <- cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20)
  expect_equal(out$value[2], qf(.95, df1 = 3, df2 = 20), tolerance = 1e-10)
  # A zero-area lower tail sits at the boundary of the support.
  expect_equal(out$value[1], 0)
  expect_equal(out$area_greater[2], .05, tolerance = 1e-10)
})

test_that("cv_f() splits alpha between the tails when asked", {
  out <- cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
              alternative = "not_equal")
  expect_equal(out$value[1], qf(.025, df1 = 3, df2 = 20), tolerance = 1e-10)
  expect_equal(out$value[2], qf(.975, df1 = 3, df2 = 20), tolerance = 1e-10)
  expect_gt(out$value[1], 0)
})

test_that("cv_f() honors alpha_lower and alpha_upper directly", {
  out <- cv_f(alpha_lower = .01, alpha_upper = .04, df_numerator = 3,
              df_denominator = 20)
  expect_equal(out$value[1], qf(.01, df1 = 3, df2 = 20), tolerance = 1e-10)
  expect_equal(out$value[2], qf(.96, df1 = 3, df2 = 20), tolerance = 1e-10)
})

test_that("cv_f() accepts a noncentral parameter", {
  out <- cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20, ncp = 5,
              verbose = FALSE)
  expect_equal(out$value[2], qf(.95, df1 = 3, df2 = 20, ncp = 5),
               tolerance = 1e-10)
  # A noncentral F is shifted up relative to the central one.
  expect_gt(out$value[2], qf(.95, df1 = 3, df2 = 20))
})

test_that("cv_f() returns the documented structure", {
  out <- cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "value", "area_less", "area_greater"))
  expect_equal(out$term, c("lower_cv", "upper_cv"))
  simple <- cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
                 verbose = FALSE)
  expect_named(simple, c("term", "value"))
})

test_that("cv_f() errors on bad input", {
  expect_error(cv_f(alpha_level = .05), "df_numerator")
  expect_error(cv_f(alpha_level = .05, df_numerator = 3), "df_denominator")
  expect_error(cv_f(df_numerator = 3, df_denominator = 20), "alpha")
  expect_error(cv_f(alpha_level = .05, df_numerator = 0, df_denominator = 20),
               "df_numerator")
  expect_error(cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 0),
               "df_denominator")
  expect_error(cv_f(alpha_level = 1, df_numerator = 3, df_denominator = 20), "alpha")
  expect_error(cv_f(alpha_level = .05, alpha_lower = .01, df_numerator = 3,
                    df_denominator = 20), "only use one approach")
})

test_that("cv_f() handles an infinite numerator df (Table A.2's last column)", {
  # Supplying ncp = 0 explicitly would send qf() down its noncentral algorithm,
  # which returns NaN for an infinite numerator df. The central algorithm gives
  # the value the table prints, F_{.05; Inf, 10} = 2.54.
  v <- cv_f(alpha_level = .05, df_numerator = Inf, df_denominator = 10,
            verbose = FALSE)$value[2]
  expect_false(is.nan(v))
  expect_equal(v, qf(.95, Inf, 10), tolerance = 1e-10)
  expect_equal(round(v, 2), 2.54)
  # The infinite-numerator, infinite-denominator corner is 1.00.
  expect_equal(round(cv_f(alpha_level = .05, df_numerator = Inf, df_denominator = Inf,
                          verbose = FALSE)$value[2], 2), 1.00)
})

test_that("cv_f() still honors a nonzero noncentral parameter", {
  v <- cv_f(alpha_level = .05, df_numerator = 2, df_denominator = 20, ncp = 5,
            verbose = FALSE)$value[2]
  expect_equal(v, qf(.95, 2, 20, ncp = 5), tolerance = 1e-10)
})

test_that("cv_f() rejects a negative or non-finite ncp instead of returning NaN", {
  # A negative ncp routes qf()/pf() into the noncentral algorithm out of domain
  # and returns a polished NaN; require a finite, non-negative value.
  expect_error(cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
                    ncp = -1), "ncp")
  expect_error(cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
                    ncp = NA_real_), "ncp")
  expect_error(cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
                    ncp = Inf), "ncp")
})
