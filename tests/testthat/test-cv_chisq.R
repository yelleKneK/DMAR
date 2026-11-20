test_that("cv_chisq() reproduces Appendix Table A.9 of MDK (2027)", {
  # Rows are df = 1 through 5; columns are the upper-tail areas the table
  # gives: .10, .05, .025, .01, .005, .001.
  alphas <- c(.10, .05, .025, .01, .005, .001)
  printed <- rbind(
    c(2.71,  3.84,  5.02,  6.63,  7.88, 10.83),
    c(4.61,  5.99,  7.38,  9.21, 10.60, 13.82),
    c(6.25,  7.81,  9.35, 11.35, 12.84, 16.27),
    c(7.78,  9.49, 11.14, 13.28, 14.86, 18.47),
    c(9.24, 11.07, 12.83, 15.09, 16.75, 20.52)
  )
  # Ten of the table's 180 entries have an exact value within 0.0005 of the
  # point where the second decimal turns over, so which way a two-decimal table
  # rounds them is arbitrary. The entry at df = 3, alpha = .01 is one of them:
  # its exact value is 11.34487 and the table prints 11.35. Agreement is
  # therefore asserted to the table's printed precision.
  for (df in 1:5) {
    got <- vapply(alphas, function(a)
      cv_chisq(alpha_level = a, df = df, verbose = FALSE)$value[2], numeric(1))
    expect_true(all(abs(got - printed[df, ]) < 0.006))
  }
})

test_that("cv_chisq() is exact, whatever a two-decimal table rounds to", {
  # The boundary entries above are a property of the tabulation, not of the
  # computation; the quantile itself is exact to machine precision.
  expect_equal(cv_chisq(alpha_level = .01, df = 3, verbose = FALSE)$value[2],
               qchisq(.99, df = 3), tolerance = 1e-12)
  expect_equal(cv_chisq(alpha_level = .01, df = 3, verbose = FALSE)$value[2],
               11.344867, tolerance = 1e-6)
})

test_that("cv_chisq() is the upper-tail chi square quantile by default", {
  out <- cv_chisq(alpha_level = .05, df = 3)
  expect_equal(out$value[2], qchisq(.95, df = 3), tolerance = 1e-10)
  # A zero-area lower tail sits at the boundary of the support.
  expect_equal(out$value[1], 0)
  expect_equal(out$area_greater[2], .05, tolerance = 1e-10)
})

test_that("cv_chisq() splits alpha between the tails when asked", {
  # The pair an interval for a variance uses.
  out <- cv_chisq(alpha_level = .05, df = 10, alternative = "not_equal")
  expect_equal(out$value[1], qchisq(.025, df = 10), tolerance = 1e-10)
  expect_equal(out$value[2], qchisq(.975, df = 10), tolerance = 1e-10)
  expect_gt(out$value[1], 0)
})

test_that("cv_chisq() honors alpha_lower and alpha_upper directly", {
  out <- cv_chisq(alpha_lower = .01, alpha_upper = .04, df = 10)
  expect_equal(out$value[1], qchisq(.01, df = 10), tolerance = 1e-10)
  expect_equal(out$value[2], qchisq(.96, df = 10), tolerance = 1e-10)
})

test_that("cv_chisq() accepts a noncentral parameter", {
  out <- cv_chisq(alpha_level = .05, df = 3, ncp = 5, verbose = FALSE)
  expect_equal(out$value[2], qchisq(.95, df = 3, ncp = 5), tolerance = 1e-10)
  expect_gt(out$value[2], qchisq(.95, df = 3))
})

test_that("cv_chisq() returns the documented structure", {
  out <- cv_chisq(alpha_level = .05, df = 3)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "value", "area_less", "area_greater"))
  expect_equal(out$term, c("lower_cv", "upper_cv"))
  expect_named(cv_chisq(alpha_level = .05, df = 3, verbose = FALSE),
               c("term", "value"))
})

test_that("cv_chisq() errors on bad input", {
  expect_error(cv_chisq(alpha_level = .05), "degrees of freedom")
  expect_error(cv_chisq(df = 3), "alpha")
  expect_error(cv_chisq(alpha_level = .05, df = 0), "positive value")
  expect_error(cv_chisq(alpha_level = 1, df = 3), "alpha")
  expect_error(cv_chisq(alpha_level = .05, alpha_lower = .01, df = 3),
               "only use one approach")
})

test_that("cv_chisq() rejects a negative or non-finite ncp instead of returning NaN", {
  # A negative ncp routes qchisq()/pchisq() into the noncentral algorithm out of
  # domain and returns a polished NaN; require a finite, non-negative value.
  expect_error(cv_chisq(alpha_level = .05, df = 3, ncp = -1), "ncp")
  expect_error(cv_chisq(alpha_level = .05, df = 3, ncp = NA_real_), "ncp")
  expect_error(cv_chisq(alpha_level = .05, df = 3, ncp = Inf), "ncp")
})
