test_that("ss_aipe_sc() returns a tidy single-row data frame with the planned sample size", {
  res <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "necessary_n_per_group")
  expect_true(is.numeric(res$value) && res$value > 0)
})

test_that("ss_aipe_sc() default conf_level = .95 produces a sensible per-group N", {
  res <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4)
  # Documented value: a moderately large per-group sample is needed for a
  # standardized contrast of .6 to be estimated to within .4 at 95%.
  expect_gt(res$value, 50)
  expect_lt(res$value, 500)
})

test_that("ss_aipe_sc() assurance path no longer errors on inner conf_limits_nct calls", {
  # Regression test: the documented example used to fail with
  # "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'"
  # because the inner conf_limits_nct() calls did not pass conf_level = NULL.
  res <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0),
                    width = .4, assurance = .90)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "necessary_n_per_group")
  # Assurance must (weakly) increase the planned N versus the expected-width
  # plan, because the constraint is more demanding.
  res_noassur <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0),
                            width = .4)
  expect_gte(res$value, res_noassur$value)
})

test_that("ss_aipe_sc() rejects illegal contrast weights", {
  expect_error(
    ss_aipe_sc(psi_standardized = .5, c_weights = c(1, 1), width = .3),
    "sum of the coefficients"
  )
  expect_error(
    ss_aipe_sc(psi_standardized = .5, c_weights = c(2, 1, -3), width = .3),
    "fractions"
  )
})

test_that("ss_aipe_sc() accepts explicit alpha_lower/alpha_upper for symmetric CI", {
  # Explicit symmetric alphas must reproduce the default conf_level = .95 plan.
  ref <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4)
  alt <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
                    conf_level = NULL,
                    alpha_lower = .025, alpha_upper = .025)
  expect_equal(alt$value, ref$value)
})

test_that("ss_aipe_sc() accepts asymmetric alpha_lower/alpha_upper", {
  # Asymmetric split with the same total alpha (here .05) generally requires
  # a different N. Just check the function runs and returns a sensible value.
  res <- ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
                    conf_level = NULL,
                    alpha_lower = .01, alpha_upper = .04)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "necessary_n_per_group")
  expect_gt(res$value, 0)
})

test_that("ss_aipe_sc() rejects mixing conf_level and alphas", {
  expect_error(
    ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = .95, alpha_lower = .025, alpha_upper = .025),
    "cannot mix them"
  )
})

test_that("ss_aipe_sc() rejects single-tail alpha (need both)", {
  expect_error(
    ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = NULL, alpha_lower = .025),
    "Supply both"
  )
})

test_that("ss_aipe_sc() rejects out-of-range alpha bounds", {
  expect_error(
    ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = NULL, alpha_lower = -.01, alpha_upper = .05),
    "strictly positive"
  )
  expect_error(
    ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = NULL, alpha_lower = .5, alpha_upper = .5),
    "sum to less than 1"
  )
})

test_that("ss_aipe_sc() rejects out-of-range assurance values", {
  expect_error(
    ss_aipe_sc(psi_standardized = .5, c_weights = c(.5, .5, -.5, -.5), width = .3,
               assurance = 1.1),
    "assurance"
  )
  expect_error(
    ss_aipe_sc(psi_standardized = .5, c_weights = c(.5, .5, -.5, -.5), width = .3,
               assurance = .3),
    "larger than 0\\.5"
  )
})
