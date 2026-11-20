test_that("cv_scheffe() matches sqrt((k-1) * qf(1-alpha, k-1, df))", {
  res <- cv_scheffe(alpha_level = .05, df_numerator = 3, df_denominator = 36)
  expected <- sqrt(3 * qf(.95, 3, 36))
  expect_equal(res$value, expected, tolerance = 1e-12)
  expect_equal(res$term, "upper_cv")
})

test_that("cv_scheffe() verbose=FALSE drops area columns", {
  res <- cv_scheffe(alpha_level = .05, df_numerator = 3, df_denominator = 36,
                    verbose = FALSE)
  expect_named(res, c("term", "value"))
})

test_that("cv_scheffe() is monotonic in alpha (smaller alpha -> larger crit)", {
  strict <- cv_scheffe(alpha_level = .01, df_numerator = 3, df_denominator = 36)$value
  lenient <- cv_scheffe(alpha_level = .10, df_numerator = 3, df_denominator = 36)$value
  expect_gt(strict, lenient)
})

test_that("cv_scheffe() is monotonic in df_numerator (more contrasts -> larger crit)", {
  small <- cv_scheffe(alpha_level = .05, df_numerator = 2, df_denominator = 100)$value
  large <- cv_scheffe(alpha_level = .05, df_numerator = 5, df_denominator = 100)$value
  expect_gt(large, small)
})

test_that("cv_scheffe() errors on bad inputs", {
  expect_error(cv_scheffe(),                                                "alpha")
  expect_error(cv_scheffe(alpha_level = 0,    df_numerator = 3, df_denominator = 36), "alpha")
  expect_error(cv_scheffe(alpha_level = 1.1,  df_numerator = 3, df_denominator = 36), "alpha")
  expect_error(cv_scheffe(alpha_level = .05,  df_numerator = 0, df_denominator = 36), "df_numerator")
  expect_error(cv_scheffe(alpha_level = .05,  df_numerator = 3, df_denominator = 0),  "df_denominator")
})
