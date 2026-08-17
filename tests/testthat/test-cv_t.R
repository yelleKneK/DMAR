test_that("cv_t() returns the expected tidy data frame for two-sided alpha", {
  result <- cv_t(alpha_level = 0.05, df = 30)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "area_less", "area_greater"))
  expect_equal(result$term, c("lower_cv", "upper_cv"))
})

test_that("cv_t() returns a dmar_tbl", {
  expect_s3_class(cv_t(alpha_level = .05, df = 13), "dmar_tbl")
  expect_s3_class(cv_t(alpha_lower = 0, alpha_upper = .05, df = 13, verbose = FALSE), "dmar_tbl")
})

test_that("cv_t() at alpha = 0.05, df = 30 matches qt(0.025, 30) and qt(0.975, 30)", {
  result <- cv_t(alpha_level = 0.05, df = 30)
  expect_equal(result$value[result$term == "lower_cv"], qt(0.025, df = 30), tolerance = 1e-12)
  expect_equal(result$value[result$term == "upper_cv"], qt(0.975, df = 30), tolerance = 1e-12)
})

test_that("cv_t() with ncp != 0 returns the noncentral t critical values", {
  result <- cv_t(alpha_level = 0.05, df = 30, ncp = 2)
  expect_equal(result$value[result$term == "lower_cv"], qt(0.025, df = 30, ncp = 2), tolerance = 1e-10)
  expect_equal(result$value[result$term == "upper_cv"], qt(0.975, df = 30, ncp = 2), tolerance = 1e-10)
})

test_that("the lower-tail alternative accepts 'less', the counterpart of 'greater'", {
  # The upper-tail accept list carried "greater" while the lower-tail list
  # carried "lesser" but not "less", so the natural counterpart of the
  # spelling that worked did not. It did not error cleanly either: the call
  # fell through both branches and died on an unset condition with
  # "missing value where TRUE/FALSE needed".
  expect_equal(cv_t(alpha_level = .05, df = 13, alternative = "less"),
               cv_t(alpha_level = .05, df = 13, alternative = "lesser"))
  expect_equal(cv_f(alpha_level = .05, df_numerator = 2, df_denominator = 20,
                    alternative = "less"),
               cv_f(alpha_level = .05, df_numerator = 2, df_denominator = 20,
                    alternative = "lesser"))
  expect_equal(cv_chisq(alpha_level = .05, df = 3, alternative = "less"),
               cv_chisq(alpha_level = .05, df = 3, alternative = "lesser"))
  # And "less" really is routed to the lower tail: it spends the whole of
  # alpha there, so its lower critical value sits closer to zero than the
  # two-sided one, and nothing is cut off above.
  one <- cv_t(alpha_level = .05, df = 13, alternative = "less")
  two <- cv_t(alpha_level = .05, df = 13, alternative = "not_equal")
  expect_gt(one$value[one$term == "lower_cv"],
            two$value[two$term == "lower_cv"])
  expect_identical(one$value[one$term == "upper_cv"], Inf)
})
