test_that("cv_z() returns the expected tidy data frame for two-sided alpha", {
  result <- cv_z(alpha_level = 0.05)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "area_less", "area_greater"))
  expect_equal(result$term, c("lower_cv", "upper_cv"))
})

test_that("cv_z() returns a dmar_tbl", {
  expect_s3_class(cv_z(alpha_level = .05), "dmar_tbl")
  expect_s3_class(cv_z(alpha_lower = 0, alpha_upper = .05, verbose = FALSE), "dmar_tbl")
})

test_that("cv_z() at alpha = 0.05 matches qnorm(0.025) and qnorm(0.975)", {
  result <- cv_z(alpha_level = 0.05)
  expect_equal(result$value[result$term == "lower_cv"], qnorm(0.025), tolerance = 1e-12)
  expect_equal(result$value[result$term == "upper_cv"], qnorm(0.975), tolerance = 1e-12)
})

test_that("cv_z() with alpha = 0.01 widens the critical values", {
  r05 <- cv_z(alpha_level = 0.05)
  r01 <- cv_z(alpha_level = 0.01)
  expect_lt(r01$value[r01$term == "lower_cv"], r05$value[r05$term == "lower_cv"])
  expect_gt(r01$value[r01$term == "upper_cv"], r05$value[r05$term == "upper_cv"])
})
