test_that("cv_tukey_hsd() returns a data frame with term='upper_cv'", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
  expect_s3_class(result, "data.frame")
  expect_equal(result$term, "upper_cv")
})

test_that("cv_tukey_hsd() verbose output includes area_less and area_greater", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3, verbose = TRUE)
  expect_named(result, c("term", "value", "area_less", "area_greater"))
  # area_greater should equal alpha for exact critical value
  expect_equal(result$area_greater, 0.05, tolerance = 1e-10)
})

test_that("cv_tukey_hsd() non-verbose output has only term and value", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("cv_tukey_hsd() matches manual qtukey calculation", {
  alpha <- 0.05
  df <- 60
  groups <- 4
  expected <- qtukey(1 - alpha, nmeans = groups, df = df) / sqrt(2)
  result <- cv_tukey_hsd(alpha_level = alpha, df = df, groups = groups)
  expect_equal(result$value, expected)
})

test_that("cv_tukey_hsd() critical value increases with more groups", {
  cv_3 <- cv_tukey_hsd(alpha_level = .05, df = 30, groups = 3)$value
  cv_5 <- cv_tukey_hsd(alpha_level = .05, df = 30, groups = 5)$value
  expect_gt(cv_5, cv_3)
})

test_that("cv_tukey_hsd() critical value decreases with more df", {
  cv_small_df <- cv_tukey_hsd(alpha_level = .05, df = 10, groups = 3)$value
  cv_large_df <- cv_tukey_hsd(alpha_level = .05, df = 100, groups = 3)$value
  expect_gt(cv_small_df, cv_large_df)
})

test_that("cv_tukey_hsd() errors on bad inputs", {
  expect_error(cv_tukey_hsd(df = 10, groups = 3), "alpha")
  expect_error(cv_tukey_hsd(alpha_level = .05, groups = 3), "degrees of freedom")
  expect_error(cv_tukey_hsd(alpha_level = .05, df = 10), "groups")
  expect_error(cv_tukey_hsd(alpha_level = .05, df = 10, groups = 1), "two or more")
  expect_error(cv_tukey_hsd(alpha_level = 0, df = 10, groups = 3), "greater than zero")
  expect_error(cv_tukey_hsd(alpha_level = .05, df = -1, groups = 3), "positive")
})
