test_that("cv_tukey_hsd() returns a tidy data.frame with term='upper_cv'", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
  expect_s3_class(result, "data.frame")
  expect_equal(result$term, "upper_cv")
  expect_true(all(c("term", "value", "area_less", "area_greater") %in% names(result)))
})

test_that("cv_tukey_hsd() value matches qtukey()/sqrt(2)", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
  expected <- qtukey(p = 1 - .05, nmeans = 3, df = 27) / sqrt(2)
  expect_equal(result$value, expected)
})

test_that("cv_tukey_hsd() verbose=FALSE drops the area columns", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3, verbose = FALSE)
  expect_equal(names(result), c("term", "value"))
})

test_that("cv_tukey_hsd() area_less + area_greater sum to 1", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
  expect_equal(result$area_less + result$area_greater, 1)
})

test_that("cv_tukey_hsd() area_greater equals alpha", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
  expect_equal(result$area_greater, .05)
})

test_that("cv_tukey_hsd() larger alpha produces smaller critical value", {
  strict <- cv_tukey_hsd(alpha_level = .01, df = 27, groups = 3)
  lenient <- cv_tukey_hsd(alpha_level = .10, df = 27, groups = 3)
  expect_gt(strict$value, lenient$value)
})

test_that("cv_tukey_hsd() more groups produces larger critical value", {
  three <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
  six   <- cv_tukey_hsd(alpha_level = .05, df = 27, groups = 6)
  expect_gt(six$value, three$value)
})

test_that("cv_tukey_hsd() errors on missing df, groups, alpha", {
  expect_error(cv_tukey_hsd(alpha_level = .05, groups = 3), "degrees of freedom")
  expect_error(cv_tukey_hsd(alpha_level = .05, df = 27), "number of groups")
  expect_error(cv_tukey_hsd(df = 27, groups = 3), "'alpha_level'")
})

test_that("cv_tukey_hsd() errors on invalid argument values", {
  expect_error(cv_tukey_hsd(alpha_level = 0,    df = 27, groups = 3), "greater than zero")
  expect_error(cv_tukey_hsd(alpha_level = 1,    df = 27, groups = 3), "less than 1")
  expect_error(cv_tukey_hsd(alpha_level = .05,  df = 0,  groups = 3), "positive value")
  expect_error(cv_tukey_hsd(alpha_level = .05,  df = 27, groups = 1), "two or more groups")
})
