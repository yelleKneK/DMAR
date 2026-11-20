test_that("cv_tukey_hsd() returns a tidy data.frame with term='upper_cv'", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3)
  expect_s3_class(result, "data.frame")
  expect_equal(result$term, "upper_cv")
  expect_true(all(c("term", "value", "area_less", "area_greater") %in% names(result)))
})

test_that("cv_tukey_hsd() value matches qtukey()/sqrt(2)", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3)
  expected <- qtukey(p = 1 - .05, nmeans = 3, df = 27) / sqrt(2)
  expect_equal(result$value, expected)
})

test_that("cv_tukey_hsd() verbose=FALSE drops the area columns", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3, verbose = FALSE)
  expect_equal(names(result), c("term", "value"))
})

test_that("cv_tukey_hsd() area_less + area_greater sum to 1", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3)
  expect_equal(result$area_less + result$area_greater, 1)
})

test_that("cv_tukey_hsd() area_greater equals alpha", {
  result <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3)
  expect_equal(result$area_greater, .05)
})

test_that("cv_tukey_hsd() larger alpha produces smaller critical value", {
  strict <- cv_tukey_hsd(alpha_level = .01, df = 27, n_groups = 3)
  lenient <- cv_tukey_hsd(alpha_level = .10, df = 27, n_groups = 3)
  expect_gt(strict$value, lenient$value)
})

test_that("cv_tukey_hsd() more groups produces larger critical value", {
  three <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 3)
  six   <- cv_tukey_hsd(alpha_level = .05, df = 27, n_groups = 6)
  expect_gt(six$value, three$value)
})

test_that("cv_tukey_hsd() errors on missing df, groups, alpha", {
  expect_error(cv_tukey_hsd(alpha_level = .05, n_groups = 3), "degrees of freedom")
  expect_error(cv_tukey_hsd(alpha_level = .05, df = 27), "number of groups")
  expect_error(cv_tukey_hsd(df = 27, n_groups = 3), "'alpha_level'")
})

test_that("cv_tukey_hsd() errors on invalid argument values", {
  expect_error(cv_tukey_hsd(alpha_level = 0,    df = 27, n_groups = 3), "strictly between 0 and 1")
  expect_error(cv_tukey_hsd(alpha_level = 1,    df = 27, n_groups = 3), "strictly between 0 and 1")
  expect_error(cv_tukey_hsd(alpha_level = .05,  df = 0,  n_groups = 3), "positive value")
  expect_error(cv_tukey_hsd(alpha_level = .05,  df = 27, n_groups = 1), "two or more groups")
})

test_that("cv_tukey_hsd() at 2 and 3 error df is the exact critical value, not qtukey()'s", {
  # Oracle: an independent quadrature of the studentized range (the 1.1.0
  # review's F3 evidence), q_.99;20,2 = 37.943462287637 and
  # q_.99;20,3 = 19.764809074501; divided by sqrt(2), 26.830 and 13.976.
  # stats::qtukey() gives 31.485 and 14.105 on this scale (17 and 0.9
  # percent too large).
  two <- cv_tukey_hsd(alpha_level = .01, df = 2, n_groups = 20)
  three <- cv_tukey_hsd(alpha_level = .01, df = 3, n_groups = 20)
  expect_equal(two$value, 37.943462287637 / sqrt(2), tolerance = 1e-8)
  expect_equal(three$value, 19.764809074501 / sqrt(2), tolerance = 1e-8)
  expect_equal(two$area_greater, 0.01, tolerance = 1e-8)
  expect_equal(three$area_less, 0.99, tolerance = 1e-8)
  expect_gt(qtukey(.99, 20, 2) / sqrt(2) - two$value, 4.6)
})

test_that("cv_tukey_hsd() works below 2 error df, where qtukey() returns NaN", {
  # Two groups: q = sqrt(2) t, so the critical value is t_{1 - alpha/2, df}.
  out <- cv_tukey_hsd(alpha_level = .05, df = 1.5, n_groups = 2)
  expect_equal(out$value, qt(0.975, 1.5), tolerance = 1e-12)
  expect_equal(out$area_greater, 0.05, tolerance = 1e-10)
  skip_on_cran()  # a spline of the range distribution for four means
  # Four groups: the independent quadrature gives q_.95;4,1.5 = 14.40883847907.
  expect_true(is.nan(suppressWarnings(qtukey(0.95, 4, 1.5))))
  expect_equal(cv_tukey_hsd(.05, df = 1.5, n_groups = 4)$value,
               14.40883847907 / sqrt(2), tolerance = 1e-8)
})

test_that("cv_tukey_hsd() at alpha .0001 is the exact critical value, not qtukey()'s", {
  skip_on_cran()  # a spline of the range distribution for ten means
  # qtukey(.9999, 10, 11) is 12.008462418, 1.17e-4 below the exact quantile
  # 12.008579776018 (the 1.1.0 review's two independent quadratures agree to
  # 1e-11); divided by sqrt(2), 8.4912652 against 8.4913482.
  out <- cv_tukey_hsd(alpha_level = 1e-4, df = 11, n_groups = 10)
  expect_lt(abs(out$value - 8.4913482), 1e-6)
  expect_equal(out$value, 12.0085797760180 / sqrt(2), tolerance = 1e-9)
  expect_equal(out$area_greater, 1e-4, tolerance = 1e-3)
  expect_gt(out$value - qtukey(0.9999, 10, 11) / sqrt(2), 8e-5)
})

test_that("cv_tukey_hsd() with 1000 groups is the exact critical value", {
  skip_on_cran()  # a spline of the range distribution for 1000 means
  # qtukey(.999, 1000, 41) is 10.299349568, 1.45e-4 below the exact quantile
  # 10.2994946166086 (both independent quadratures); divided by sqrt(2),
  # 7.2827399 against 7.2828425.
  out <- cv_tukey_hsd(alpha_level = .001, df = 41, n_groups = 1000)
  expect_equal(out$value, 10.2994946166086 / sqrt(2), tolerance = 1e-9)
  expect_equal(out$area_greater, 0.001, tolerance = 1e-6)
  expect_gt(out$value - qtukey(0.999, 1000, 41) / sqrt(2), 1e-4)
})
