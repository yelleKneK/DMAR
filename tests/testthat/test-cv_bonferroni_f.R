test_that("cv_bonferroni_f() reproduces Appendix Table A.3 of MDK (2027)", {
  # The table is for 1 numerator degree of freedom and a family-wise alpha of
  # .05. Rows are the denominator (error) df; columns are C = 1, 2, ..., 10.
  printed <- rbind(
    c(161.45, 647.79, 1458.36, 2593.16, 4052.18, 5835.43, 7942.91, 10374.62,
      13130.56, 16210.72),
    c( 18.51,  38.51,   58.50,   78.50,   98.50,  118.50,  138.50,   158.50,
        178.50,   198.50),
    c( 10.13,  17.44,   23.59,   29.07,   34.12,   38.83,   43.29,    47.54,
         51.62,    55.55),
    c(  7.71,  12.22,   15.69,   18.62,   21.20,   23.53,   25.68,    27.68,
         29.56,    31.33)
  )
  for (df in 1:4) {
    got <- vapply(1:10, function(C)
      cv_bonferroni_f(alpha_level = .05, df_denominator = df, n_comparisons = C,
                      verbose = FALSE)$value, numeric(1))
    expect_equal(round(got, 2), printed[df, ])
  }
})

test_that("cv_bonferroni_f() is the F critical value read at alpha / C", {
  b <- cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5,
                       verbose = FALSE)$value
  f <- cv_f(alpha_level = .05 / 5, df_numerator = 1, df_denominator = 20,
            verbose = FALSE)$value[2]
  expect_equal(b, f, tolerance = 1e-10)
})

test_that("cv_bonferroni_f() reports the per-comparison rate in area_greater", {
  out <- cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5)
  expect_equal(out$area_greater, .05 / 5, tolerance = 1e-10)
  expect_equal(out$term, "upper_cv")
})

test_that("cv_bonferroni_f() with one comparison leaves alpha alone", {
  expect_equal(
    cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 1,
                    verbose = FALSE)$value,
    qf(.95, df1 = 1, df2 = 20), tolerance = 1e-10)
})

test_that("cv_bonferroni_f() grows with the size of the family", {
  v <- vapply(1:10, function(C)
    cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = C,
                    verbose = FALSE)$value, numeric(1))
  expect_true(all(diff(v) > 0))
})

test_that("cv_bonferroni_f() allows a larger numerator df", {
  expect_equal(
    cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 4,
                    df_numerator = 3, verbose = FALSE)$value,
    qf(1 - .05 / 4, df1 = 3, df2 = 20), tolerance = 1e-10)
})

test_that("cv_bonferroni_f() returns the documented structure", {
  out <- cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "value", "area_less", "area_greater"))
  expect_named(cv_bonferroni_f(alpha_level = .05, df_denominator = 20,
                               n_comparisons = 5, verbose = FALSE),
               c("term", "value"))
})

test_that("cv_bonferroni_f() errors on bad input", {
  expect_error(cv_bonferroni_f(alpha_level = .05), "df_denominator")
  expect_error(cv_bonferroni_f(alpha_level = .05, df_denominator = 20),
               "n_comparisons")
  expect_error(cv_bonferroni_f(alpha_level = .05, df_denominator = 0,
                               n_comparisons = 5), "df_denominator")
  expect_error(cv_bonferroni_f(alpha_level = .05, df_denominator = 20,
                               n_comparisons = 0), "n_comparisons")
  expect_error(cv_bonferroni_f(alpha_level = .05, df_denominator = 20,
                               n_comparisons = 2.5), "n_comparisons")
  expect_error(cv_bonferroni_f(alpha_level = 1, df_denominator = 20,
                               n_comparisons = 5), "alpha")
})
