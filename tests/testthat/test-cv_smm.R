test_that("cv_smm() with m = 1 reduces exactly to the two-sided t critical value", {
  res <- cv_smm(alpha_level = .05, df = 36, n_comparisons = 1)
  expect_equal(res$value, qt(.975, df = 36), tolerance = 1e-8)
})

test_that("cv_smm() at df = Inf matches the closed-form normal maximum modulus", {
  # P(max|Z_i| <= c) = (2 Phi(c) - 1)^m = 1 - alpha  ->  c = qnorm((1+(1-a)^(1/m))/2)
  m <- 6
  expect_equal(cv_smm(.05, df = Inf, n_comparisons = m)$value,
               qnorm((1 + 0.95^(1 / m)) / 2), tolerance = 1e-8)
})

test_that("cv_smm() reproduces published SMM values (Appendix Table A.5)", {
  # 'Number of groups' a = 4, 5, 6 -> m = a(a-1)/2 pairwise comparisons,
  # alpha = .05, large df: printed 2.63, 2.80, 2.93.
  got <- vapply(c(6, 10, 15),
                function(m) cv_smm(.05, df = Inf, n_comparisons = m,
                                   verbose = FALSE)$value, numeric(1))
  expect_equal(round(got, 2), c(2.63, 2.80, 2.93))
})

test_that("cv_smm() is deterministic (no Monte Carlo)", {
  expect_identical(cv_smm(.05, df = 36, n_comparisons = 5)$value,
                   cv_smm(.05, df = 36, n_comparisons = 5)$value)
})

test_that("cv_smm() returns the documented data.frame structure", {
  res <- cv_smm(alpha_level = .05, df = 36, n_comparisons = 5)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value", "area_less", "area_greater"))
  expect_equal(res$term, "upper_cv")
  expect_equal(res$area_less + res$area_greater, 1)
})

test_that("cv_smm() critical value increases with more comparisons", {
  v3 <- cv_smm(alpha_level = .05, df = 36, n_comparisons = 3)$value
  v8 <- cv_smm(alpha_level = .05, df = 36, n_comparisons = 8)$value
  expect_gt(v8, v3)
})

test_that("cv_smm() errors on bad inputs", {
  expect_error(cv_smm(),                                            "alpha")
  expect_error(cv_smm(alpha_level = .05),                                 "df")
  expect_error(cv_smm(alpha_level = .05, df = 36),                        "n_comparisons")
  expect_error(cv_smm(alpha_level = 0,   df = 36, n_comparisons = 3),     "alpha")
  expect_error(cv_smm(alpha_level = .05, df = 0,  n_comparisons = 3),     "df")
  expect_error(cv_smm(alpha_level = .05, df = 36, n_comparisons = 0),     "n_comparisons")
  expect_error(cv_smm(alpha_level = .05, df = 36, n_comparisons = 1.5),   "n_comparisons")
})

test_that("cv_smm() works at a large df (the chi-squared bump is not missed)", {
  # Integrating the error variate over (0, Inf) makes an adaptive rule miss the
  # chi-squared mass once df is large, which used to leave the root finder
  # without a bracket and error out. df = 200 is an ordinary sample size.
  for (d in c(200, 500, 1000)) {
    v <- cv_smm(.05, df = d, n_comparisons = 6, verbose = FALSE)$value
    expect_true(is.finite(v))
    # The value must sit between the single-t and the Sidak bound, and approach
    # the df = Inf limit from above as df grows.
    expect_gt(v, qt(.975, d))
    expect_lt(v, qnorm((1 + 0.95^(1/6)) / 2) * 1.05)
  }
  # Large df converges to the known-variance limit.
  expect_equal(cv_smm(.05, df = 1e6, n_comparisons = 6, verbose = FALSE)$value,
               cv_smm(.05, df = Inf, n_comparisons = 6, verbose = FALSE)$value,
               tolerance = 1e-4)
})
