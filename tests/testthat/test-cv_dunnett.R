test_that("cv_dunnett() returns a tidy data.frame and the value lies between the t and Tukey crit", {
  res <- cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value", "area_less", "area_greater"))
  # Dunnett crit should be larger than unadjusted two-sided t...
  expect_gt(res$value, qt(.975, df = 36))
  # ...but smaller than Tukey HSD with 4 groups (since Dunnett tests fewer pairs).
  tukey_val <- qtukey(.95, nmeans = 4, df = 36) / sqrt(2)
  expect_lt(res$value, tukey_val)
})

test_that("cv_dunnett() with one comparison reduces exactly to the t critical value", {
  # A single treatment-vs-control comparison is an ordinary t test.
  expect_equal(cv_dunnett(.05, df = 36, n_comparisons = 1)$value,
               qt(.975, df = 36), tolerance = 1e-8)                  # two-sided
  expect_equal(cv_dunnett(.05, df = 36, n_comparisons = 1,
                          alternative = "greater")$value,
               qt(.95, df = 36), tolerance = 1e-8)                   # one-sided
})

test_that("cv_dunnett() at df = Inf reduces to the normal limit for one comparison", {
  expect_equal(cv_dunnett(.05, df = Inf, n_comparisons = 1)$value,
               qnorm(.975), tolerance = 1e-8)
})

test_that("cv_dunnett() two-sided value matches Maxwell-Delaney-Kelley Appendix A.6", {
  # df = 60, m = 3 (a = 4 groups), alpha = .05 -> 2.41 in Table A.6.
  res <- cv_dunnett(alpha_level = .05, df = 60, n_comparisons = 3)
  expect_equal(round(res$value, 2), 2.41)
})

test_that("cv_dunnett() one-sided greater matches Appendix A.7", {
  # df = 60, m = 3, alpha = .05 one-sided -> 2.10 in Table A.7.
  res <- cv_dunnett(alpha_level = .05, df = 60, n_comparisons = 3,
                    alternative = "greater")
  expect_equal(round(res$value, 2), 2.10)
})

test_that("cv_dunnett() is deterministic (no Monte Carlo)", {
  expect_identical(cv_dunnett(.05, df = 36, n_comparisons = 4)$value,
                   cv_dunnett(.05, df = 36, n_comparisons = 4)$value)
})

test_that("cv_dunnett() one-sided less is the negative of greater", {
  g <- cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3, alternative = "greater")
  l <- cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3, alternative = "less")
  expect_equal(l$value, -g$value)
  expect_equal(l$term, "lower_cv")
  expect_equal(g$term, "upper_cv")
})

test_that("cv_dunnett() critical value increases with more comparisons", {
  v2 <- cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 2)$value
  v6 <- cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 6)$value
  expect_gt(v6, v2)
})

test_that("cv_dunnett() errors on bad inputs", {
  expect_error(cv_dunnett(),                                                       "alpha")
  expect_error(cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 3,
                          alternative = "weird"),                                  "alternative")
  expect_error(cv_dunnett(alpha_level = 0,   df = 36, n_comparisons = 3),                "alpha")
  expect_error(cv_dunnett(alpha_level = .05, df = 0,  n_comparisons = 3),                "df")
  expect_error(cv_dunnett(alpha_level = .05, df = 36, n_comparisons = 0),                "n_comparisons")
})

test_that("cv_dunnett() works at a large df (the chi-squared bump is not missed)", {
  skip_on_cran()  # a critical value at four large df; the A.6 and A.7 anchors run on CRAN
  # See the companion test in test-cv_smm.R: over (0, Inf) the error-variate
  # integral silently returned zero once df was large, and the root finder had
  # no bracket. df = 200 is an ordinary sample size.
  for (d in c(200, 500, 1000)) {
    v <- cv_dunnett(.05, df = d, n_comparisons = 3, verbose = FALSE)$value
    expect_true(is.finite(v))
    expect_gt(v, qt(.975, d))
  }
  # Large df converges to the known-variance limit.
  expect_equal(cv_dunnett(.05, df = 1e6, n_comparisons = 3, verbose = FALSE)$value,
               cv_dunnett(.05, df = Inf, n_comparisons = 3, verbose = FALSE)$value,
               tolerance = 1e-4)
})
