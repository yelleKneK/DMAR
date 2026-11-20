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

# On CRAN the file stops here, after the anchor above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

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

test_that("cv_smm() values are unchanged wherever they were right", {
  # Computed by cv_smm() as released in DMAR 1.0.0, before its integral moved
  # to R/smm_internals.R; on a grid of 600 cells every value whose tail area
  # was correct came back bit for bit, and these are compared to 1e-12.
  v <- function(a, d, m) cv_smm(a, df = d, n_comparisons = m, verbose = FALSE)$value
  expect_equal(v(.05, 1.7, 3),    7.74708390064474, tolerance = 1e-12)
  expect_equal(v(.05, 19.625, 3), 2.59833865147556, tolerance = 1e-12)
  expect_equal(v(.05, 36, 6),     2.77485303213572, tolerance = 1e-12)
  expect_equal(v(.05, 1.7, 10),  10.60132304822876, tolerance = 1e-12)
})

test_that("cv_smm() is right at very small df, where DMAR 1.0.0 failed", {
  skip_on_cran()
  # With df below 1 at conventional levels (below about .8 at an alpha_level
  # of .05, .6 at .10), or df of about 2 or less with an alpha_level of .02 or
  # less, the chi squared integral's integrand rises from 0 to 1 in a sliver
  # by the origin that the adaptive rule missed. In DMAR 1.0.0 the calls
  # below at df = 1 and less stopped with "the integral is probably
  # divergent", and the one at df = 2 and alpha_level = .001 returned 32.80,
  # whose tail area is half as large again as alpha_level; the df = 1.25 cell,
  # at the edge of the region, was right then and must stay right. Each value
  # is checked by a trapezoid rule on a fine grid of log S, independent of
  # both integrals in R/smm_internals.R. Below the grid the integrand is 1 to
  # machine precision, so the mass of S there is added in closed form; at
  # df = .3 it is about 2e-4, and at df = .5 and less the check would fail
  # without it.
  tail_brute <- function(q, df, m) {
    s_lo  <- 1e-12
    log_s <- seq(log(s_lo), log(50), length.out = 2e5)
    s <- exp(log_s)
    f <- -expm1(m * log1p(-2 * stats::pnorm(-q * s))) *
      stats::dchisq(df * s^2, df) * 2 * df * s^2
    sum((f[-1] + f[-length(f)]) / 2 * diff(log_s)) +
      stats::pchisq(df * s_lo^2, df)
  }
  cells <- data.frame(alpha = c(.01, .001, .005, .02, .05, .10, .20),
                      df    = c(1,   2,    1.25, 1,   .5,  .5,  .3),
                      m     = c(2,   2,    10,   45,  6,   3,   45))
  for (k in seq_len(nrow(cells))) {
    v <- cv_smm(cells$alpha[k], df = cells$df[k],
                n_comparisons = cells$m[k], verbose = FALSE)$value
    expect_equal(tail_brute(v, cells$df[k], cells$m[k]), cells$alpha[k],
                 tolerance = 1e-6)
  }
  expect_equal(cv_smm(.001, df = 2, n_comparisons = 2, verbose = FALSE)$value,
               40.4342099, tolerance = 1e-8)
})
