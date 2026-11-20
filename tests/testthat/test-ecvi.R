test_that("ecvi() matches the Browne-Cudeck / lavaan point estimate", {
  # (chisq + 2*npar)/n.
  res <- ecvi(chisq = 24.361, df = 8, npar = 13, multiplier = 301)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(res$value[res$term == "ecvi"], (24.361 + 2 * 13) / 301)
  # The interval brackets the point estimate.
  lo <- res$value[res$term == "ci_lower"]
  hi <- res$value[res$term == "ci_upper"]
  est <- res$value[res$term == "ecvi"]
  expect_lte(lo, est)
  expect_gte(hi, est)
  expect_identical(attr(res, "conf_level"), 0.95)
})

test_that("ecvi(fit) agrees with lavaan::fitMeasures", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(
    "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
     verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  res <- ecvi(fit)
  oracle <- unname(lavaan::fitMeasures(fit, "ecvi"))
  expect_equal(res$value[res$term == "ecvi"], oracle, tolerance = 1e-6)
})

# On CRAN the file stops here, after the anchor above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

test_that("ecvi(fit) divides by lavaan's N - G under the Wishart likelihood and GLS", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  model <- "visual =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6
            speed =~ x7 + x8 + x9"
  d <- lavaan::HolzingerSwineford1939
  # lavaan's values: 0.4234070491 (Wishart), 0.7862329 (Wishart, two
  # schools, divisor 299).
  fits <- list(
    wishart = lavaan::cfa(model, data = d, likelihood = "wishart"),
    wishart_groups = lavaan::cfa(model, data = d, likelihood = "wishart",
                                 group = "school"),
    gls = lavaan::cfa(model, data = d, estimator = "GLS"),
    normal = lavaan::cfa(model, data = d))
  for (f in fits) {
    res <- ecvi(f)
    est <- res$value[res$term == "ecvi"]
    expect_equal(est, unname(lavaan::fitMeasures(f, "ecvi")), tolerance = 1e-10)
    expect_lte(res$value[res$term == "ci_lower"], est)
    expect_gte(res$value[res$term == "ci_upper"], est)
  }
  expect_equal(ecvi(fits$wishart)$value[1], 0.4234070491, tolerance = 1e-9)
  expect_equal(ecvi(fits$wishart_groups)$value[1], 0.7862329, tolerance = 1e-6)
  fm <- lavaan::fitMeasures(fits$wishart_groups, c("chisq", "npar"))
  expect_equal(ecvi(fits$wishart_groups)$value[1],
               (fm[["chisq"]] + 2 * fm[["npar"]]) / 299)
  # Under the normal likelihood N ECVI = AIC + 2 l_s.
  fm <- lavaan::fitMeasures(fits$normal, c("aic", "unrestricted.logl"))
  expect_equal(301 * ecvi(fits$normal)$value[1],
               fm[["aic"]] + 2 * fm[["unrestricted.logl"]])
})

test_that("ecvi() validates input", {
  expect_error(ecvi(chisq = 10, df = 0, npar = 5, multiplier = 100), "at least 1")
  expect_error(ecvi(chisq = 10, df = 5, npar = 5, multiplier = -1), "non-negative")
  expect_error(ecvi(chisq = 10, df = 5, npar = 5, multiplier = 100, conf_level = 1),
               "\\(0, 1\\)")
})

test_that("a small statistic gives no stray warning, and its estimate lies below the interval", {
  # T = 5 on 8 df lies below the upper .025 quantile: the lower noncentrality
  # limit is 0 and the lower ECVI limit is (df + 2q)/n, above the untruncated
  # estimate. ci_nc_chisq()'s warning about a 'prob_greater' row ecvi() does
  # not return is not passed on.
  expect_no_warning(x <- ecvi(chisq = 5, df = 8, npar = 13, multiplier = 301))
  expect_equal(x$value[1:2], c((5 + 26) / 301, (8 + 26) / 301))
  expect_lt(x$value[1], x$value[2])
  expect_gt(x$value[3], x$value[2])
  # T = 1.5 lies below the lower .025 quantile (2.18): no noncentrality is
  # consistent with it, and the upper limit is NA with ecvi()'s own warning.
  expect_warning(y <- ecvi(chisq = 1.5, df = 8, npar = 13, multiplier = 301),
                 "upper limit of the ECVI is NA")
  expect_equal(y$value[2], (8 + 26) / 301)
  expect_true(is.na(y$value[3]))
})

test_that("ecvi() refuses a zero divisor and an inadmissible solution", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(ecvi(chisq = 10, df = 5, npar = 5, multiplier = 0), "'multiplier' must be positive")
  skip_if_not_installed("lavaan")
  # A fit with a negative variance, which sem_fit_measures() refuses; ecvi()
  # had returned a value for such a fit (1.116 [0.799, 1.598] for a 60 case
  # two-factor fit of the Holzinger and Swineford tests with a latent
  # variance of -72.9).
  P <- diag(5)
  P[1, 2] <- P[2, 1] <- P[1, 3] <- P[3, 1] <- 0.8
  P[2, 3] <- P[3, 2] <- 0.3
  P[4, 5] <- P[5, 4] <- 0.5
  dimnames(P) <- list(paste0("y", 1:5), paste0("y", 1:5))
  heywood <- suppressWarnings(lavaan::cfa("f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5",
                                          sample.cov = P, sample.nobs = 300))
  # (lavaan's own check warns as it finds the negative variance.)
  expect_false(suppressWarnings(lavaan::lavInspect(heywood, "post.check")))
  expect_error(suppressWarnings(ecvi(heywood)), "inadmissible")
  expect_error(ecvi(stats::lm(t2_cubes ~ t1_visual_perception, data = holzinger_swineford)),
               "lavaan model")
})

test_that("a statistic in the millions gets ECVI limits around the estimate", {
  skip_on_cran()
  # The limits once collapsed to 29.83, below the estimate 30.004.
  x <- ecvi(chisq = 3e6, df = 3000, npar = 200, multiplier = 1e5, conf_level = 0.90)
  expect_lt(x$value[2], x$value[1])
  expect_gt(x$value[3], x$value[1])
})
