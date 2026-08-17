# Regression tests for ci_sc(). The ncp branch previously set the
# noncentrality but never the standardized contrast, so the closing
# data.frame() errored with three terms against two values; and a scalar n
# recycled by length(means), which is zero on the psi-only and ncp-only
# paths, so those calls stopped on the n / c_weights length check.

test_that("ci_sc() works when the noncentrality parameter is supplied directly", {
  c_weights <- c(.5, -.5, -.5, .5)
  res <- ci_sc(ncp = 2.5, s_anova = .8, c_weights = c_weights,
               n = rep(3, 4), N = 12)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_equal(res$term, c("lower_limit", "std_contrast", "upper_limit"))
  expect_lt(res$value[res$term == "lower_limit"],
            res$value[res$term == "std_contrast"])
  expect_gt(res$value[res$term == "upper_limit"],
            res$value[res$term == "std_contrast"])
  # Every branch defines the noncentrality as the standardized contrast
  # divided by sqrt(sum(c_weights^2 / n)), so the standardized contrast a
  # supplied ncp encodes is ncp * sqrt(sum(c_weights^2 / n)) (hand value:
  # 2.5 * sqrt(1/3)).
  expect_equal(res$value[res$term == "std_contrast"],
               2.5 * sqrt(sum(c_weights^2 / rep(3, 4))),
               tolerance = 1e-12)
})

test_that("ci_sc() means, psi, and ncp parameterizations of one effect agree", {
  means <- c(2, 4, 9, 13)
  s_anova <- .8
  c_weights <- c(.5, -.5, -.5, .5)
  n <- rep(3, 4)
  N <- 12

  # Hand computation of the noncentrality the means imply: the
  # unstandardized contrast is sum(c_weights * means) = 1, the standardized
  # contrast is 1 / .8 = 1.25, and lambda = 1.25 / sqrt(1/3).
  psi_unstd <- sum(c_weights * means)
  lambda <- (psi_unstd / s_anova) / sqrt(sum(c_weights^2 / n))

  via_means <- ci_sc(means = means, s_anova = s_anova, c_weights = c_weights,
                     n = n, N = N)
  via_psi <- ci_sc(psi = psi_unstd, s_anova = s_anova, c_weights = c_weights,
                   n = n, N = N)
  via_ncp <- ci_sc(ncp = lambda, s_anova = s_anova, c_weights = c_weights,
                   n = n, N = N)

  expect_equal(via_psi$value, via_means$value, tolerance = 1e-10)
  expect_equal(via_ncp$value, via_means$value, tolerance = 1e-10)
  expect_equal(via_ncp$value[via_ncp$term == "std_contrast"], 1.25,
               tolerance = 1e-10)
})

test_that("ci_sc() recycles a scalar n on the psi-only and ncp-only paths", {
  c_weights <- c(.5, -.5, -.5, .5)

  vec_psi <- ci_sc(psi = 1, s_anova = .8, c_weights = c_weights,
                   n = rep(3, 4), N = 12)
  scl_psi <- ci_sc(psi = 1, s_anova = .8, c_weights = c_weights,
                   n = 3, N = 12)
  expect_equal(scl_psi$value, vec_psi$value, tolerance = 1e-12)

  vec_ncp <- ci_sc(ncp = 2.5, s_anova = .8, c_weights = c_weights,
                   n = rep(3, 4), N = 12)
  scl_ncp <- ci_sc(ncp = 2.5, s_anova = .8, c_weights = c_weights,
                   n = 3, N = 12)
  expect_equal(scl_ncp$value, vec_ncp$value, tolerance = 1e-12)
})

test_that("ci_sc() ncp path covers the true standardized contrast at the nominal rate", {
  # Monte Carlo confirmation of the ncp branch. The fast anchor that stays on
  # CRAN is the cross-parameterization equivalence test above (means, psi,
  # and ncp inputs of one effect return identical tables).
  skip_on_cran()
  set.seed(113)
  c_weights <- c(.5, -.5, -.5, .5)
  n <- rep(10, 4)
  N <- 40
  mu <- c(2, 4, 9, 13)
  sigma <- 4
  true_std <- sum(c_weights * mu) / sigma
  part_of_se <- sqrt(sum(c_weights^2 / n))
  G <- 2000
  hit <- logical(G)
  for (g in seq_len(G)) {
    y <- lapply(seq_along(n), function(j) rnorm(n[j], mu[j], sigma))
    m <- vapply(y, mean, numeric(1))
    s2_pooled <- sum(vapply(y, function(v) sum((v - mean(v))^2), numeric(1))) / (N - 4)
    t_obs <- sum(c_weights * m) / (sqrt(s2_pooled) * part_of_se)
    ci <- ci_sc(ncp = t_obs, s_anova = sqrt(s2_pooled), c_weights = c_weights,
                n = n, N = N)
    hit[g] <- ci$value[ci$term == "lower_limit"] <= true_std &&
      true_std <= ci$value[ci$term == "upper_limit"]
  }
  # Under this seed the realized coverage is 0.951; the band is about three
  # Monte Carlo standard errors around the nominal .95.
  expect_gt(mean(hit), 0.935)
  expect_lt(mean(hit), 0.965)
})

test_that("ci_sc() rejects a means vector whose length differs from c_weights", {
  # With scalar n now recycled by length(c_weights), a short means vector
  # would otherwise recycle silently inside sum(c_weights * means).
  expect_error(
    ci_sc(means = c(2, 4, 9), s_anova = .8, c_weights = c(.5, .5, -.5, -.5),
          n = 3, N = 12),
    "lengths of 'means' and 'c_weights'"
  )
})
