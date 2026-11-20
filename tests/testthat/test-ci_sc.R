# Regression tests for ci_sc(). The ncp branch previously set the
# noncentrality but never the standardized contrast, so the closing
# data.frame() errored with three terms against two values; and a scalar n
# recycled by length(means), which is zero on the contrast-only and ncp-only
# paths, so those calls stopped on the n / c_weights length check.

test_that("ci_sc() works when the noncentrality parameter is supplied directly", {
  c_weights <- c(.5, -.5, -.5, .5)
  res <- ci_sc(ncp = 2.5, s_anova = .8, c_weights = c_weights,
               n_per_group = rep(3, 4), N = 12)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_equal(res$term, c("ci_lower", "std_contrast", "ci_upper"))
  expect_lt(res$value[res$term == "ci_lower"],
            res$value[res$term == "std_contrast"])
  expect_gt(res$value[res$term == "ci_upper"],
            res$value[res$term == "std_contrast"])
  # Every branch defines the noncentrality as the standardized contrast
  # divided by sqrt(sum(c_weights^2 / n)), so the standardized contrast a
  # supplied ncp encodes is ncp * sqrt(sum(c_weights^2 / n)) (hand value:
  # 2.5 * sqrt(1/3)).
  expect_equal(res$value[res$term == "std_contrast"],
               2.5 * sqrt(sum(c_weights^2 / rep(3, 4))),
               tolerance = 1e-12)
})

test_that("ci_sc() means, contrast, and ncp parameterizations of one effect agree", {
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
                     n_per_group = n, N = N)
  via_psi <- ci_sc(contrast = psi_unstd, s_anova = s_anova, c_weights = c_weights,
                   n_per_group = n, N = N)
  via_ncp <- ci_sc(ncp = lambda, s_anova = s_anova, c_weights = c_weights,
                   n_per_group = n, N = N)

  expect_equal(via_psi$value, via_means$value, tolerance = 1e-10)
  expect_equal(via_ncp$value, via_means$value, tolerance = 1e-10)
  expect_equal(via_ncp$value[via_ncp$term == "std_contrast"], 1.25,
               tolerance = 1e-10)
})

test_that("ci_sc() recycles a scalar n on the contrast-only and ncp-only paths", {
  c_weights <- c(.5, -.5, -.5, .5)

  vec_psi <- ci_sc(contrast = 1, s_anova = .8, c_weights = c_weights,
                   n_per_group = rep(3, 4), N = 12)
  scl_psi <- ci_sc(contrast = 1, s_anova = .8, c_weights = c_weights,
                   n_per_group = 3, N = 12)
  expect_equal(scl_psi$value, vec_psi$value, tolerance = 1e-12)

  vec_ncp <- ci_sc(ncp = 2.5, s_anova = .8, c_weights = c_weights,
                   n_per_group = rep(3, 4), N = 12)
  scl_ncp <- ci_sc(ncp = 2.5, s_anova = .8, c_weights = c_weights,
                   n_per_group = 3, N = 12)
  expect_equal(scl_ncp$value, vec_ncp$value, tolerance = 1e-12)
})

test_that("ci_sc() ncp path covers the true standardized contrast at the nominal rate", {
  # Monte Carlo confirmation of the ncp branch. The fast anchor that stays on
  # CRAN is the cross-parameterization equivalence test above (means, contrast,
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
                n_per_group = n, N = N)
    hit[g] <- ci$value[ci$term == "ci_lower"] <= true_std &&
      true_std <= ci$value[ci$term == "ci_upper"]
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
          n_per_group = 3, N = 12),
    "lengths of 'means' and 'c_weights'"
  )
})

test_that("ci_sc() requires a positive s_anova on every path", {
  # Through 1.0.1 s_anova was checked only on the contrast and ncp paths: with
  # means and no s_anova the call failed with "argument is of length zero",
  # and an s_anova of 0 or below was not refused (a negative one reversed
  # the sign of the standardized contrast).
  args <- list(means = c(1, 2, 3), c_weights = c(1, -1, 0),
               n_per_group = c(10, 10, 10))
  expect_error(do.call(ci_sc, args), "standard deviation of the errors")
  for (bad in list(0, -1, NA_real_, c(1, 2))) {
    expect_error(do.call(ci_sc, c(args, list(s_anova = bad))),
                 "single positive, finite number", info = format(bad))
  }
  # None of means, contrast, and ncp failed with "object 'lambda' not found".
  expect_error(ci_sc(s_anova = 1, c_weights = c(1, -1), n_per_group = 10),
               "one of the vector of means")
  ok <- do.call(ci_sc, c(args, list(s_anova = 2)))
  expect_equal(ok$value[ok$term == "std_contrast"], -1 / 2,
               tolerance = 1e-12)
})

test_that("ci_sc() judges contrast weights by value, so integer and double weights agree", {
  # Through 1.0.x the positive weights were compared with 1 by identical(),
  # so c(1L, -1L) and -1:1 were refused although their values pass.
  for (k in list(list(int = c(1L, -1L), dbl = c(1, -1), means = c(5, 3)),
                 list(int = -1:1, dbl = c(-1, 0, 1), means = c(2, 4, 9)))) {
    n <- rep(8, length(k$int))
    via <- function(cw) list(
      ci_sc(means = k$means, s_anova = 2, c_weights = cw, n_per_group = n),
      ci_sc(contrast = sum(k$dbl * k$means), s_anova = 2, c_weights = cw,
            n_per_group = n),
      ci_sc(ncp = 1.7, s_anova = 2, c_weights = cw, n_per_group = n,
            conf_level = NULL, alpha_lower = .01, alpha_upper = .04))
    int <- via(k$int)
    dbl <- via(k$dbl)
    for (i in seq_along(int)) expect_identical(int[[i]]$value, dbl[[i]]$value)
  }
  # Hand value for the two-group means path: the standardized contrast is
  # (5 - 3) / 2 = 1.
  expect_identical(
    ci_sc(means = c(5, 3), s_anova = 2, c_weights = c(1L, -1L),
          n_per_group = 8)$value[2], 1)
})

test_that("ci_sc(), ss_aipe_sc(), and ss_aipe_sc_sensitivity() apply one weights rule", {
  run <- function(f) tryCatch({ f(); "accepted" }, error = function(e)
    paste(class(e)[1], conditionMessage(e)))
  accepted <- list(c(1L, -1L), -1:1, c(.5, .5, -1), c(1/3, 1/3, 1/3, -1),
                   c(.7, .2, .1, -1), c(.1, .2, .7, -1))
  refused <- list(c(1, 1, -2), c(2L, -2L), c(.25, .25, -.5), c(.5, .5, -.5),
                  c(.7, .2, .1 + 1e-6, -1 - 1e-6), c(0, 0, 0), c("1", "-1"),
                  c(1, NA, -1))
  for (cw in c(accepted, refused)) {
    by_ci <- run(function() ci_sc(means = seq_along(cw), s_anova = 1,
                                  c_weights = cw, n_per_group = 10))
    by_plan <- run(function() ss_aipe_sc(true_std_contrast = .5,
                                         c_weights = cw, width = 2))
    by_sens <- run(function() ss_aipe_sc_sensitivity(
      true_std_contrast = .5, specified_n_per_group = 5, c_weights = cw,
      width = 2, G = 2, print_iter = FALSE, seed = 1))
    expect_identical(by_plan, by_ci, info = format(cw))
    expect_identical(by_sens, by_ci, info = format(cw))
    if (by_ci != "accepted")
      expect_match(by_ci, "^dmar_argument_invalid", info = format(cw))
  }
  # The positive weights sum to 1 within 1e-8 in any order: 1.0.x refused
  # c(.7, .2, .1, -1), whose floating-point sum is 1 - 1.1e-16, and
  # accepted the same weights in the order c(.1, .2, .7, -1).
  for (cw in accepted) expect_identical(
    run(function() ci_sc(means = seq_along(cw), s_anova = 1, c_weights = cw,
                         n_per_group = 10)), "accepted", info = format(cw))
  for (cw in refused) expect_false(identical(
    run(function() ci_sc(means = seq_along(cw), s_anova = 1, c_weights = cw,
                         n_per_group = 10)), "accepted"), info = format(cw))
  expect_equal(
    ci_sc(means = c(4, 2, 3, 1), s_anova = 1, c_weights = c(.7, .2, .1, -1),
          n_per_group = 10)$value,
    ci_sc(means = c(3, 2, 4, 1), s_anova = 1, c_weights = c(.1, .2, .7, -1),
          n_per_group = 10)$value, tolerance = 1e-12)
})
