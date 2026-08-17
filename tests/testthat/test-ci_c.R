# Regression test for ci_c(): the df_error argument is documented for
# multi-factor designs but previously left df_2 undefined when supplied,
# so the qt() calls errored. Supplying df_error must now work and must
# change the interval relative to the default one-way df.

test_that("ci_c() honors a user-supplied df_error", {
  means <- c(10, 12, 14)
  args  <- list(means = means, s_anova = 4, c_weights = c(-1, 0, 1),
                n = c(20, 20, 20))
  expect_error(do.call(ci_c, c(args, list(df_error = 50))), NA)

  default_ci <- do.call(ci_c, args)
  custom_ci  <- do.call(ci_c, c(args, list(df_error = 50)))
  # Default one-way df is N - groups = 60 - 3 = 57; df_error = 50 is smaller,
  # so the interval is wider.
  w_default <- diff(range(default_ci$value[default_ci$term %in%
                                           c("lower_limit", "upper_limit")]))
  w_custom  <- diff(range(custom_ci$value[custom_ci$term %in%
                                          c("lower_limit", "upper_limit")]))
  expect_gt(w_custom, w_default)
})

test_that("ci_c() works on the psi-only path with a scalar n", {
  # Regression: a scalar n was recycled by length(means), which is zero when
  # only psi is supplied, so the documented psi-only call stopped on the
  # n / c_weights length check.
  c_weights <- c(.5, -.5, -.5, .5)
  res <- ci_c(psi = 3, s_anova = .8, c_weights = c_weights, n = 3, N = 12)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, c("lower_limit", "contrast", "upper_limit"))

  # Hand anchor: psi +/- qt(.975, df = N - 4) * s_anova * sqrt(sum(c^2 / n)),
  # with df = 12 - 4 = 8.
  se <- .8 * sqrt(sum(c_weights^2 / rep(3, 4)))
  expect_equal(res$value[res$term == "contrast"], 3, tolerance = 1e-12)
  expect_equal(res$value[res$term == "lower_limit"], 3 - qt(.975, 8) * se,
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "upper_limit"], 3 + qt(.975, 8) * se,
               tolerance = 1e-10)

  # A scalar n agrees with the equivalent explicit vector n.
  vec <- ci_c(psi = 3, s_anova = .8, c_weights = c_weights,
              n = rep(3, 4), N = 12)
  expect_equal(res$value, vec$value, tolerance = 1e-12)
})

test_that("ci_c() rejects a means vector whose length differs from c_weights", {
  # With scalar n now recycled by length(c_weights), a short means vector
  # would otherwise recycle silently inside sum(c_weights * means).
  expect_error(
    ci_c(means = c(10, 12, 14), s_anova = 4, c_weights = c(.5, .5, -.5, -.5),
         n = 20),
    "lengths of 'means' and 'c_weights'"
  )
})

test_that("ci_c() reproduces the hand t-based interval psi +/- qt * SE", {
  means <- c(10, 12, 14); s_anova <- 4
  c_weights <- c(-1, 0, 1); n <- c(20, 20, 20)

  psi     <- sum(c_weights * means)                 # 4
  N       <- sum(n)                                 # 60
  df_e    <- N - length(c_weights)                  # 57
  se_part <- sqrt(sum(c_weights^2 / n))
  cv      <- qt(0.975, df = df_e)
  lo_hand <- psi - cv * se_part * s_anova
  hi_hand <- psi + cv * se_part * s_anova

  res <- ci_c(means = means, s_anova = s_anova,
              c_weights = c_weights, n = n)
  expect_equal(res$value[res$term == "contrast"], psi, tolerance = 1e-12)
  expect_equal(res$value[res$term == "lower_limit"], lo_hand,
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "upper_limit"], hi_hand,
               tolerance = 1e-10)
})
