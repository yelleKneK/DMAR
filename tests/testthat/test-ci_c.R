# Regression test for ci_c(): the df_error argument is documented for
# multi-factor designs but previously left df_2 undefined when supplied,
# so the qt() calls errored. Supplying df_error must now work and must
# change the interval relative to the default one-way df.

test_that("ci_c() honors a user-supplied df_error", {
  means <- c(10, 12, 14)
  args  <- list(means = means, s_anova = 4, c_weights = c(-1, 0, 1),
                n_per_group = c(20, 20, 20))
  expect_error(do.call(ci_c, c(args, list(df_error = 50))), NA)

  default_ci <- do.call(ci_c, args)
  custom_ci  <- do.call(ci_c, c(args, list(df_error = 50)))
  # Default one-way df is N - groups = 60 - 3 = 57; df_error = 50 is smaller,
  # so the interval is wider.
  w_default <- diff(range(default_ci$value[default_ci$term %in%
                                           c("ci_lower", "ci_upper")]))
  w_custom  <- diff(range(custom_ci$value[custom_ci$term %in%
                                          c("ci_lower", "ci_upper")]))
  expect_gt(w_custom, w_default)
})

test_that("ci_c() works on the contrast-only path with a scalar n", {
  # Regression: a scalar n was recycled by length(means), which is zero when
  # only contrast is supplied, so the documented contrast-only call stopped on the
  # n / c_weights length check.
  c_weights <- c(.5, -.5, -.5, .5)
  res <- ci_c(contrast = 3, s_anova = .8, c_weights = c_weights, n_per_group = 3, N = 12)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, c("ci_lower", "contrast", "ci_upper"))

  # Hand anchor: psi +/- qt(.975, df = N - 4) * s_anova * sqrt(sum(c^2 / n)),
  # with df = 12 - 4 = 8.
  se <- .8 * sqrt(sum(c_weights^2 / rep(3, 4)))
  expect_equal(res$value[res$term == "contrast"], 3, tolerance = 1e-12)
  expect_equal(res$value[res$term == "ci_lower"], 3 - qt(.975, 8) * se,
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "ci_upper"], 3 + qt(.975, 8) * se,
               tolerance = 1e-10)

  # A scalar n agrees with the equivalent explicit vector n.
  vec <- ci_c(contrast = 3, s_anova = .8, c_weights = c_weights,
              n_per_group = rep(3, 4), N = 12)
  expect_equal(res$value, vec$value, tolerance = 1e-12)
})

test_that("ci_c() rejects a means vector whose length differs from c_weights", {
  # With scalar n now recycled by length(c_weights), a short means vector
  # would otherwise recycle silently inside sum(c_weights * means).
  expect_error(
    ci_c(means = c(10, 12, 14), s_anova = 4, c_weights = c(.5, .5, -.5, -.5),
         n_per_group = 20),
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
              c_weights = c_weights, n_per_group = n)
  expect_equal(res$value[res$term == "contrast"], psi, tolerance = 1e-12)
  expect_equal(res$value[res$term == "ci_lower"], lo_hand,
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "ci_upper"], hi_hand,
               tolerance = 1e-10)
})

test_that("ci_c() requires a positive s_anova on the means path as on the contrast path", {
  # Through 1.0.1 s_anova was checked only when the contrast was given: with means
  # and no s_anova each limit was psi plus a product with NULL, which
  # data.frame() recycled, so the table gave a zero-width interval at the
  # estimate (-1, -1, -1 here) labeled "Confidence level: 95%" and a
  # sentence "contrast = -1.00, 95% CI [-1.00, -1.00]".
  args <- list(means = c(1, 2, 3), c_weights = c(1, -1, 0), n_per_group = c(10, 10, 10))
  expect_error(do.call(ci_c, args), "standard deviation of the errors")
  expect_error(do.call(ci_c, c(args, list(alpha_lower = .01,
                                          alpha_upper = .04))),
               "standard deviation of the errors")
  expect_error(ci_c(contrast = 1, c_weights = c(1, -1), n_per_group = 10),
               "standard deviation of the errors")
  # An s_anova of 0 gave the same zero-width interval, and a negative one
  # swapped the limits; neither is a standard deviation.
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(do.call(ci_c, c(args, list(s_anova = bad))),
                 "single positive, finite number", info = format(bad))
  }
  # Neither means nor contrast failed with "arguments imply differing number of
  # rows: 3, 0".
  expect_error(ci_c(s_anova = 1, c_weights = c(1, -1), n_per_group = 10),
               "either the vector of means")
  expect_error(ci_c(means = c(1, 2), s_anova = 1, n_per_group = 10),
               "contrast weights")
  # A valid s_anova still gives the interval it did.
  ok <- do.call(ci_c, c(args, list(s_anova = 2)))
  se <- 2 * sqrt(2 / 10)
  expect_equal(ok$value, -1 + c(-1, 0, 1) * stats::qt(.975, 27) * se,
               tolerance = 1e-12)
})

test_that("ci_c() refuses one tail area given without the other", {
  # DMAR 1.0.0 stopped inside data.frame() ("arguments imply differing
  # number of rows: 3, 2") or qt() ("Non-numeric argument to mathematical
  # function").
  for (tails in list(list(alpha_upper = .025), list(alpha_lower = .025),
                     list(alpha_upper = 1e-20))) {
    e <- expect_error(do.call(ci_c, c(list(means = c(1, 0), s_anova = 1,
                                          c_weights = c(1, -1),
                                          n_per_group = c(20, 20),
                                          conf_level = NULL), tails)),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, c("alpha_lower", "alpha_upper"))
    expect_match(conditionMessage(e), paste0("^Only '", names(tails), "' was given"))
  }
})
