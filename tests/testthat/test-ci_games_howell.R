# The blood pressure data of Maxwell, Delaney, and Kelley (2027, Chapter 5,
# Table 5.2): 24 mild hypertensives randomly assigned to four treatments. The
# chapter uses these data to illustrate Tukey's HSD and its modifications for
# unequal variances, which is what ci_games_howell() computes.
bp <- data.frame(
  y = c(84, 95, 93, 104, 99, 106,     # 1 drug therapy
        81, 84, 92, 101, 80, 108,     # 2 biofeedback
        98, 95, 86, 87, 94, 101,      # 3 diet
        91, 78, 85, 80, 81, 76),      # 4 combination
  g = factor(rep(1:4, each = 6))
)

test_that("the blood pressure group summaries match MDK (2027) Table 5.2", {
  expect_equal(round(as.numeric(tapply(bp$y, bp$g, mean)), 3),
               c(96.833, 91.000, 93.500, 81.833))
  expect_equal(round(as.numeric(tapply(bp$y, bp$g, var)), 3),
               c(64.567, 132.000, 35.500, 29.367))
})

test_that("ci_games_howell() uses the Welch-Satterthwaite df of MDK Eq. 5.14", {
  out <- ci_games_howell(bp$y, bp$g)
  # Hand-compute the df for the 2 - 1 contrast from Equation 5.14.
  v1 <- 64.567 / 6; v2 <- 132.000 / 6
  df_expected <- (v1 + v2)^2 / (v1^2 / 5 + v2^2 / 5)
  expect_equal(out$df[out$contrast == "2 - 1"], df_expected, tolerance = 1e-3)
})

test_that("ci_games_howell() with two groups is exactly Welch's t test", {
  # q_{alpha; 2, df} = sqrt(2) * t_{1-alpha/2, df}, so the Games-Howell
  # interval and p-value coincide with t.test(var.equal = FALSE).
  d  <- droplevels(subset(bp, g %in% c("1", "4")))
  gh <- ci_games_howell(d$y, d$g)
  tt <- stats::t.test(y ~ g, data = d, var.equal = FALSE)
  # t.test reports mean_1 - mean_4; ci_games_howell reports 4 - 1.
  expect_equal(gh$mean_difference, -as.numeric(diff(rev(tt$estimate))),
               tolerance = 1e-8)
  expect_equal(sort(c(gh$lower_limit, gh$upper_limit)),
               sort(as.numeric(-tt$conf.int)), tolerance = 1e-6)
  # The identity holds to about 1e-10 in absolute terms; the tolerance here is
  # relative, and what limits the agreement is stats::ptukey()'s own precision,
  # not the procedure.
  expect_equal(gh$p_adjusted, tt$p.value, tolerance = 1e-6)
  expect_equal(gh$df, as.numeric(tt$parameter), tolerance = 1e-8)
})

test_that("ci_games_howell() returns the documented structure", {
  out <- ci_games_howell(bp$y, bp$g)
  expect_s3_class(out, "dmar_post_hoc_ci")
  expect_named(out, c("contrast", "mean_difference", "se", "df",
                      "q_statistic", "lower_limit", "upper_limit",
                      "p_adjusted"))
  expect_equal(nrow(out), 6L)                    # a(a-1)/2 with a = 4
  expect_true(all(out$lower_limit < out$upper_limit))
  expect_equal(attr(out, "conf_level"), 0.95)
})

test_that("ci_games_howell() accepts a fitted model and matches the vectors", {
  a <- ci_games_howell(bp$y, bp$g)
  b <- ci_games_howell(stats::aov(y ~ g, data = bp))
  expect_equal(a$mean_difference, b$mean_difference)
  expect_equal(a$p_adjusted, b$p_adjusted)
})

test_that("ci_games_howell() differs from Tukey-Kramer when variances differ", {
  # Group 2's variance is about 4.5 times group 4's, so the pooled error term
  # is not interchangeable with the separate ones.
  gh <- ci_games_howell(bp$y, bp$g)
  tk <- ci_tukey_kramer(bp$y, bp$g)
  expect_equal(gh$mean_difference, tk$mean_difference)   # same point estimates
  expect_false(isTRUE(all.equal(gh$p_adjusted, tk$p_adjusted)))
  # Per-pair Welch df are smaller than the pooled 20 error df.
  expect_true(all(gh$df < 20))
})

test_that("ci_games_howell() tidies and glances through the shared methods", {
  out <- ci_games_howell(bp$y, bp$g)
  td  <- generics::tidy(out)
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper",
                     "p_adjusted", "conf_level"))
  expect_equal(nrow(td), 6L)
  expect_equal(generics::glance(out)$n_contrasts, 6L)
})

test_that("ci_games_howell() rejects bad input", {
  expect_error(ci_games_howell(bp$y, bp$g, conf_level = 0), "conf_level")
  expect_error(ci_games_howell(bp$y, bp$g[-1]), "same length")
  expect_error(ci_games_howell(bp$y), "numeric vector")
  expect_error(ci_games_howell(c(1, 2), factor(c("a", "b"))), "two observations")
})

test_that("ci_games_howell() errors on a pair of constant groups, naming it", {
  # Both groups essentially constant makes the pairwise sampling variance zero,
  # so the standard error is zero and the Welch df an undefined 0/0. Base R's
  # t.test() errors "data are essentially constant" on the same inputs; match
  # that spirit and name the offending contrast.
  expect_error(t.test(c(1, 1, 1), c(2, 2, 2)), "constant")   # the oracle
  y <- c(1, 1, 1, 2, 2, 2)
  g <- factor(rep(c("a", "b"), each = 3))
  expect_error(ci_games_howell(y, g), "constant")
  expect_error(ci_games_howell(y, g), "b - a")
})

test_that("ci_games_howell() still works when only one group is constant", {
  # Welch's t handles a single constant group, so this pair is not degenerate.
  y <- c(1, 1, 1, 2, 3, 4)
  g <- factor(rep(c("a", "b"), each = 3))
  out <- ci_games_howell(y, g)
  expect_equal(out$se, sqrt((0 / 3 + var(c(2, 3, 4)) / 3) / 2), tolerance = 1e-10)
  expect_false(is.nan(out$df))
})
