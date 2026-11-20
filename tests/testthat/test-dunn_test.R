bp <- data.frame(
  y = c(84, 95, 93, 104, 99, 106,
        81, 84, 92, 101, 80, 108,
        98, 95, 86, 87, 94, 101,
        91, 78, 85, 80, 81, 76),
  g = factor(rep(1:4, each = 6))
)

test_that("dunn_test() reproduces a hand computation of Dunn (1964)", {
  # These data contain ties (84, 95, 101, 80, and 81 each occur twice), so the
  # tie correction applies and the variance term is below the untied N(N+1)/12.
  out <- dunn_test(bp$y, bp$g, method = "none")
  r  <- rank(bp$y)
  mr <- tapply(r, bp$g, mean)
  N  <- 24
  tc <- as.numeric(table(bp$y))
  s2 <- N * (N + 1) / 12 - sum(tc^3 - tc) / (12 * (N - 1))
  se_expected <- sqrt(s2 * (1 / 6 + 1 / 6))
  z_expected  <- (mr[[4]] - mr[[1]]) / se_expected
  row <- out[out$contrast == "4 - 1", ]
  expect_equal(row$se, se_expected, tolerance = 1e-10)
  expect_equal(row$z_statistic, z_expected, tolerance = 1e-10)
  expect_equal(row$p_value, 2 * pnorm(-abs(z_expected)), tolerance = 1e-12)
})

test_that("dunn_test() with two groups matches the Wilcoxon normal approximation", {
  # With a = 2 and no ties, Dunn's z is the standardized rank-sum statistic,
  # so it agrees with the normal approximation to the Wilcoxon rank-sum test
  # (without the continuity correction wilcox.test applies by default).
  d  <- droplevels(subset(bp, g %in% c("1", "4")))
  out <- dunn_test(d$y, d$g, method = "none")
  r   <- rank(d$y)
  W   <- sum(r[d$g == "1"])
  n1  <- 6; n2 <- 6; N <- 12
  z_w <- (W - n1 * (N + 1) / 2) / sqrt(n1 * n2 * (N + 1) / 12)
  expect_equal(abs(out$z_statistic), abs(z_w), tolerance = 1e-10)
})

test_that("dunn_test() applies the tie correction", {
  # Duplicate values force ties; the tie term must shrink the variance, so the
  # standard error is smaller than the untied N(N+1)/12 form.
  tied <- data.frame(y = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3),
                     g = factor(rep(1:2, each = 6)))
  out  <- dunn_test(tied$y, tied$g, method = "none")
  N <- 12
  se_untied <- sqrt((N * (N + 1) / 12) * (1 / 6 + 1 / 6))
  expect_lt(out$se, se_untied)
})

test_that("dunn_test() returns the documented structure", {
  out <- dunn_test(bp$y, bp$g)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("contrast", "mean_rank_difference", "se",
                      "z_statistic", "p_value", "p_adjusted"))
  expect_equal(nrow(out), 6L)
  expect_equal(attr(out, "adjustment_method"), "holm")
})

test_that("dunn_test() adjusts p-values as p.adjust does", {
  raw <- dunn_test(bp$y, bp$g, method = "none")
  adj <- dunn_test(bp$y, bp$g, method = "holm")
  expect_equal(adj$p_adjusted, p.adjust(raw$p_value, "holm"))
  expect_true(all(adj$p_adjusted >= raw$p_value))
  # 'none' leaves the values alone.
  expect_equal(raw$p_adjusted, raw$p_value)
})

test_that("dunn_test() accepts a fitted model and matches the vectors", {
  a <- dunn_test(bp$y, bp$g)
  b <- dunn_test(stats::aov(y ~ g, data = bp))
  expect_equal(a$z_statistic, b$z_statistic)
})

test_that("dunn_test() rejects bad input", {
  expect_error(dunn_test(bp$y, bp$g, method = "bogus"), "p.adjust.methods")
  expect_error(dunn_test(bp$y, bp$g[-1]), "same length")
  expect_error(dunn_test(bp$y), "numeric vector")
})

test_that("dunn_test() errors on all-tied data instead of returning NaN", {
  # Every observation tied at one value makes the tie correction cancel the
  # untied variance exactly, so sigma2 = 0 and every z and p-value would be an
  # undefined 0/0. Refuse the degenerate input with an informative error.
  all_tied <- data.frame(y = rep(5, 12), g = factor(rep(1:3, each = 4)))
  expect_error(dunn_test(all_tied$y, all_tied$g), "tie-corrected|tied")
  # A single distinct value across the whole sample is the trigger, even with
  # only two groups.
  expect_error(dunn_test(rep(0, 8), factor(rep(1:2, each = 4))), "tied")
})
