## ss_aipe_pcm() -- sample size for a polynomial change model (linear trend by default).

test_that("ss_aipe_pcm() returns a tidy data.frame with necessary_n_per_group", {
  res <- ss_aipe_pcm(variance_trend = 1, error_variance = 1,
                     duration = 6, frequency = 3, width = 0.20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_per_group" %in% res$term)
  expect_gt(res$value[res$term == "necessary_n_per_group"], 0)
})

test_that("ss_aipe_pcm() requires more N for a tighter width target", {
  wide  <- ss_aipe_pcm(1, 1, duration = 6, frequency = 3, width = 0.30)$value[1]
  tight <- ss_aipe_pcm(1, 1, duration = 6, frequency = 3, width = 0.10)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_pcm() floors an impossibly wide target at the admissible minimum", {
  # A huge target width previously crashed the fixed-point search ("missing
  # value where TRUE/FALSE needed") because the per-group n dropped to 1, making
  # the error df 2n - 2 = 0. The smallest admissible design is n = 2 per group.
  r <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                   duration = 4, frequency = 1, width = 100)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})

test_that("ss_aipe_pcm() validates its boundary inputs", {
  expect_error(ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                           duration = 4, frequency = 1, width = 0.025,
                           conf_level = 1.2), "conf_level")
  expect_error(ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                           duration = 4, frequency = 1, width = 0), "width")
  expect_error(ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                           duration = 4, frequency = 1, width = -1), "width")
  expect_error(ss_aipe_pcm(variance_trend = -1, error_variance = 0.0262,
                           duration = 4, frequency = 1, width = 0.025),
               "variance_trend")
  expect_error(ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                           duration = 4, frequency = 1, width = 0.025,
                           assurance = 1.5), "assurance")
})

test_that("ss_aipe_pcm() requires more N for larger error variance", {
  low_err  <- ss_aipe_pcm(1, 0.5, duration = 6, frequency = 3, width = 0.20)$value[1]
  high_err <- ss_aipe_pcm(1, 2.0, duration = 6, frequency = 3, width = 0.20)$value[1]
  expect_gt(high_err, low_err)
})

test_that("ss_aipe_pcm() reproduces Kelley & Rausch (2011) Table 1 (expected width)", {
  # Tolerance-of-antisocial-thinking example (Kelley & Rausch, 2011, Table 1,
  # p. 400): National Youth Survey data with level-one error variance 0.0262
  # and slope variance 0.003, planning the expected width of a 95% confidence
  # interval on the group-by-time slope difference. T is the number of
  # measurement occasions, M = frequency * duration + 1. The tabled values are
  # per-group sample sizes.
  n <- function(D, f, w) {
    r <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                     duration = D, frequency = f, width = w, conf_level = .95)
    r$value[r$term == "necessary_n_per_group"]
  }
  # T = 3 (D = 2, f = 1)
  expect_equal(n(2, 1, 0.025), 793)
  expect_equal(n(2, 1, 0.050), 200)
  # T = 5 (D = 4, f = 1)
  expect_equal(n(4, 1, 0.025), 278)
  expect_equal(n(4, 1, 0.050), 71)
  # T = 10 (D = 9, f = 1)
  expect_equal(n(9, 1, 0.025), 165)
  expect_equal(n(9, 1, 0.050), 43)
})

test_that("ss_aipe_pcm() reproduces Kelley & Rausch (2011) Table 2 (85% assurance)", {
  # Same design as Table 1, now requiring 85% assurance that the realized 95%
  # confidence interval will be no wider than the target (Kelley & Rausch,
  # 2011, Table 2, p. 401). Per-group sample sizes.
  n <- function(D, f, w) {
    r <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                     duration = D, frequency = f, width = w, conf_level = .95,
                     assurance = .85)
    r$value[r$term == "necessary_n_per_group"]
  }
  # T = 3 (D = 2, f = 1)
  expect_equal(n(2, 1, 0.025), 822)
  expect_equal(n(2, 1, 0.050), 214)
  # T = 5 (D = 4, f = 1)
  expect_equal(n(4, 1, 0.025), 295)
  expect_equal(n(4, 1, 0.050), 79)
  # T = 10 (D = 9, f = 1)
  expect_equal(n(9, 1, 0.025), 178)
  expect_equal(n(9, 1, 0.050), 49)
  # An assurance plan never asks for fewer subjects than the expected-width plan.
  expect_gte(n(4, 1, 0.025),
             ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                         duration = 4, frequency = 1, width = 0.025,
                         conf_level = .95)$value[1])
})

test_that("ss_aipe_pcm() uses the correct cubic K_p constant", {
  # The degree-p change coefficient has sampling variance
  #   Var = error_variance * frequency^(2p) / sum_c2_pm,
  #   sum_c2_pm = K_p * (M+p)! / (M-p-1)!,
  # with K_p = (p!)^2 / [(2p)! (2p+1)!] (Raudenbush & Liu, 2001, p. 392). This
  # test fixes frequency = 1, so frequency^(2p) = 1 and only K_p is under
  # examination. For the cubic, K_3 = (3!)^2 / (6! 7!) = 36 / 3628800 =
  # 1/100800; the buggy literal 1/1000800 (an extra zero in the denominator) is
  # about one tenth of that. Because K_p sits in the denominator of the
  # variance, the typo inflates the cubic variance, and the planned sample
  # size, roughly tenfold. The constant is recovered here directly from an OLS
  # fit, with no reference to K_p: R&L's cubic change coefficient is 3! times
  # the raw t^3 monomial coefficient, so its sampling variance is (3!)^2 times
  # the (X'X)^{-1} cubic diagonal element.
  p <- 3L
  K3 <- factorial(p)^2 / (factorial(2 * p) * factorial(2 * p + 1))
  expect_equal(1 / K3, 100800)
  recovered_sum_c2 <- function(M) {
    t <- 0:(M - 1)
    1 / (factorial(p)^2 * solve(crossprod(outer(t, 0:p, `^`)))[p + 1, p + 1])
  }
  for (M in c(7L, 8L, 10L, 12L)) {
    expect_equal(recovered_sum_c2(M),
                 K3 * factorial(M + p) / factorial(M - p - 1))
  }

  # The planner must reproduce the AIPE sample size implied by that constant.
  hand_n <- function(M, width, conf_level = .95,
                     variance_trend = 0.003, error_variance = 0.0262) {
    vtm <- error_variance / (K3 * factorial(M + p) / factorial(M - p - 1))
    n <- M + 2L
    repeat {
      tc <- stats::qt(1 - (1 - conf_level) / 2, df = 2 * n - 2)
      n_new <- ceiling(8 * (variance_trend + vtm) * tc^2 / width^2)
      if (n_new == n) break
      n <- n_new
    }
    n
  }
  M <- 8L  # duration 7, frequency 1 -> M = f*D + 1 = 8
  expect_equal(
    ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                duration = M - 1L, frequency = 1, width = 0.05,
                trend = "cubic")$value[1],
    hand_n(M, 0.05))
})

test_that("ss_aipe_pcm() accepts variance_true_minus_estimated_trend directly", {
  # For a linear trend with M = 5 (duration 4, frequency 1),
  #   sum_c2_pm = (1/12) * 6! / 3! = 10,  so vtm = error_variance / 10.
  via_error <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                           duration = 4, frequency = 1, width = 0.025)$value[1]
  via_vtm <- ss_aipe_pcm(variance_trend = 0.003,
                         variance_true_minus_estimated_trend = 0.0262 / 10,
                         duration = 4, frequency = 1, width = 0.025)$value[1]
  expect_equal(via_vtm, via_error)
  expect_equal(via_vtm, 278)  # Kelley & Rausch (2011, Table 1, T = 5)
})

test_that("ss_aipe_pcm() validates the error_variance / vtm pair", {
  # Supplying both, consistently, is allowed and matches either alone.
  both_ok <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                         variance_true_minus_estimated_trend = 0.0262 / 10,
                         duration = 4, frequency = 1, width = 0.025)$value[1]
  expect_equal(both_ok, 278)
  # Supplying both, inconsistently, is an error.
  expect_error(
    ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                variance_true_minus_estimated_trend = 0.005,
                duration = 4, frequency = 1, width = 0.025),
    "inconsistent")
  # Supplying neither is an error.
  expect_error(
    ss_aipe_pcm(variance_trend = 0.003, duration = 4, frequency = 1,
                width = 0.025),
    "Specify either")
})

test_that("ss_aipe_pcm() scales the error variance of the slope by frequency^(2p)", {
  # The change coefficient is a per-unit-time rate, so its sampling variance
  # carries a frequency^(2p) factor (Raudenbush & Liu, 2001, p. 392; the same
  # factor ss_power_pcm() puts in V). Every Kelley & Rausch (2011) table uses
  # frequency = 1, where the factor is 1, so the benchmarks are unaffected; the
  # factor only bites for frequency != 1. The anchor here is non-circular: the
  # per-unit-time OLS slope variance for one subject is error_variance / Sxx,
  # with Sxx the corrected sum of squares of the actual time grid.
  hand_n <- function(D, f, width, conf_level = .95,
                     variance_trend = 0.003, error_variance = 0.0262) {
    M <- f * D + 1
    t <- seq(0, D, length.out = M)
    vtm <- error_variance / sum((t - mean(t))^2)  # = error_variance * f^2 / sum_c2_pm
    n <- M + 2
    repeat {
      tc <- stats::qt(1 - (1 - conf_level) / 2, df = 2 * n - 2)
      n_new <- ceiling(8 * (variance_trend + vtm) * tc^2 / width^2)
      if (n_new == n) break
      n <- n_new
    }
    n
  }
  # frequency = 2 (short span) and frequency = 0.5 (long span), both M = 5.
  expect_equal(ss_aipe_pcm(0.003, 0.0262, duration = 2, frequency = 2,
                           width = 0.05)$value[1], hand_n(2, 2, 0.05))
  expect_equal(ss_aipe_pcm(0.003, 0.0262, duration = 8, frequency = 0.5,
                           width = 0.05)$value[1], hand_n(8, 0.5, 0.05))

  # Same M = 5, different frequency: a slope estimated over a longer span needs
  # fewer subjects. Before the frequency fix these were identical.
  n_short <- ss_aipe_pcm(0.003, 0.0262, duration = 2, frequency = 2,
                         width = 0.05)$value[1]
  n_long  <- ss_aipe_pcm(0.003, 0.0262, duration = 8, frequency = 0.5,
                         width = 0.05)$value[1]
  expect_gt(n_short, n_long)

  # Cross-function: ss_power_pcm's reported error variance of the slope (V)
  # equals the variance ss_aipe_pcm now uses, and both equal the direct OLS
  # per-unit-time slope variance.
  pw <- ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
                     frequency = 2, duration = 2, N = 238)
  V_slope <- pw$value[pw$term == "error_var_of_slopes"]
  expect_equal(V_slope, 0.0262 * 2^2 / 10)
  expect_equal(V_slope, 0.0262 / sum((seq(0, 2, length.out = 5) - 1)^2))
})
