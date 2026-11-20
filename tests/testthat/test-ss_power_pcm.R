test_that("ss_power_pcm() reproduces the Raudenbush & Liu (2001) NYS benchmark", {
  # National Youth Survey example (Raudenbush & Liu, 2001, p. 393, Table 1):
  # frequency = 1, duration = 4 (M = 5), level-1 variance 0.0262, slope
  # variance tau = 0.003, standardized slope difference -0.40, n = 238. The
  # paper's hand computation gives sampling variance V = 0.00262, slope
  # reliability 0.53, noncentrality 5.046, and power 0.61. DMAR matches
  # MBESS::ss.power.pcm to six decimals (Actual.Power = 0.612236).
  res <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003, true_var_error = .0262,
                      frequency = 1, duration = 4, specified_n_per_group = 119)
  v <- function(t) res$value[res$term == t]
  expect_equal(round(v("true_var_trend_error"), 5), 0.00262)
  # The oracle by hand: Kelley and Rausch (2011, Equation 19) give the
  # sampling-error variance of one individual's estimated linear change
  # coefficient as sigma^2_e / sum(c^2), with the linear weights -2:2 at
  # M = 5 (sum of squares 10).
  expect_equal(v("true_var_trend_error"), .0262 / sum((-2:2)^2), tolerance = 1e-15)
  expect_equal(round(v("reliability"), 2), 0.53)
  expect_equal(round(v("actual_power"), 2), 0.61)
  expect_equal(v("total_N"), 238)
  expect_equal(v("n_occasions"), 5)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_power_pcm() returns a tidy data.frame with the documented columns", {
  res <- ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_per_group", "necessary_n_per_group", "total_N", "actual_power",
                    "n_occasions") %in% res$term))
})

test_that("ss_power_pcm() actual_power at returned n meets desired_power", {
  res <- ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  expect_gte(res$value[res$term == "actual_power"], 0.80 - 1e-6)
})

test_that("ss_power_pcm() total_N = ss_c + ss_t", {
  res <- ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  ss_c <- res$value[res$term == "necessary_n_per_group"]
  ss_t <- res$value[res$term == "necessary_n_per_group"]
  expect_equal(res$value[res$term == "total_N"], ss_c + ss_t)
})

test_that("ss_power_pcm() n_occasions = frequency * duration + 1", {
  res <- ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  expect_equal(res$value[res$term == "n_occasions"], 5 * 4 + 1)
})

test_that("ss_power_pcm() necessary sample size grows as beta shrinks", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  big_eff   <- ss_power_pcm(true_std_reg_coef = 0.4, true_var_trend = 0.05, true_var_error = 1,
                            frequency = 5, duration = 4,
                            desired_power = 0.80)$value[2]
  small_eff <- ss_power_pcm(true_std_reg_coef = 0.1, true_var_trend = 0.05, true_var_error = 1,
                            frequency = 5, duration = 4,
                            desired_power = 0.80)$value[2]
  expect_gt(small_eff, big_eff)
})

test_that("ss_power_pcm() necessary sample size grows as duration shrinks (less follow-up)", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  long  <- ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                        frequency = 5, duration = 4,
                        desired_power = 0.80)$value[2]
  short <- ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                        frequency = 5, duration = 1,
                        desired_power = 0.80)$value[2]
  expect_gt(short, long)
})

test_that("ss_power_pcm() requires exactly one of specified_n_per_group or desired_power", {
  # Neither supplied: nothing to solve for. Without the guard this used to die
  # with an obscure length-zero error inside the sample-size search loop.
  expect_error(
    ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                 frequency = 5, duration = 4),
    "Specify one of 'specified_n_per_group' or 'desired_power'"
  )
  # Both supplied: the problem is over-determined (the unspecified one is the
  # quantity solved for), so this is rejected rather than silently ignored.
  expect_error(
    ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                 frequency = 5, duration = 4, specified_n_per_group = 50, desired_power = 0.80),
    "Specify only one of 'specified_n_per_group' or 'desired_power'"
  )
})

test_that("ss_power_pcm() solves N at 0.80 power for the NYS design (R&L 2001, Table 2)", {
  res <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003, true_var_error = .0262,
                      frequency = 1, duration = 4, desired_power = .80)
  v <- function(t) res$value[res$term == t]
  expect_equal(v("total_N"), 370)                 # 185 per group
  expect_gte(v("actual_power"), 0.80 - 1e-6)
  # Their Table 2 cells at D = 4, f = 1: n = 300 -> .71, n = 400 -> .83.
  pwr <- function(n) {
    r <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003, true_var_error = .0262,
                      frequency = 1, duration = 4, specified_n_per_group = n / 2)
    r$value[r$term == "actual_power"]
  }
  expect_equal(round(pwr(300), 2), 0.71)
  expect_equal(round(pwr(400), 2), 0.83)
})

test_that("ss_power_pcm two-sided power counts both tails (HIGH-02)", {
  # Nondirectional power must equal the base-R noncentral-t two-tail value, and
  # at a null slope (beta = 0) must equal the nominal alpha, not alpha/2.
  two_tail <- function(ncp, df, a = 0.05) {
    cv <- stats::qt(1 - a / 2, df)
    stats::pt(-cv, df, ncp) + stats::pt(cv, df, ncp, lower.tail = FALSE)
  }
  r0 <- ss_power_pcm(true_reg_coef = 0, true_var_trend = 0.05, true_var_error = 1, frequency = 5,
                     duration = 4, specified_n_per_group = 50)
  expect_equal(r0$value[r0$term == "actual_power"], 0.05, tolerance = 1e-10)

  r1 <- ss_power_pcm(true_reg_coef = 0.05, true_var_trend = 0.05, true_var_error = 1, frequency = 5,
                     duration = 4, specified_n_per_group = 50)
  ncp <- r1$value[r1$term == "noncentral_t_parm"]
  expect_equal(r1$value[r1$term == "actual_power"],
               two_tail(ncp, 100 - 2), tolerance = 1e-10)
})

test_that("ss_power_pcm errors instead of hanging on an unreachable target (HIGH-04)", {
  # A null slope has power alpha at every N, so no sample size reaches a higher
  # target; the search must fail fast with an informative error, not loop.
  expect_error(
    ss_power_pcm(true_reg_coef = 0, true_var_trend = 0.05, true_var_error = 1, frequency = 5,
                 duration = 4, desired_power = 0.8),
    "unreachable|Could not reach")
  # A reachable target still resolves.
  r <- ss_power_pcm(true_reg_coef = 0.3, true_var_trend = 0.05, true_var_error = 1, frequency = 5,
                    duration = 4, desired_power = 0.8)
  expect_gte(r$value[r$term == "actual_power"], 0.8)
})

test_that("ss_power_pcm() refuses a desired_power outside (0, 1)", {
  # No range check was made before 1.1.0; the rule and the classed
  # refusal are those of every power planner.
  for (bad in c(0, 1, 1.2)) {
    expect_error(
      ss_power_pcm(true_std_reg_coef = 0.2, true_var_trend = 0.05, true_var_error = 1,
                   frequency = 5, duration = 4, desired_power = bad),
      class = "dmar_desired_power_invalid")
  }
})

test_that("ss_power_pcm() takes the difference on one scale, true_reg_coef or true_std_reg_coef, and reports both", {
  skip_on_cran()
  # Ken's ruling of 2026-11-11: the name of the argument says the scale.
  # true_std_reg_coef is the group difference divided by sqrt(true_var_trend)
  # (Raudenbush & Liu, 2001, p. 391); the table echoes the one given and
  # holds the same difference on the other scale. DMAR 1.0.x took beta with a
  # flag, standardized, and named the rows unstd_coefficient and
  # std_coefficient.
  std <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003,
                      true_var_error = .0262, frequency = 1, duration = 4,
                      specified_n_per_group = 119)
  raw <- ss_power_pcm(true_reg_coef = -.4 * sqrt(.003), true_var_trend = .003,
                      true_var_error = .0262, frequency = 1, duration = 4,
                      specified_n_per_group = 119)
  v <- function(x, t) x$value[x$term == t]
  expect_identical(v(std, "true_std_reg_coef"), -.4)
  expect_equal(v(std, "true_reg_coef"), -.4 * sqrt(.003), tolerance = 1e-15)
  expect_identical(v(raw, "true_reg_coef"), -.4 * sqrt(.003))
  expect_equal(v(raw, "true_std_reg_coef"), -.4, tolerance = 1e-15)
  # The two scales describe one design, so the power is the same, and it is
  # the power computed here from Raudenbush and Liu's noncentrality,
  # n beta^2 [tau / (tau + V)] / 4 (their Equations 12, 14, and 15), on
  # n - 2 degrees of freedom: .612249 (their Table 1 gives .61).
  V <- .0262 * 12 / (5 * 24)
  lam <- sqrt(238 * .16 * .003 / (.003 + V) / 4)
  cv <- stats::qt(.975, 236)
  pw <- stats::pt(cv, 236, lam, lower.tail = FALSE) + stats::pt(-cv, 236, lam)
  expect_equal(v(std, "actual_power"), pw, tolerance = 1e-12)
  expect_equal(v(raw, "actual_power"), pw, tolerance = 1e-12)
  expect_equal(v(std, "noncentral_t_parm"), lam, tolerance = 1e-12)
  # Reversing the sign of the difference leaves the power unchanged.
  flip <- ss_power_pcm(true_std_reg_coef = .4, true_var_trend = .003,
                       true_var_error = .0262, frequency = 1, duration = 4,
                       specified_n_per_group = 119)
  expect_equal(v(flip, "actual_power"), pw, tolerance = 1e-12)
  # The planned size is the same on either scale.
  n_std <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003,
                        true_var_error = .0262, frequency = 1, duration = 4,
                        desired_power = .8)
  n_raw <- ss_power_pcm(true_reg_coef = -.4 * sqrt(.003), true_var_trend = .003,
                        true_var_error = .0262, frequency = 1, duration = 4,
                        desired_power = .8)
  expect_equal(v(n_std, "total_N"), 370)
  expect_equal(v(n_raw, "total_N"), 370)
  expect_false(any(c("unstd_coefficient", "std_coefficient") %in%
                     c(std$term, raw$term)))
  # The order of the rows is that of DMAR 1.0.x.
  expect_identical(std$term[8:9], c("true_reg_coef", "true_std_reg_coef"))
  expect_identical(std$term[12], "true_var_trend_error")
  expect_false("error_var_of_slopes" %in% std$term)
})

test_that("ss_power_pcm() refuses both scales, neither, and the 1.0.x flag", {
  args <- list(true_var_trend = .003, true_var_error = .0262, frequency = 1,
               duration = 4, specified_n_per_group = 119)
  e <- tryCatch(do.call(ss_power_pcm, c(list(true_reg_coef = -.02,
                                             true_std_reg_coef = -.4), args)),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, c("true_reg_coef", "true_std_reg_coef"))
  expect_match(conditionMessage(e), "both were given")
  e <- tryCatch(do.call(ss_power_pcm, args), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_match(conditionMessage(e), "neither was given")
  expect_error(do.call(ss_power_pcm, c(list(true_std_reg_coef = NA), args)),
               class = "dmar_argument_invalid")
  # The flag of DMAR 1.0.x and its argument name stop with R's own error.
  for (old in list(list(true_reg_coef = -.4, standardized = TRUE),
                   list(true_std_reg_coef = -.4, standardized = FALSE),
                   list(beta = -.4))) {
    expect_error(do.call(ss_power_pcm, c(old, args)), "unused argument")
  }
  expect_false("standardized" %in% names(formals(ss_power_pcm)))
  # A positional 1.0.x call (beta, tau, level_1_variance first) gives both
  # scales and is refused rather than read on the other scale.
  expect_error(ss_power_pcm(-.4, .003, .0262, 1, 4, specified_n_per_group = 119),
               class = "dmar_argument_invalid")
})
