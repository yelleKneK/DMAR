test_that("ss_power_pcm() returns a tidy data.frame with the documented columns", {
  res <- ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_per_group", "necessary_n_per_group", "total_N", "actual_power",
                    "measurement_occasions") %in% res$term))
})

test_that("ss_power_pcm() actual_power at returned n meets desired_power", {
  res <- ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  expect_gte(res$value[res$term == "actual_power"], 0.80 - 1e-6)
})

test_that("ss_power_pcm() total_N = ss_c + ss_t", {
  res <- ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  ss_c <- res$value[res$term == "necessary_n_per_group"]
  ss_t <- res$value[res$term == "necessary_n_per_group"]
  expect_equal(res$value[res$term == "total_N"], ss_c + ss_t)
})

test_that("ss_power_pcm() measurement_occasions = frequency * duration + 1", {
  res <- ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                      frequency = 5, duration = 4, desired_power = 0.80)
  expect_equal(res$value[res$term == "measurement_occasions"], 5 * 4 + 1)
})

test_that("ss_power_pcm() necessary sample size grows as beta shrinks", {
  big_eff   <- ss_power_pcm(beta = 0.4, tau = 0.05, level_1_variance = 1,
                            frequency = 5, duration = 4,
                            desired_power = 0.80)$value[2]
  small_eff <- ss_power_pcm(beta = 0.1, tau = 0.05, level_1_variance = 1,
                            frequency = 5, duration = 4,
                            desired_power = 0.80)$value[2]
  expect_gt(small_eff, big_eff)
})

test_that("ss_power_pcm() necessary sample size grows as duration shrinks (less follow-up)", {
  long  <- ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                        frequency = 5, duration = 4,
                        desired_power = 0.80)$value[2]
  short <- ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                        frequency = 5, duration = 1,
                        desired_power = 0.80)$value[2]
  expect_gt(short, long)
})

test_that("ss_power_pcm() requires exactly one of N or desired_power", {
  # Neither supplied: nothing to solve for. Without the guard this used to die
  # with an obscure length-zero error inside the sample-size search loop.
  expect_error(
    ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                 frequency = 5, duration = 4),
    "Specify one of 'N' or 'desired_power'"
  )
  # Both supplied: the problem is over-determined (the unspecified one is the
  # quantity solved for), so this is rejected rather than silently ignored.
  expect_error(
    ss_power_pcm(beta = 0.2, tau = 0.05, level_1_variance = 1,
                 frequency = 5, duration = 4, N = 100, desired_power = 0.80),
    "Specify only one of 'N' or 'desired_power'"
  )
})

test_that("ss_power_pcm() reproduces the Raudenbush & Liu (2001) NYS benchmark", {
  # National Youth Survey example (Raudenbush & Liu, 2001, p. 393, Table 1):
  # frequency = 1, duration = 4 (M = 5), level-1 variance 0.0262, slope
  # variance tau = 0.003, standardized slope difference -0.40, n = 238. The
  # paper's hand computation gives sampling variance V = 0.00262, slope
  # reliability 0.53, noncentrality 5.046, and power 0.61. DMAR matches
  # MBESS::ss.power.pcm to six decimals (Actual.Power = 0.612236).
  res <- ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
                      frequency = 1, duration = 4, N = 238)
  v <- function(t) res$value[res$term == t]
  expect_equal(round(v("error_var_of_slopes"), 5), 0.00262)
  expect_equal(round(v("reliability"), 2), 0.53)
  expect_equal(round(v("actual_power"), 2), 0.61)
  expect_equal(v("total_N"), 238)
  expect_equal(v("measurement_occasions"), 5)
})

test_that("ss_power_pcm() solves N at 0.80 power for the NYS design (R&L 2001, Table 2)", {
  res <- ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
                      frequency = 1, duration = 4, desired_power = .80)
  v <- function(t) res$value[res$term == t]
  expect_equal(v("total_N"), 370)                 # 185 per group
  expect_gte(v("actual_power"), 0.80 - 1e-6)
  # Their Table 2 cells at D = 4, f = 1: n = 300 -> .71, n = 400 -> .83.
  pwr <- function(n) {
    r <- ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
                      frequency = 1, duration = 4, N = n)
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
  r0 <- ss_power_pcm(beta = 0, tau = 0.05, level_1_variance = 1, frequency = 5,
                     duration = 4, N = 100, standardized = FALSE)
  expect_equal(r0$value[r0$term == "actual_power"], 0.05, tolerance = 1e-10)

  r1 <- ss_power_pcm(beta = 0.05, tau = 0.05, level_1_variance = 1, frequency = 5,
                     duration = 4, N = 100, standardized = FALSE)
  ncp <- r1$value[r1$term == "noncentral_t_parm"]
  expect_equal(r1$value[r1$term == "actual_power"],
               two_tail(ncp, 100 - 2), tolerance = 1e-10)
})

test_that("ss_power_pcm errors instead of hanging on an unreachable target (HIGH-04)", {
  # A null slope has power alpha at every N, so no sample size reaches a higher
  # target; the search must fail fast with an informative error, not loop.
  expect_error(
    ss_power_pcm(beta = 0, tau = 0.05, level_1_variance = 1, frequency = 5,
                 duration = 4, desired_power = 0.8, standardized = FALSE),
    "unreachable|Could not reach")
  # A reachable target still resolves.
  r <- ss_power_pcm(beta = 0.3, tau = 0.05, level_1_variance = 1, frequency = 5,
                    duration = 4, desired_power = 0.8, standardized = FALSE)
  expect_gte(r$value[r$term == "actual_power"], 0.8)
})
