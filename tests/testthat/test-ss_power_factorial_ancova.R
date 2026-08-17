v <- function(tab, t) tab$value[tab$term == t]

test_that("ss_power_factorial_ancova() with R2 = 0 and q = 0 matches the ANOVA planner", {
  anc <- ss_power_factorial_ancova(c(2, 4, 3), effect_indices = 1, f = 0.10,
                                   covariate_R2 = 0, n_covariates = 0,
                                   desired_power = 0.80)
  ano <- ss_power_factorial_anova(c(2, 4, 3), effect_indices = 1, f = 0.10,
                                  desired_power = 0.80)
  expect_equal(v(anc, "necessary_n_per_cell"),
               ano$value[ano$term == "necessary_n_per_cell"])
  expect_equal(v(anc, "actual_power"),
               ano$value[ano$term == "actual_power"], tolerance = 1e-10)
})

test_that("covariates reduce the required sample size as the textbook says", {
  plain <- ss_power_factorial_ancova(c(2, 4, 3), 1, f = 0.10,
                                     desired_power = 0.80)
  r10 <- ss_power_factorial_ancova(c(2, 4, 3), 1, f = 0.10,
                                   covariate_R2 = 0.10, n_covariates = 2,
                                   desired_power = 0.80)
  r25 <- ss_power_factorial_ancova(c(2, 4, 3), 1, f = 0.10,
                                   covariate_R2 = 0.25, n_covariates = 2,
                                   desired_power = 0.80)
  expect_gt(v(plain, "total_N"), v(r10, "total_N"))
  expect_gt(v(r10, "total_N"), v(r25, "total_N"))
  # The adjusted effect size is f / sqrt(1 - R2).
  expect_equal(v(r25, "f_adjusted"), 0.10 / sqrt(0.75))
  # And the error df accounting: N - cells - q.
  expect_equal(v(r25, "df_error"), v(r25, "total_N") - 24 - 2)
})

test_that("realized power at a given n matches a direct noncentral F computation", {
  res <- ss_power_factorial_ancova(c(2, 4, 3), c(1, 2, 3), f = 0.15,
                                   covariate_R2 = 0.25, n_covariates = 2,
                                   n_per_cell = 6)
  N <- 6 * 24
  lambda <- N * (0.15 / sqrt(0.75))^2
  df_h <- (2 - 1) * (4 - 1) * (3 - 1)
  df_e <- N - 24 - 2
  expect_equal(v(res, "df_effect"), df_h)
  expect_equal(v(res, "actual_power"),
               1 - pf(qf(0.95, df_h, df_e), df_h, df_e, ncp = lambda))
})

test_that("ss_power_factorial_ancova() validates its arguments", {
  expect_error(ss_power_factorial_ancova(c(2, 4, 3), 1, f = 0.1,
                                         covariate_R2 = 0.2),
               "how many")
  expect_error(ss_power_factorial_ancova(c(2, 4, 3), 1, f = 0.1,
                                         covariate_R2 = 1),
               "\\[0, 1\\)")
  expect_error(ss_power_factorial_ancova(c(2, 4, 3), 4, f = 0.1),
               "indices into")
  expect_error(ss_power_factorial_ancova(c(2, 4, 3), 1), "exactly one")
  expect_error(ss_power_factorial_ancova(c(2, 4, 3), 1, f = 0.1,
                                         partial_eta_squared = 0.02),
               "exactly one")
})

test_that("ss_power_factorial_ancova errors instead of hanging on a tiny effect (HIGH-04)", {
  expect_error(
    ss_power_factorial_ancova(factor_levels = c(2, 4, 3), effect_indices = 1,
                              f = 1e-12, desired_power = 0.80),
    "too small|Could not reach|did not converge")
  r <- ss_power_factorial_ancova(factor_levels = c(2, 4, 3), effect_indices = 1,
                                 f = 0.25, desired_power = 0.80)
  expect_gte(r$value[r$term == "actual_power"], 0.80)
})
