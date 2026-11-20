## Smoke + sanity tests for the new MDK power-planning functions.

# ---- ss_power_smd ----

test_that("ss_power_smd() necessary_n is monotonic in desired_power", {
  small <- ss_power_smd(true_smd = 0.5, desired_power = 0.70)
  large <- ss_power_smd(true_smd = 0.5, desired_power = 0.95)
  expect_lt(small[small$term == "necessary_n_per_group", "value"],
            large[large$term == "necessary_n_per_group", "value"])
})

test_that("ss_power_smd() necessary_n is monotonic in effect size", {
  big_eff   <- ss_power_smd(true_smd = 0.8, desired_power = 0.80)
  small_eff <- ss_power_smd(true_smd = 0.2, desired_power = 0.80)
  expect_lt(big_eff[big_eff$term == "necessary_n_per_group", "value"],
            small_eff[small_eff$term == "necessary_n_per_group", "value"])
})

test_that("ss_power_smd() power at returned n equals or exceeds desired_power", {
  res <- ss_power_smd(true_smd = 0.5, desired_power = 0.80)
  n   <- res[res$term == "necessary_n_per_group", "value"]
  pwr <- ss_power_smd(true_smd = 0.5, n_1 = n)[3, 2]
  expect_gte(pwr, 0.80)
})

test_that("ss_power_smd() validates inputs", {
  expect_error(ss_power_smd(true_smd = NA, desired_power = 0.8),    "'true_smd' must be a single finite numeric")
  expect_error(ss_power_smd(true_smd = 0.5, alpha_level = 1.5),     "'alpha_level' must be a single number in")
})

# ---- ss_power_one_way_anova ----

test_that("ss_power_one_way_anova() necessary_N decreases with larger f", {
  small <- ss_power_one_way_anova(a = 3, true_f = 0.10, desired_power = 0.80)
  large <- ss_power_one_way_anova(a = 3, true_f = 0.40, desired_power = 0.80)
  expect_lt(large[large$term == "necessary_N", "value"],
            small[small$term == "necessary_N", "value"])
})

test_that("ss_power_one_way_anova() f and eta_squared give the same answer", {
  eta2 <- 0.0588  # corresponds to f = 0.25
  f_res    <- ss_power_one_way_anova(a = 3, true_f = sqrt(eta2 / (1 - eta2)), desired_power = 0.80)
  eta2_res <- ss_power_one_way_anova(a = 3, true_eta_squared = eta2, desired_power = 0.80)
  expect_equal(f_res$value, eta2_res$value)
})

test_that("ss_power_one_way_anova() errors when neither true_f nor true_eta_squared given", {
  expect_error(ss_power_one_way_anova(a = 3, desired_power = 0.80),
               "Specify exactly one of 'true_f' or 'true_eta_squared'")
})

# ---- ss_power_c and ss_power_sc ----

test_that("ss_power_c() with true_sd = 1 and true_contrast equal to true_std_contrast matches ss_power_sc()", {
  psi_val <- 0.5
  cw <- c(0.5, 0.5, -0.5, -0.5)
  c_res  <- ss_power_c(true_contrast = psi_val, c_weights = cw, true_sd = 1, desired_power = 0.80)
  sc_res <- ss_power_sc(true_std_contrast = psi_val, c_weights = cw, desired_power = 0.80)
  expect_equal(c_res$value, sc_res$value)
})

test_that("ss_power_sc() necessary_n decreases with larger effect", {
  small <- ss_power_sc(true_std_contrast = 0.2, c_weights = c(0.5, 0.5, -0.5, -0.5), desired_power = 0.80)
  large <- ss_power_sc(true_std_contrast = 0.8, c_weights = c(0.5, 0.5, -0.5, -0.5), desired_power = 0.80)
  expect_lt(large[large$term == "necessary_n_per_group", "value"],
            small[small$term == "necessary_n_per_group", "value"])
})

test_that("ss_power_c() rejects c_weights that do not sum to zero", {
  expect_error(ss_power_c(true_contrast = 0.5, c_weights = c(1, 1, 0, 0), true_sd = 1, desired_power = 0.80),
               "must sum to zero")
})

# ---- ss_power_c_ancova ----

test_that("ss_power_c_ancova() needs fewer subjects than ss_power_c when true_cor > 0", {
  cw <- c(0.5, 0.5, -0.5, -0.5)
  anova_res  <- ss_power_c(       true_contrast = 0.5, c_weights = cw, true_sd = 1,            desired_power = 0.80)
  ancova_res <- ss_power_c_ancova(true_contrast = 0.5, c_weights = cw, true_sd = 1, true_cor = 0.5, desired_power = 0.80)
  expect_lt(ancova_res[ancova_res$term == "necessary_n_per_group", "value"],
            anova_res[anova_res$term  == "necessary_n_per_group", "value"])
})

test_that("ss_power_c_ancova() with true_cor = 0 closely matches ss_power_c() (off by one due to df adjustment)", {
  cw <- c(0.5, 0.5, -0.5, -0.5)
  c_res      <- ss_power_c(       true_contrast = 0.5, c_weights = cw, true_sd = 1,            desired_power = 0.80)
  ancova_res <- ss_power_c_ancova(true_contrast = 0.5, c_weights = cw, true_sd = 1, true_cor = 0,   desired_power = 0.80)
  diff <- abs(ancova_res[ancova_res$term == "necessary_n_per_group", "value"] -
              c_res[c_res$term == "necessary_n_per_group", "value"])
  expect_lte(diff, 1)
})

test_that("ss_power_c_ancova() rejects true_cor outside (-1, 1)", {
  expect_error(ss_power_c_ancova(true_contrast = 0.5, c_weights = c(.5, -.5), true_sd = 1, true_cor = 1.0, desired_power = 0.80),
               "'true_cor' must be a single numeric value in")
})

# ---- ss_power_factorial_anova ----

test_that("ss_power_factorial_anova() with one-factor design matches ss_power_one_way_anova()", {
  one_way   <- ss_power_one_way_anova(a = 3, true_f = 0.25, desired_power = 0.80)
  factorial <- ss_power_factorial_anova(factor_levels = 3, effect_indices = 1, true_f = 0.25, desired_power = 0.80)
  expect_equal(
    one_way[one_way$term == "necessary_N",   "value"],
    factorial[factorial$term == "total_N", "value"]
  )
})

test_that("ss_power_factorial_anova() interaction df_effect is product of (levels - 1)", {
  res <- ss_power_factorial_anova(factor_levels = c(2, 3, 4), effect_indices = c(1, 2, 3),
                                  true_f = 0.20, desired_power = 0.80)
  expect_equal(res[res$term == "df_effect", "value"], (2 - 1) * (3 - 1) * (4 - 1))
})

test_that("ss_power_factorial_anova() rejects bad indices", {
  expect_error(ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = c(1, 1),
                                        true_f = 0.25, desired_power = 0.80),
               "must be unique integer indices")
  expect_error(ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 5,
                                        true_f = 0.25, desired_power = 0.80),
               "must be unique integer indices")
})

# ---- ss_power_rm_anova ----

test_that("ss_power_rm_anova() power increases with larger true_cor (positive within-subject correlation)", {
  res_low  <- ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.0, desired_power = 0.80)
  res_high <- ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.6, desired_power = 0.80)
  expect_lt(res_high[res_high$term == "necessary_N", "value"],
            res_low[res_low$term  == "necessary_N", "value"])
})

test_that("ss_power_rm_anova() Greenhouse-Geisser epsilon < 1 increases necessary n", {
  no_eps  <- ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, true_epsilon = 1.0,  desired_power = 0.80)
  with_eps <- ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, true_epsilon = 0.6, desired_power = 0.80)
  expect_gte(with_eps[with_eps$term == "necessary_N", "value"],
             no_eps[no_eps$term      == "necessary_N", "value"])
})

# ---- ss_power_r ----

test_that("ss_power_r() necessary_N decreases with larger rho", {
  small <- ss_power_r(true_cor = 0.10, desired_power = 0.80)
  large <- ss_power_r(true_cor = 0.50, desired_power = 0.80)
  expect_lt(large[large$term == "necessary_N", "value"],
            small[small$term == "necessary_N", "value"])
})

test_that("ss_power_r() power at returned N equals or exceeds desired_power", {
  res <- ss_power_r(true_cor = 0.30, desired_power = 0.80)
  N   <- res[res$term == "necessary_N", "value"]
  pwr <- ss_power_r(true_cor = 0.30, specified_N = N)[2, 2]
  expect_gte(pwr, 0.80)
})

test_that("ss_power_r() power at N - 1 is below desired_power (smallest sufficient N)", {
  res <- ss_power_r(true_cor = 0.30, desired_power = 0.80)
  N   <- res[res$term == "necessary_N", "value"]
  pwr_below <- ss_power_r(true_cor = 0.30, specified_N = N - 1)[2, 2]
  expect_lt(pwr_below, 0.80)
})

test_that("ss_power_r() rejects degenerate rho == rho_0", {
  expect_error(ss_power_r(true_cor = 0.30, null_cor = 0.30, desired_power = 0.80),
               "no effect to detect")
})
