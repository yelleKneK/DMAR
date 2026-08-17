# Regression tests for the dmar_ss_power broom methods on the user-supplied-N
# ("what power does this N give me?") path. tidy() / glance() previously
# returned NA for the sample size of ss_power_r() and ss_power_smd() because
# the specified_* term names were missing from the lookup vectors.

test_that("tidy() recovers the user-supplied N on the realized-power path", {
  expect_equal(generics::tidy(ss_power_r(rho = 0.30, N = 100))$estimate, 100)
  expect_equal(generics::tidy(ss_power_smd(smd = 0.5, n_1 = 30))$estimate, 30)
})

test_that("the necessary-N (planning) path still reports the planned N", {
  expect_equal(generics::tidy(ss_power_r(rho = 0.30, desired_power = 0.80))$estimate,
               ss_power_r(rho = 0.30, desired_power = 0.80)$value[
                 ss_power_r(rho = 0.30, desired_power = 0.80)$term == "necessary_N"])
})

test_that("for an unbalanced smd the reported N is n_1", {
  expect_equal(generics::tidy(ss_power_smd(smd = 0.5, n_1 = 30, n_2 = 50))$estimate, 30)
})

test_that("glance() exposes power on the realized-power path without NA sample size", {
  g <- generics::glance(ss_power_r(rho = 0.30, N = 100))
  expect_true(is.data.frame(g))
  expect_false(anyNA(generics::tidy(ss_power_r(rho = 0.30, N = 100))$estimate))
})

test_that("the sample-size and power term vectors are the single source of truth", {
  # glance() must exclude exactly the rows tidy() already reports, or a row
  # would be both summarized and appended as though it were a planning input.
  # The two vectors drifting apart is how the sample-size list came to hold
  # necessary_n_per_group but not specified_n_per_group.
  expect_true(all(c("necessary_n_per_group", "specified_n_per_group") %in%
                    DMAR:::.SS_POWER_SIZE_TERMS))
  expect_true("composite_power" %in% DMAR:::.SS_POWER_POWER_TERMS)
})

test_that("a per-group row is preferred to its total-N counterpart", {
  # Priority order is the order of .SS_POWER_SIZE_TERMS, and it must resolve the
  # same way whether the size was planned or specified: tidy() cannot mean
  # per-group in one branch and total in the other for the same design.
  planned <- data.frame(term = c("necessary_n_per_group", "necessary_N"),
                        value = c(95, 190))
  specified <- data.frame(term = c("specified_n_per_group", "specified_N"),
                          value = c(95, 190))
  expect_equal(DMAR:::.ss_power_sample_size(planned), 95)
  expect_equal(DMAR:::.ss_power_sample_size(specified), 95)
})

test_that("composite_power is recognized as the power of a design", {
  x <- data.frame(term = c("specified_n_per_group", "composite_power"),
                  value = c(95, 0.8))
  expect_equal(DMAR:::.ss_power_actual_power(x), 0.8)
})

test_that("a table with no recognized size or power row yields NA, not an error", {
  x <- data.frame(term = c("something", "else"), value = c(1, 2))
  expect_true(is.na(DMAR:::.ss_power_sample_size(x)))
  expect_true(is.na(DMAR:::.ss_power_actual_power(x)))
})

test_that("the contrast planners ss_power_c and ss_power_equivalence_c summarize", {
  # ss_power_c reports one per-group size and one power, with term names already
  # in the lookup vectors, so tidy() resolves the per-group size directly.
  pc <- ss_power_c(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                   sigma = 1, desired_power = 0.80)
  expect_s3_class(pc, "dmar_ss_power")
  expect_equal(generics::tidy(pc)$estimate,
               pc$value[pc$term == "necessary_n_per_group"])
  expect_equal(generics::tidy(pc)$power,
               pc$value[pc$term == "actual_power"])

  # ss_power_equivalence_c carries necessary_n_per_group, N, actual_power;
  # tidy() must report the per-group size, never the total N (= J n).
  pe <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                               delta_upper = 5, desired_power = 0.90)
  expect_s3_class(pe, "dmar_ss_power")
  per_group <- pe$value[pe$term == "necessary_n_per_group"]
  expect_equal(generics::tidy(pe)$estimate, per_group)
  expect_false(isTRUE(all.equal(generics::tidy(pe)$estimate,
                                pe$value[pe$term == "total_N"])))
  expect_equal(generics::glance(pe)$estimate, per_group)
})

test_that("ss_power_c_ancova and ss_power_sc summarize the per-group size", {
  for (res in list(
    ss_power_c_ancova(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                      sigma = 1, rho = 0.5, desired_power = 0.80),
    ss_power_sc(psi_standardized = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                desired_power = 0.80))) {
    expect_s3_class(res, "dmar_ss_power")
    expect_equal(generics::tidy(res)$estimate,
                 res$value[res$term == "necessary_n_per_group"])
    expect_equal(generics::tidy(res)$power,
                 res$value[res$term == "actual_power"])
  }
})

test_that("ss_power_one_way_anova summarizes the total N in both branches", {
  # Planning branch: rows include both necessary_N and n_per_group; tidy() must
  # report the total (necessary_N), consistently with the specified-N branch.
  plan <- ss_power_one_way_anova(a = 3, f = 0.25, desired_power = 0.80)
  expect_s3_class(plan, "dmar_ss_power")
  expect_equal(generics::tidy(plan)$estimate,
               plan$value[plan$term == "necessary_N"])
  expect_equal(generics::tidy(plan)$power,
               plan$value[plan$term == "actual_power"])
  # Specified-N branch reports specified_N.
  got <- ss_power_one_way_anova(a = 3, f = 0.25, N = 150)
  expect_equal(generics::tidy(got)$estimate,
               got$value[got$term == "specified_N"])
})

test_that("the ANOVA and cluster planners summarize their planning-unit size", {
  # Each planner's tidy() estimate is the design's planning-unit size, taken
  # from the row named below; the power is actual_power. The size-row names are
  # design-specific (per cell, per subject, per cluster) and were added to
  # .SS_POWER_SIZE_TERMS so the broom lookup resolves them ahead of the total.
  cases <- list(
    list(res = ss_power_factorial_anova(factor_levels = c(2, 3),
                                        effect_indices = 2, f = 0.25,
                                        desired_power = 0.80),
         size = "necessary_n_per_cell"),
    list(res = ss_power_factorial_ancova(factor_levels = c(2, 4, 3),
                                         effect_indices = 1, f = 0.10,
                                         covariate_R2 = 0, n_covariates = 0,
                                         desired_power = 0.80),
         size = "necessary_n_per_cell"),
    list(res = ss_power_split_plot_anova(a = 2, b = 4, effect = "between",
                                    f = 0.25, rho = 0.5, desired_power = 0.80),
         size = "necessary_n_per_group"),
    list(res = ss_power_rm_anova(a = 4, f = 0.25, rho = 0.5,
                                 desired_power = 0.80),
         size = "necessary_n_subjects"),
    list(res = ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05,
                                      desired_power = 0.80),
         size = "necessary_J_per_arm"),
    list(res = ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                                 mu = c(90, 92, 88, 81), sigma_squared = 144,
                                 desired_power = 0.90),
         size = "necessary_n_per_group")
  )
  for (case in cases) {
    expect_s3_class(case$res, "dmar_ss_power")
    expect_equal(generics::tidy(case$res)$estimate,
                 case$res$value[case$res$term == case$size])
    expect_equal(generics::tidy(case$res)$power,
                 case$res$value[case$res$term == "actual_power"])
    expect_equal(nrow(generics::glance(case$res)), 1L)
  }
})

test_that("the realized-power branch reports the specified planning-unit size", {
  r1 <- ss_power_rm_anova(a = 4, f = 0.25, rho = 0.5, n = 20)
  expect_equal(generics::tidy(r1)$estimate, r1$value[r1$term == "specified_n_subjects"])
  r2 <- ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, J = 25)
  expect_equal(generics::tidy(r2)$estimate, r2$value[r2$term == "specified_J_per_arm"])
  r3 <- ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1), mu = c(90, 92, 88, 81),
                          sigma_squared = 144, n_per_group = 20)
  expect_equal(generics::tidy(r3)$estimate, r3$value[r3$term == "specified_n_per_group"])
})

test_that("ss_power_rc forwards ss_power_reg_coef's tidy tagged output", {
  r <- ss_power_rc(alpha_level = .05, cohen_f2 = 0.2130898, p = 5,
                   directional = FALSE, desired_power = .80)
  expect_s3_class(r, "dmar_ss_power")
  expect_equal(generics::tidy(r)$estimate, r$value[r$term == "necessary_N"])
  expect_equal(generics::tidy(r)$power, r$value[r$term == "actual_power"])
  # Realized-power branch reports the specified total.
  rp <- ss_power_rc(alpha_level = .05, cohen_f2 = 0.2130898, p = 5,
                    specified_N = 25, directional = FALSE)
  expect_equal(generics::tidy(rp)$estimate, rp$value[rp$term == "specified_N"])
})

test_that("ss_power_indirect_effect summarizes N and the joint indirect power", {
  r <- ss_power_indirect_effect(a = .39, b = .39, desired_power = .80)
  expect_s3_class(r, "dmar_ss_power")
  # tidy() reports the sample size and the joint (indirect-effect) power.
  expect_equal(generics::tidy(r)$estimate, r$value[r$term == "necessary_N"])
  expect_equal(generics::tidy(r)$power, r$value[r$term == "actual_power"])
  # glance() carries the diagnostic component powers as extra columns.
  g <- generics::glance(r)
  expect_true(all(c("power_a", "power_b") %in% names(g)))
  expect_equal(g$power_a, r$value[r$term == "power_a"])
})

test_that("the sensitivity siblings report empirical vs analytic power", {
  skip_on_cran()
  set.seed(113)
  s <- ss_power_R2_sensitivity(true_R2 = 0.30, estimated_R2 = 0.30,
                               desired_power = 0.80, p = 5,
                               random_predictors = FALSE,
                               generate_random_predictors = TRUE,
                               G = 200, print_iter = FALSE)
  expect_s3_class(s, "dmar_ss_power_sensitivity")
  td <- generics::tidy(s)
  expect_equal(td$sample_size, s$value[s$term == "total_N"])
  expect_equal(td$empirical_power, s$value[s$term == "empirical_power"])
  expect_equal(td$analytic_power, s$value[s$term == "analytic_power"])
  g <- generics::glance(s)
  expect_equal(nrow(g), 1L)
  expect_true("mean_R2" %in% names(g))       # simulated distribution carried
  expect_false("total_N" %in% names(g))    # head rows not duplicated as inputs
})

test_that("unbalanced designs keep a non-missing size and retain every size (MEDIUM-01)", {
  # Unbalanced contrast: the single per-group size is NA, so tidy() must fall
  # through to the total N rather than reporting estimate = NA.
  ct <- generics::tidy(ss_power_contrast(c_weights = c(1, -1, 0), sigma_squared = 1,
                                         psi = 0.5, n_per_group = c(20, 40, 30)))
  expect_false(is.na(ct$estimate))
  expect_equal(ct$estimate, 90)          # total N
  expect_false(is.na(ct$power))

  # Unbalanced SMD: glance() must retain the second group's size.
  g <- generics::glance(ss_power_smd(smd = 0.5, n_1 = 30, n_2 = 50))
  expect_true("specified_n_2" %in% names(g))
  expect_equal(g$specified_n_2, 50)
})
