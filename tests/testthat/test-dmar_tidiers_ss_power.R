# Regression tests for the dmar_ss_power broom methods on the user-supplied-N
# ("what power does this N give me?") path. tidy() / glance() previously
# returned NA for the sample size of ss_power_r() and ss_power_smd() because
# the specified_* term names were missing from the lookup vectors.

test_that("tidy() recovers the user-supplied N on the realized-power path", {
  expect_equal(generics::tidy(ss_power_r(true_cor = 0.30, specified_N = 100))$estimate, 100)
  expect_equal(generics::tidy(ss_power_smd(true_smd = 0.5, n_1 = 30))$estimate, 30)
})

test_that("the necessary-N (planning) path still reports the planned N", {
  expect_equal(generics::tidy(ss_power_r(true_cor = 0.30, desired_power = 0.80))$estimate,
               ss_power_r(true_cor = 0.30, desired_power = 0.80)$value[
                 ss_power_r(true_cor = 0.30, desired_power = 0.80)$term == "necessary_N"])
})

test_that("for an unbalanced smd the reported N is n_1", {
  expect_equal(generics::tidy(ss_power_smd(true_smd = 0.5, n_1 = 30, n_2 = 50))$estimate, 30)
})

test_that("glance() exposes power on the realized-power path without NA sample size", {
  g <- generics::glance(ss_power_r(true_cor = 0.30, specified_N = 100))
  expect_true(is.data.frame(g))
  expect_false(anyNA(generics::tidy(ss_power_r(true_cor = 0.30, specified_N = 100))$estimate))
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
  pc <- ss_power_c(true_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                   true_sd = 1, desired_power = 0.80)
  expect_s3_class(pc, "dmar_ss_power")
  expect_equal(generics::tidy(pc)$estimate,
               pc$value[pc$term == "necessary_n_per_group"])
  expect_equal(generics::tidy(pc)$power,
               pc$value[pc$term == "actual_power"])

  # ss_power_equivalence_c carries necessary_n_per_group, N, actual_power;
  # tidy() must report the per-group size, never the total N (= J n).
  pe <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                               contrast_upper = 5, desired_power = 0.90)
  expect_s3_class(pe, "dmar_ss_power")
  per_group <- pe$value[pe$term == "necessary_n_per_group"]
  expect_equal(generics::tidy(pe)$estimate, per_group)
  expect_false(isTRUE(all.equal(generics::tidy(pe)$estimate,
                                pe$value[pe$term == "total_N"])))
  expect_equal(generics::glance(pe)$estimate, per_group)
})

test_that("ss_power_c_ancova and ss_power_sc summarize the per-group size", {
  for (res in list(
    ss_power_c_ancova(true_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                      true_sd = 1, true_cor = 0.5, desired_power = 0.80),
    ss_power_sc(true_std_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
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
  plan <- ss_power_one_way_anova(a = 3, true_f = 0.25, desired_power = 0.80)
  expect_s3_class(plan, "dmar_ss_power")
  expect_equal(generics::tidy(plan)$estimate,
               plan$value[plan$term == "necessary_N"])
  expect_equal(generics::tidy(plan)$power,
               plan$value[plan$term == "actual_power"])
  # Specified-N branch reports specified_N.
  got <- ss_power_one_way_anova(a = 3, true_f = 0.25, specified_N = 150)
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
                                        effect_indices = 2, true_f = 0.25,
                                        desired_power = 0.80),
         size = "necessary_n_per_cell"),
    list(res = ss_power_factorial_ancova(factor_levels = c(2, 4, 3),
                                         effect_indices = 1, true_f = 0.10,
                                         true_covariate_R2 = 0, n_covariates = 0,
                                         desired_power = 0.80),
         size = "necessary_n_per_cell"),
    list(res = ss_power_split_plot_anova(a = 2, b = 4, effect = "between",
                                    true_f = 0.25, true_cor = 0.5, desired_power = 0.80),
         size = "necessary_n_per_group"),
    list(res = ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5,
                                 desired_power = 0.80),
         size = "necessary_N"),
    list(res = ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05,
                                      desired_power = 0.80),
         size = "necessary_n_clusters_per_arm"),
    list(res = ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                                 true_means = c(90, 92, 88, 81), true_var_error = 144,
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
  r1 <- ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, specified_N = 20)
  expect_equal(generics::tidy(r1)$estimate, r1$value[r1$term == "specified_N"])
  r2 <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05, specified_n_clusters_per_arm = 25)
  expect_equal(generics::tidy(r2)$estimate, r2$value[r2$term == "specified_n_clusters_per_arm"])
  r3 <- ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1), true_means = c(90, 92, 88, 81),
                          true_var_error = 144, specified_n_per_group = 20)
  expect_equal(generics::tidy(r3)$estimate, r3$value[r3$term == "specified_n_per_group"])
})

test_that("ss_power_rc forwards ss_power_reg_coef's tidy tagged output", {
  r <- ss_power_rc(alpha_level = .05, true_f2 = 0.2130898, p = 5,
                   directional = FALSE, desired_power = .80)
  expect_s3_class(r, "dmar_ss_power")
  expect_equal(generics::tidy(r)$estimate, r$value[r$term == "necessary_N"])
  expect_equal(generics::tidy(r)$power, r$value[r$term == "actual_power"])
  # Realized-power branch reports the specified total.
  rp <- ss_power_rc(alpha_level = .05, true_f2 = 0.2130898, p = 5,
                    specified_N = 25, directional = FALSE)
  expect_equal(generics::tidy(rp)$estimate, rp$value[rp$term == "specified_N"])
})

test_that("ss_power_indirect_effect summarizes N and the joint indirect power", {
  r <- ss_power_indirect_effect(true_std_reg_coef_a = .39, true_std_reg_coef_b = .39, desired_power = .80)
  expect_s3_class(r, "dmar_ss_power")
  # tidy() reports the sample size and the joint (indirect-effect) power,
  # the row power_ab.
  expect_equal(generics::tidy(r)$estimate, r$value[r$term == "necessary_N"])
  expect_equal(generics::tidy(r)$power, r$value[r$term == "power_ab"])
  # glance() carries the diagnostic component powers as extra columns.
  g <- generics::glance(r)
  expect_true(all(c("power_a", "power_b") %in% names(g)))
  expect_equal(g$power_a, r$value[r$term == "power_a"])
})

test_that("glance() leaves out only the power row tidy() reads, so a composite table keeps power_ab", {
  # power_ab joined the power registry for ss_power_indirect_effect(), and
  # glance() then left out every registry row, so a composite table with
  # parameters labeled a, b, and ab lost its marginal power power_ab.
  x <- data.frame(term = c("specified_N", "composite_power", "power_a",
                           "power_b", "power_ab", "alpha_level"),
                  value = c(100, .5, .9, .8, .6, .05))
  class(x) <- c("dmar_ss_power", "data.frame")
  expect_equal(generics::tidy(x)$power, .5)
  g <- generics::glance(x)
  expect_identical(names(g), c("term", "estimate", "power", "power_a", "power_b",
                               "power_ab", "alpha_level"))
  expect_equal(c(g$power, g$power_a, g$power_b, g$power_ab), c(.5, .9, .8, .6))
  # ss_power_indirect_effect(): power_ab is the design's power, read by
  # tidy() and reported by glance() as `power`, not again as a column.
  r <- ss_power_indirect_effect(true_std_reg_coef_a = .3, true_std_reg_coef_b = .3,
                                specified_N = 100)
  gi <- generics::glance(r)
  expect_false("power_ab" %in% names(gi))
  expect_equal(gi$power, r$value[r$term == "power_ab"])
  expect_true(all(c("power_a", "power_b") %in% names(gi)))
  # A table whose power row is one of the other registry names keeps its
  # remaining rows, as before.
  y <- data.frame(term = c("necessary_n_per_group", "total_N", "actual_power",
                           "alpha_level"),
                  value = c(50, 100, .81, .05))
  class(y) <- c("dmar_ss_power", "data.frame")
  expect_identical(names(generics::glance(y)),
                   c("term", "estimate", "power", "total_N", "alpha_level"))
})

test_that("glance() of ss_power_composite_sem() with parameters a, b, and ab keeps power_ab", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  pop <- "
    x ~~ 1*x
    m ~ 0.4*x
    m ~~ 0.84*m
    y ~ 0.35*m + 0.15*x
    y ~~ 0.813*y
  "
  mod <- "
    m ~ a*x
    y ~ b*m + cp*x
    ab := a*b
  "
  r <- ss_power_composite_sem(model = mod, true_model = pop,
                              parameters = c("a", "b", "ab"),
                              specified_N = 100, G = 25, seed = 113)
  expect_true(all(c("composite_power", "power_a", "power_b", "power_ab") %in% r$term))
  expect_equal(generics::tidy(r)$power, r$value[r$term == "composite_power"])
  g <- generics::glance(r)
  expect_true(all(c("power_a", "power_b", "power_ab") %in% names(g)))
  expect_equal(g$power_ab, r$value[r$term == "power_ab"])
  expect_equal(g$power, r$value[r$term == "composite_power"])
})

test_that("the sensitivity siblings report empirical vs analytic power", {
  skip_on_cran()
  set.seed(113)
  s <- ss_power_R2_sensitivity(true_R2 = 0.30, supposed_R2 = 0.30,
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
  ct <- generics::tidy(ss_power_contrast(c_weights = c(1, -1, 0), true_var_error = 1,
                                         true_contrast = 0.5, specified_n_per_group = c(20, 40, 30)))
  expect_false(is.na(ct$estimate))
  expect_equal(ct$estimate, 90)          # total N
  expect_false(is.na(ct$power))

  # Unbalanced SMD: glance() must retain the second group's size.
  g <- generics::glance(ss_power_smd(true_smd = 0.5, n_1 = 30, n_2 = 50))
  expect_true("specified_n_2" %in% names(g))
  expect_equal(g$specified_n_2, 50)
})
