# The 1.1.0 names of the contrast, analysis of variance, and analysis of
# covariance functions: each function takes its one name for a quantity,
# and the names it had before, those of DMAR 1.0.x and of the earlier
# 1.1.0 builds, stop. The old names are given to do.call() as data, so that
# a later rename sweep does not change them.

.renamed <- list(
  ci_c = list(new = "contrast", old = "psi"),
  ci_sc_ancova = list(new = "contrast", old = "psi"),
  equivalence_c = list(new = c("contrast", "contrast_lower", "contrast_upper"),
                       old = c("psi_hat", "psi", "delta_lower", "delta_upper")),
  power_equivalence_c = list(
    new = c("contrast_lower", "contrast_upper", "true_contrast"),
    old = c("delta_lower", "delta_upper", "true_psi")),
  ss_power_equivalence_c = list(
    new = c("contrast_lower", "contrast_upper", "true_contrast"),
    old = c("delta_lower", "delta_upper", "true_psi")),
  simulate_ancova_data = list(
    new = c("true_mean_Y", "true_mean_X", "true_cor"),
    old = c("mu_y", "mu_x", "rho")),
  simulate_ancova_factorial_data = list(
    new = c("true_mean_Y", "true_mean_X", "true_cor_YX", "true_cor_XX"),
    old = c("mu_y", "mu_x", "rho_y_x", "rho_x_x")),
  simulate_anova_data = list(new = "true_mean", old = "mu"),
  ss_aipe_c = list(new = "true_var_error",
                   old = c("error_variance", "true_error_variance",
                           "true_error_var")),
  ss_aipe_c_ancova = list(
    new = c("true_var_error_ancova", "true_var_error_anova", "true_cor"),
    old = c("error_var_ancova", "error_var_anova", "rho",
            "true_error_var_ancova", "true_error_var_anova", "true_rho")),
  ss_aipe_c_ancova_sensitivity = list(
    new = c("true_var_error_ancova", "true_var_error_anova",
            "supposed_var_error_ancova", "supposed_var_error_anova",
            "true_cor", "supposed_cor", "true_mean_Y", "true_mean_X"),
    old = c("true_error_var_ancova", "true_error_var_anova",
            "est_error_var_ancova", "est_error_var_anova", "rho", "est_rho",
            "mu_y", "mu_x", "estimated_error_var_ancova",
            "estimated_error_var_anova", "supposed_error_var_ancova",
            "supposed_error_var_anova", "true_rho", "estimated_rho",
            "true_mu_Y", "true_mu_X")),
  ss_aipe_c_sensitivity = list(
    new = c("true_var_error", "supposed_var_error", "true_contrast"),
    old = c("true_error_variance", "estimated_error_variance", "true_psi",
            "true_error_var", "supposed_error_var")),
  ss_aipe_pcm = list(new = c("true_var_trend", "true_var_error",
                             "true_var_trend_error"),
                     old = c("variance_trend", "error_variance",
                             "true_variance_trend", "true_error_variance",
                             "true_error_var",
                             "variance_true_minus_estimated_trend")),
  ss_aipe_pcm_sensitivity = list(
    new = c("true_var_trend", "true_var_error", "supposed_var_trend",
            "supposed_var_error"),
    old = c("true_variance_trend", "true_error_variance",
            "estimated_variance_trend", "estimated_error_variance",
            "true_error_var", "supposed_error_var")),
  ss_aipe_sc = list(new = "true_std_contrast",
                    old = c("psi_standardized", "true_psi_standardized")),
  ss_aipe_sc_ancova = list(
    new = c("true_contrast", "true_std_contrast", "true_cor"),
    old = c("psi", "psi_standardized", "rho", "true_psi",
            "true_psi_standardized", "true_rho")),
  ss_aipe_sc_ancova_sensitivity = list(
    new = c("true_std_contrast", "supposed_std_contrast", "true_mean_X",
            "true_cor"),
    old = c("true_psi", "estimated_psi", "mu_x", "rho",
            "true_psi_standardized", "estimated_psi_standardized",
            "true_mu_X", "true_rho")),
  ss_aipe_sc_sensitivity = list(
    new = c("true_std_contrast", "supposed_std_contrast"),
    old = c("true_psi", "estimated_psi", "true_psi_standardized",
            "estimated_psi_standardized")),
  ss_power_c = list(new = "true_contrast", old = c("psi", "true_psi")),
  ss_power_c_ancova = list(new = c("true_contrast", "true_cor"),
                           old = c("psi", "rho", "true_psi", "true_rho")),
  ss_power_composite_ancova = list(
    new = "true_cor", old = c("correlations", "true_correlations")),
  ss_power_composite_ancova_2group = list(new = "true_cor",
                                          old = c("rho", "true_rho")),
  ss_power_composite_factorial_ancova_het = list(
    new = "true_cor", old = c("correlations", "true_correlations")),
  ss_power_contrast = list(
    new = c("true_contrast", "true_var_error"),
    old = c("psi", "mu", "sigma_squared", "true_psi",
            "true_error_variance", "true_error_var")),
  ss_power_pcm = list(
    new = c("true_reg_coef", "true_std_reg_coef", "true_var_trend",
            "true_var_error"),
    old = c("beta", "tau", "level_1_variance", "true_beta",
            "true_variance_trend", "true_error_variance", "true_error_var",
            "standardized")),
  ss_power_rm_anova = list(new = "true_cor", old = c("rho", "true_rho")),
  ss_power_sc = list(new = "true_std_contrast",
                     old = c("psi_standardized", "true_psi_standardized")),
  ss_power_split_plot_anova = list(new = "true_cor",
                                   old = c("rho", "true_rho"))
)

test_that("each function has its 1.1.0 names and none of the names it had before", {
  for (fn in names(.renamed)) {
    fm <- names(formals(get(fn, envir = asNamespace("DMAR"))))
    expect_true(all(.renamed[[fn]]$new %in% fm), info = fn)
    expect_false(any(.renamed[[fn]]$old %in% fm), info = fn)
  }
  for (fn in c("ci_c_ancova", "ci_sc")) {
    fm <- names(formals(get(fn, envir = asNamespace("DMAR"))))
    expect_true("contrast" %in% fm, info = fn)
    expect_false("psi" %in% fm, info = fn)
  }
})

test_that("a name a function had before stops with R's unused argument error", {
  for (fn in names(.renamed)) {
    f <- get(fn, envir = asNamespace("DMAR"))
    for (old in .renamed[[fn]]$old) {
      expect_error(do.call(f, stats::setNames(list(1), old)),
                   "unused argument", info = paste(fn, old))
    }
  }
})

test_that("ci_c_ancova() and ci_sc() take ... and refuse the 1.0.x name psi by class", {
  # ci_c_ancova() passes ... to ci_c_ancova_bp() under the Bryant-Paulson
  # procedure and ignored it under "t", so psi = would have been dropped.
  args <- list(s_ancova = sqrt(29), c_weights = c(1, -1, 0), n_per_group = 10,
               cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
  err <- expect_error(do.call(ci_c_ancova, c(list(psi = -4.5), args)),
                      class = "dmar_argument_invalid")
  expect_identical(err$arg, "psi")
  expect_match(conditionMessage(err), "'contrast'")
  expect_error(do.call(ci_c_ancova, c(list(psi = -4.5, adj_means = c(7.5, 12, 14),
                                           procedure = "bryant_paulson"), args)),
               class = "dmar_argument_invalid")
  # ci_sc() passes ... to ci_nc_t().
  err <- expect_error(do.call(ci_sc, list(psi = 3, s_anova = .8,
                                          c_weights = c(.5, -.5, -.5, .5),
                                          n_per_group = 3, N = 12)),
                      class = "dmar_argument_invalid")
  expect_identical(err$arg, "psi")
  expect_match(conditionMessage(err), "'contrast'")
  # The new name gives the same interval as the means it summarizes.
  expect_equal(
    ci_sc(contrast = 3, s_anova = .8, c_weights = c(.5, -.5, -.5, .5),
          n_per_group = 3, N = 12)$value,
    ci_sc(means = c(3, 0, 0, 3), s_anova = .8,
          c_weights = c(.5, -.5, -.5, .5), n_per_group = 3, N = 12)$value,
    tolerance = 1e-12)
})

test_that("the rows that echo a renamed input carry its new name", {
  eq <- equivalence_c(contrast = 1.2, se = 1.1, df_error = 120,
                      contrast_lower = 3, contrast_upper = 8)
  expect_true(all(c("contrast", "contrast_lower", "contrast_upper") %in% eq$term))
  expect_false(any(c("delta_lower", "delta_upper") %in% eq$term))
  expect_equal(eq$value[eq$term == "contrast_lower"], -3)
  expect_equal(eq$value[eq$term == "contrast_upper"], 8)

  pcm <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003,
                      true_var_error = .0262, frequency = 1, duration = 4,
                      specified_n_per_group = 119)
  expect_true(all(c("true_var_error", "true_var_trend") %in% pcm$term))
  expect_false(any(c("true_error_variance", "true_variance_trend",
                     "true_error_var") %in% pcm$term))

  # Correlations of equal magnitude keep the powers exact (no approximation
  # warning) and still give one echo row per group.
  comp <- ss_power_composite_ancova_2group(true_smd = 0.5, true_cor = c(0.4, -0.4),
                                           specified_n_per_group = 40)
  expect_true(all(c("true_cor_group_1", "true_cor_group_2") %in% comp$term))
  expect_false(any(grepl("rho", comp$term)))

  het <- ss_power_composite_factorial_ancova_het(
    factor_levels = 2,
    effects = list(list(type = "mean", factors = 1),
                   list(type = "slope", factors = 1)),
    true_means = c(-0.25, 0.25), true_cor = c(0.3, -0.3), true_sd = 1,
    specified_n_per_cell = 30)
  expect_equal(as.vector(attr(het, "true_cor")), c(0.3, -0.3))
  expect_null(attr(het, "true_correlations", exact = TRUE))
})

test_that("the sensitivity functions echo true_ and supposed_ values and summarize in English names", {
  skip_on_cran()
  c_sens <- ss_aipe_c_sensitivity(true_var_error = 4, supposed_var_error = 4,
                                  c_weights = c(-1, 0, 1), width = 1,
                                  true_contrast = 0.5, G = 5, seed = 113)
  expect_true(all(c("mean_contrast", "median_contrast", "sd_contrast",
                    "true_var_error", "supposed_var_error", "true_contrast")
                  %in% c_sens$term))
  expect_equal(c_sens$value[c_sens$term == "true_contrast"], 0.5)
  # A population variance is true_var_<what> and its planning value
  # supposed_var_<what> (Ken's rule of 2026-11-11), so no row is spelled
  # with error_var.
  expect_false(any(grepl("error_var", c_sens$term)))

  pcm_sens <- ss_aipe_pcm_sensitivity(true_var_trend = 0.003,
                                      true_var_error = 0.0262,
                                      supposed_var_trend = 0.003,
                                      supposed_var_error = 0.0262,
                                      duration = 4, frequency = 1, width = 0.05,
                                      G = 3, seed = 113)
  expect_true(all(c("true_var_trend", "true_var_error", "supposed_var_trend",
                    "supposed_var_error") %in% pcm_sens$term))
  expect_false(any(grepl("^estimated_|variance_trend|error_var",
                         pcm_sens$term)))

  sc_sens <- suppressMessages(ss_aipe_sc_sensitivity(
    true_std_contrast = 0.5, supposed_std_contrast = 0.5,
    c_weights = c(-1, 0, 1), width = 0.8, G = 5, print_iter = FALSE,
    seed = 113))
  expect_true(all(c("mean_std_contrast", "median_std_contrast",
                    "sd_std_contrast", "true_std_contrast",
                    "supposed_std_contrast") %in% sc_sens$term))
  expect_false(any(grepl("psi|^estimated_", sc_sens$term)))

  c_anc <- suppressMessages(ss_aipe_c_ancova_sensitivity(
    true_var_error_ancova = 0.8, supposed_var_error_anova = 1, supposed_cor = 0.4,
    true_cor = 0.4, G = 5, true_mean_Y = c(10, 11, 12), true_sd_Y = 1,
    true_mean_X = 5, true_sd_X = 1, c_weights = c(-1, 0, 1), width = 0.8,
    seed = 113))
  expect_true(all(c("mean_contrast", "median_contrast", "sd_contrast",
                    "true_contrast", "supposed_var_error_ancova", "true_cor")
                  %in% c_anc$term))
  expect_equal(c_anc$value[c_anc$term == "supposed_var_error_ancova"],
               1 * (1 - 0.4^2))
  expect_false(any(grepl("psi|rho|^estimated_|error_var", c_anc$term)))
})

test_that("the limits of an interval named for its estimate end in _ci_lower and _ci_upper", {
  vc <- variance_components_mls(ms_between = 40, ms_within = 10,
                                df_between = 4, df_within = 45,
                                n_per_group = 10)
  expect_identical(vc$term,
                   c("sigma2_between", "sigma2_between_ci_lower",
                     "sigma2_between_ci_upper", "sigma2_within",
                     "sigma2_within_ci_lower", "sigma2_within_ci_upper",
                     "icc", "icc_ci_lower", "icc_ci_upper"))
  pyg <- pygmalion
  pyg$grade_band <- factor(ifelse(pyg$grade <= 2, "lower", "upper"))
  fa <- factorial_anova(iq_gain ~ treatment * grade_band, data = pyg)
  expect_true(all(c("eta_squared_partial_ci_lower", "eta_squared_partial_ci_upper",
                    "omega_squared_partial_ci_lower",
                    "omega_squared_partial_ci_upper") %in% names(fa)))
  expect_false(any(grepl("partial_(lower|upper)$", names(fa))))
})

test_that("tidy() and results_sentence() read the limits named for their estimate", {
  # Renamed to <estimate>_ci_lower and <estimate>_ci_upper, the limits are
  # read as every such pair is: tidy() marks the rows of each interval with
  # its level in a conf_level column (NA on the other rows), and
  # results_sentence() reports the first estimate that has such a pair, or
  # the one term = names. Through b53fee6 tidy() gave term and estimate
  # alone and results_sentence() stopped, as the names were not read.
  a <- ancova(outcome = "iq_8", treatment = "treatment",
              covariates = "iq_pre", data = pygmalion)
  td <- generics::tidy(a)
  expect_identical(names(td), c("term", "estimate", "conf_level"))
  ci_rows <- td$term %in% c("eta_squared_partial", "eta_squared_partial_ci_lower",
                            "eta_squared_partial_ci_upper",
                            "omega_squared_partial",
                            "omega_squared_partial_ci_lower",
                            "omega_squared_partial_ci_upper")
  expect_true(all(td$conf_level[ci_rows] == 0.95))
  expect_true(all(is.na(td$conf_level[!ci_rows])))
  expect_identical(results_sentence(a),
                   "eta_squared_partial = 0.02, 95% CI [0.00, 0.06]")
  expect_identical(results_sentence(td),
                   "eta_squared_partial = 0.02, 95% CI [0.00, 0.06]")
  expect_identical(results_sentence(a, term = "omega_squared_partial"),
                   "omega_squared_partial = 0.01, 95% CI [0.00, 0.06]")

  vc <- variance_components_mls(ms_between = 40, ms_within = 10,
                                df_between = 4, df_within = 45,
                                n_per_group = 10)
  expect_identical(names(generics::tidy(vc)), c("term", "estimate", "conf_level"))
  expect_true(all(generics::tidy(vc)$conf_level == 0.95))
  expect_identical(results_sentence(vc), "sigma2_between = 3.00, 95% CI [0.35, 31.99]")
  expect_identical(results_sentence(vc, term = "icc"), "icc = 0.23, 95% CI [0.02, 0.82]")

  # A row of factorial_anova() is read by the estimate term = names; the
  # table as a whole asks for a row and for the estimate.
  pyg <- pygmalion
  pyg$grade_band <- factor(ifelse(pyg$grade <= 2, "lower", "upper"))
  fa <- suppressWarnings(factorial_anova(iq_gain ~ treatment * grade_band,
                                         data = pyg))
  expect_error(results_sentence(fa), "several estimates")
  expect_identical(
    results_sentence(fa[fa$effect == "treatment", ], term = "eta_squared_partial"),
    "eta_squared_partial = 0.03, 95% CI [0.00, 0.07], unadjusted p = 0.0044")

  # The other functions whose limits were renamed in the same pass.
  te <- teacher_expectancy
  v <- (te$n_experimental + te$n_control) / (te$n_experimental * te$n_control) +
    te$d^2 / (2 * (te$n_experimental + te$n_control))
  m <- meta_es(te$d, v)
  expect_identical(results_sentence(m),
                   "estimate = 0.05, 95% CI [-0.02, 0.13], unadjusted p = 0.1401")
  expect_identical(results_sentence(m, term = "tau2"), "tau2 = 0.00, 95% CI [0.00, 0.05]")
  tm <- generics::tidy(m)
  expect_true(all(tm$conf_level[tm$term %in% c("tau2", "tau2_ci_lower", "tau2_ci_upper",
                                               "I2", "I2_ci_lower", "I2_ci_upper")] == 0.95))
  expect_true(all(is.na(tm$conf_level[tm$term %in% c("prediction_lower",
                                                     "prediction_upper", "tau")])))
  nn <- nnt_from_smd(0.5, n_1 = 40, n_2 = 40)
  expect_identical(results_sentence(nn), "srd = 0.28, 95% CI [0.03, 0.50]")
  expect_identical(results_sentence(nn, term = "nnt"), "nnt = 3.62, 95% CI [2.02, 33.30]")
  expect_identical(results_sentence(cles(smd = 0.5, n_1 = 40, n_2 = 40)),
                   "cles = 0.64, 95% CI [0.52, 0.75]")
  set.seed(11)
  xx <- rnorm(80); mm <- 0.5 * xx + rnorm(80); yy <- 0.4 * mm + 0.2 * xx + rnorm(80)
  med <- mediate(data.frame(x = xx, m = mm, y = yy), "x", "m", "y",
                 ci_method = "sobel")
  expect_identical(results_sentence(med), "indirect_effect = 0.10, 95% CI [0.00, 0.20]")
  expect_match(results_sentence(med, term = "total_effect"), "^total_effect = ")
})

test_that("the refusal of psi reads the name from the call and evaluates nothing in ...", {
  # The old name is data and its value an unevaluated call that stops if
  # it is forced: do.call() passes it on as an expression.
  old <- stats::setNames(list(quote(stop("an argument in ... was evaluated"))),
                         "psi")
  err <- expect_error(
    do.call(ci_c_ancova, c(old, list(s_ancova = sqrt(29), c_weights = c(1, -1, 0),
                                     n_per_group = 10, cov_means = c(17, 17.7, 17.4),
                                     SSwithin_x = 752.5))),
    class = "dmar_argument_invalid")
  expect_identical(err$arg, "psi")
  err <- expect_error(
    do.call(ci_sc, c(old, list(s_anova = .8, c_weights = c(.5, -.5, -.5, .5),
                               n_per_group = 3, N = 12))),
    class = "dmar_argument_invalid")
  expect_identical(err$arg, "psi")
})
