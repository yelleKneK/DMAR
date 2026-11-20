# The argument and row renames of DMAR 1.1.0 for the standardized mean
# difference, equivalence, means, conversion, and multilevel sensitivity
# functions: English quantity names, true_ for a population value,
# supposed_ for the planning value of a sensitivity analysis, null_ for the
# value of a null hypothesis, and explicit lower and upper bounds. There is
# no shim: a DMAR 1.0.x name, or the name an argument had earlier in the
# 1.1.0 cycle, stops with R's "unused argument" error, or with a refusal
# that names it where R would otherwise match it to another argument or pass
# it on through `...` (n in bayes_one_sample_t(), mu_0 in
# ci_mahalanobis()). The calls are built with do.call() so that each old
# name is data here and not a call a rename sweep would change.

unused <- function(f, args, info = NULL) {
  expect_error(do.call(f, args), "unused argument", info = info)
}

test_that("each renamed function has its new names and none of the old ones", {
  renamed <- list(
    # The 1.0.x name n is the first letter of null_mean, so R would match
    # it; the function reads it from the call and refuses it (tested below).
    bayes_one_sample_t = list(new = c("null_mean", "N"), old = "mu_0",
                              refused = "n"),
    ci_mahalanobis = list(new = "null_mean", old = "mu_0"),
    convert_delta_lambda = list(new = "smd", old = "delta"),
    convert_lambda_delta = list(new = "ncp", old = "lambda"),
    convert_lambda_R2 = list(new = "ncp", old = "lambda"),
    equivalence_proportion_paired = list(
      new = c("prop_diff_lower", "prop_diff_upper"),
      old = c("delta_lower", "delta_upper")),
    equivalence_smd = list(new = c("smd_lower", "smd_upper"),
                           old = c("delta_lower", "delta_upper")),
    expected_smd = list(new = "true_smd", old = "delta"),
    plot_equivalence = list(new = c("effect_lower", "effect_upper"),
                            old = c("delta_lower", "delta_upper")),
    power_density_equivalence_md = list(
      new = c("lower_bound", "upper_bound", "true_effect", "df_error"),
      old = c("theta1", "theta2", "diff", "nu", "power_sigma", "sigma")),
    power_equivalence_md = list(
      new = c("lower_bound", "upper_bound", "true_effect", "df_error"),
      old = c("ltheta1", "ltheta2", "ldiff", "nu", "sigma")),
    power_equivalence_md_plot = list(
      new = c("lower_bound", "upper_bound", "df_error"),
      old = c("theta1", "theta2", "nu", "sigma")),
    ss_aipe_cliff_delta_sensitivity = list(new = "supposed_delta",
                                           old = "estimated_delta"),
    ss_aipe_cv_sensitivity = list(new = "supposed_cv",
                                  old = c("estimated_cv", "mean")),
    ss_aipe_equivalence_smd_sensitivity = list(
      new = c("supposed_smd", "smd_lower", "smd_upper"),
      old = c("estimated_smd", "delta_lower", "delta_upper", "n_per_group")),
    ss_aipe_mean_sensitivity = list(new = "supposed_sd", old = "estimated_sd"),
    ss_aipe_omega_squared_sensitivity = list(
      new = "supposed_omega_squared", old = "estimated_omega_squared"),
    ss_aipe_sm_sensitivity = list(new = "supposed_sm",
                                  old = c("estimated_sm", "desired_width")),
    ss_aipe_smd_sensitivity = list(
      new = "supposed_smd",
      old = c("estimated_smd", "estimated_delta", "true_delta",
              "desired_width", "n_per_group")),
    summary_t_test = list(new = "null_mean_diff", old = "mu"),
    var_smd = list(new = "true_smd", old = "delta"),
    var_cv = list(new = "true_cv", old = "cv"),
    welch_t = list(new = "null_mean_diff", old = "mu"),
    var_omega_squared = list(new = "true_omega_squared",
                             old = "population_omega_squared"),
    ss_aipe_crossed_sensitivity = list(
      new = c("supposed_var_unit", "supposed_var_rater",
              "supposed_var_residual"),
      old = c("estimated_var_unit", "estimated_var_rater",
              "estimated_var_residual")),
    ss_aipe_icc_sensitivity = list(
      new = "supposed_icc",
      old = c("estimated_icc", "estimated_rho", "true_rho", "k")),
    ss_aipe_mixed_effects_sensitivity = list(
      new = c("supposed_total_var", "supposed_var_X", "supposed_icc",
              "supposed_icc_X", "true_reg_coef"),
      old = c("estimated_total_var", "estimated_var_X", "estimated_icc",
              "estimated_icc_X", "true_beta", "true_sigma2_y",
              "estimated_sigma2_y", "true_sigma2_x", "estimated_sigma2_x",
              "specified_K"))
  )
  for (fn in names(renamed)) {
    fm <- names(formals(get(fn)))
    expect_true(all(renamed[[fn]]$new %in% fm), info = fn)
    expect_false(any(c(renamed[[fn]]$old, renamed[[fn]]$refused) %in% fm),
                 info = fn)
    # No old name is the beginning of a new one, so R does not match it
    # partially to another argument.
    for (old in renamed[[fn]]$old) {
      expect_false(any(startsWith(fm, old)), info = paste(fn, old))
    }
    # An old name that is the beginning of a new one would be matched by R;
    # the function refuses it by name instead.
    for (old in renamed[[fn]]$refused) {
      expect_true(any(startsWith(fm, old)), info = paste(fn, old))
    }
  }
})

test_that("bayes_one_sample_t() refuses n, which R would match to null_mean", {
  set.seed(1)
  dat <- stats::rnorm(25, 0.3)
  # (suppressWarnings: under warnPartialMatchArgs R warns of the partial
  # match before the function refuses it; the refusal is the property under
  # test.)
  # The 1.0.x summary call.
  err <- expect_error(suppressWarnings(
    do.call(bayes_one_sample_t, list(mean = .5, sd = 1, n = 30))),
    "N = ", class = "dmar_argument_invalid")
  expect_identical(err$arg, "n")
  expect_match(conditionMessage(err), "'n' is not an argument", fixed = TRUE)
  # Raw data with n: without the check R would take n as null_mean and
  # return an analysis against a comparison value of 25.
  err <- expect_error(suppressWarnings(
    do.call(bayes_one_sample_t, list(x = dat, n = 25))),
    class = "dmar_argument_invalid")
  expect_identical(err$arg, "n")
  # Written directly in the call rather than built.
  expect_error(suppressWarnings(bayes_one_sample_t(dat, n = 25)),
               class = "dmar_argument_invalid")
  # The new names run.
  a <- bayes_one_sample_t(mean = .5, sd = 1, N = 30)
  expect_identical(a$value[a$term == "N"], 30)
  b <- bayes_one_sample_t(x = dat, null_mean = 0)
  expect_identical(b$value, bayes_one_sample_t(x = dat)$value)
  # bayes_paired_t(), which calls it, is not affected.
  expect_identical(bayes_paired_t(x = dat, y = dat / 2)$value,
                   bayes_one_sample_t(dat - dat / 2)$value)
})

test_that("a DMAR 1.0.x name or an earlier 1.1.0 name stops with 'unused argument'", {
  unused(bayes_one_sample_t, list(x = c(1, 2, 3), mu_0 = 0))
  unused(convert_delta_lambda, list(delta = .5, n_1 = 20, n_2 = 20))
  unused(convert_lambda_delta, list(lambda = 2, n_1 = 20, n_2 = 20))
  unused(convert_lambda_R2, list(lambda = 5, N = 100))
  unused(equivalence_proportion_paired,
         list(n_10 = 32, n_01 = 27, N = 344, delta_lower = .05))
  unused(equivalence_proportion_paired,
         list(n_10 = 32, n_01 = 27, N = 344, delta_upper = .05))
  unused(equivalence_smd, list(smd = .1, n_1 = 50, n_2 = 50, delta_upper = .5))
  unused(equivalence_smd, list(smd = .1, n_1 = 50, n_2 = 50, delta_lower = .5,
                               smd_upper = .5))
  unused(expected_smd, list(delta = .5, n_1 = 10))
  unused(var_smd, list(delta = .5, n_1 = 20))
  unused(var_cv, list(cv = .2, N = 30))
  unused(plot_equivalence, list(estimate = 1, lower = 0, upper = 2,
                                delta_upper = 5))
  unused(plot_equivalence, list(estimate = 1, lower = 0, upper = 2,
                                delta_lower = 5, effect_upper = 5))
  for (old in c("ltheta1", "ltheta2", "ldiff", "nu", "sigma")) {
    unused(power_equivalence_md,
           c(list(alpha_level = .05, logscale = FALSE),
             stats::setNames(list(.2), old)), info = old)
  }
  for (old in c("theta1", "theta2", "diff", "nu", "power_sigma", "sigma")) {
    unused(power_density_equivalence_md,
           c(list(sd_values = .1, alpha_level = .05),
             stats::setNames(list(.2), old)), info = old)
  }
  for (old in c("theta1", "theta2", "nu", "sigma")) {
    unused(power_equivalence_md_plot,
           c(list(alpha_level = .05, logscale = FALSE),
             stats::setNames(list(.2), old)), info = old)
  }
  unused(summary_t_test, list(mean_1 = 10, sd_1 = 2, n_1 = 20, mean_2 = 9,
                              sd_2 = 2, n_2 = 20, mu = 0))
  unused(welch_t, list(x = c(1, 2, 3, 5), y = c(2, 4, 4, 6), mu = 0))
  unused(var_omega_squared, list(population_omega_squared = .1,
                                 df_effect = 2, df_error = 57, N = 60))
  unused(ss_aipe_cliff_delta_sensitivity,
         list(true_delta = .3, estimated_delta = .3, width = .3, G = 2))
  for (old in c("estimated_cv", "mean")) {
    unused(ss_aipe_cv_sensitivity,
           c(list(true_cv = .25, width = .1, G = 2),
             stats::setNames(list(.25), old)), info = old)
  }
  for (old in c("estimated_smd", "delta_lower", "delta_upper", "n_per_group")) {
    unused(ss_aipe_equivalence_smd_sensitivity,
           c(list(true_smd = 0, width = .3, G = 2),
             stats::setNames(list(.2), old)), info = old)
  }
  unused(ss_aipe_mean_sensitivity,
         list(true_sd = 10, estimated_sd = 10, width = 4, G = 2))
  unused(ss_aipe_omega_squared_sensitivity,
         list(true_omega_squared = .1, estimated_omega_squared = .1,
              df_effect = 2, width = .1, G = 2))
  for (old in c("estimated_sm", "desired_width")) {
    unused(ss_aipe_sm_sensitivity,
           c(list(true_sm = 1.5, G = 2), stats::setNames(list(.5), old)),
           info = old)
  }
  for (old in c("estimated_smd", "estimated_delta", "true_delta",
                "desired_width", "n_per_group")) {
    unused(ss_aipe_smd_sensitivity,
           c(list(true_smd = .5, G = 2), stats::setNames(list(.5), old)),
           info = old)
  }
  unused(ss_aipe_crossed_sensitivity,
         list(true_var_unit = .5, true_var_rater = .1, true_var_residual = .4,
              estimated_var_unit = .5, estimated_var_rater = .1,
              estimated_var_residual = .4, width = .4, G = 2))
  for (old in c("estimated_icc", "estimated_rho", "true_rho", "k")) {
    unused(ss_aipe_icc_sensitivity,
           c(list(true_icc = .7, n_raters = 3, width = .4, G = 2),
             stats::setNames(list(.7), old)), info = old)
  }
  for (old in c("estimated_total_var", "estimated_var_X", "estimated_icc",
                "estimated_icc_X", "true_beta", "true_sigma2_y",
                "estimated_sigma2_y", "true_sigma2_x", "estimated_sigma2_x",
                "specified_K")) {
    unused(ss_aipe_mixed_effects_sensitivity,
           c(list(true_total_var = 1, true_var_X = 1, true_icc = .1,
                  width = .4, G = 2), stats::setNames(list(.1), old)),
           info = old)
  }
})

test_that("ci_mahalanobis() refuses mu_0, which its ... would otherwise pass on", {
  set.seed(113)
  g <- matrix(stats::rnorm(60, mean = 1), 20, 3)
  err <- expect_error(do.call(ci_mahalanobis,
                              list(group_1 = g, mu_0 = c(0, 0, 0))),
                      "null_mean", class = "dmar_argument_invalid")
  expect_identical(err$arg, "mu_0")
  # Under its new name the value is used.
  a <- ci_mahalanobis(group_1 = g, null_mean = c(0, 0, 0))
  b <- ci_mahalanobis(group_1 = g, null_mean = c(0.2, 0.2, 0.2))
  expect_false(isTRUE(all.equal(a$D2, b$D2)))
})

test_that("the rows that echo a renamed input take its new name", {
  eq <- equivalence_smd(smd = .1, n_1 = 50, n_2 = 50, smd_lower = .4,
                        smd_upper = .5)
  expect_identical(eq$value[eq$term == "smd_lower"], -.4)
  expect_identical(eq$value[eq$term == "smd_upper"], .5)
  expect_false(any(c("delta_lower", "delta_upper") %in% eq$term))

  pp <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                      prop_diff_lower = .05,
                                      prop_diff_upper = .08)
  expect_identical(pp$value[pp$term == "prop_diff_lower"], -.05)
  expect_identical(pp$value[pp$term == "prop_diff_upper"], .08)
  expect_false(any(c("delta_lower", "delta_upper") %in% pp$term))

  es <- expected_smd(true_smd = c(.2, .5), n_1 = 20)
  expect_identical(names(es)[1], "true_smd")
  expect_identical(es$true_smd, c(.2, .5))
  expect_false("delta" %in% names(es))
  # A whole true_smd prints to the decimals of the rest, as delta did.
  expect_identical(format(expected_smd(true_smd = 1, n_1 = 20))$true_smd,
                   "1.000")

  echoes <- list(
    ss_aipe_smd_sensitivity = list(
      call = quote(ss_aipe_smd_sensitivity(true_smd = .5, supposed_smd = .5,
                                           width = .8, G = 2, seed = 1)),
      new = "supposed_smd", old = c("estimated_smd", "estimated_delta")),
    ss_aipe_cliff_delta_sensitivity = list(
      call = quote(ss_aipe_cliff_delta_sensitivity(
        true_delta = .3, supposed_delta = .3, width = .5, G = 2, seed = 1)),
      new = "supposed_delta", old = "estimated_delta"),
    ss_aipe_cv_sensitivity = list(
      call = quote(ss_aipe_cv_sensitivity(true_cv = .25, supposed_cv = .25,
                                          width = .2, G = 2, seed = 1)),
      new = "supposed_cv", old = "estimated_cv"),
    ss_aipe_equivalence_smd_sensitivity = list(
      call = quote(ss_aipe_equivalence_smd_sensitivity(
        true_smd = 0, supposed_smd = 0, width = .6, smd_lower = .4,
        smd_upper = .5, G = 2, seed = 1)),
      new = c("supposed_smd", "smd_lower", "smd_upper"),
      old = c("estimated_smd", "delta_lower", "delta_upper")),
    ss_aipe_mean_sensitivity = list(
      call = quote(ss_aipe_mean_sensitivity(true_sd = 2, supposed_sd = 2,
                                            width = 1, G = 2, seed = 1)),
      new = "supposed_sd", old = "estimated_sd"),
    ss_aipe_omega_squared_sensitivity = list(
      call = quote(ss_aipe_omega_squared_sensitivity(
        true_omega_squared = .2, supposed_omega_squared = .2, df_effect = 2,
        width = .3, G = 2, seed = 1)),
      new = "supposed_omega_squared", old = "estimated_omega_squared"),
    ss_aipe_sm_sensitivity = list(
      call = quote(ss_aipe_sm_sensitivity(true_sm = 5, supposed_sm = 5,
                                          width = 1, G = 2, seed = 1,
                                          print_iter = FALSE)),
      new = "supposed_sm", old = "estimated_sm"),
    ss_aipe_crossed_sensitivity = list(
      call = quote(ss_aipe_crossed_sensitivity(
        true_var_unit = .5, true_var_rater = .1, true_var_residual = .4,
        supposed_var_unit = .5, supposed_var_rater = .1,
        supposed_var_residual = .4, width = .6, G = 2, seed = 1)),
      new = c("supposed_var_unit", "supposed_var_rater",
              "supposed_var_residual"),
      old = c("estimated_var_unit", "estimated_var_rater",
              "estimated_var_residual")),
    ss_aipe_icc_sensitivity = list(
      call = quote(ss_aipe_icc_sensitivity(true_icc = .8, supposed_icc = .8,
                                           n_raters = 4, width = .4, G = 2,
                                           seed = 1)),
      new = "supposed_icc", old = c("estimated_icc", "estimated_rho")),
    ss_aipe_mixed_effects_sensitivity = list(
      call = quote(ss_aipe_mixed_effects_sensitivity(
        true_total_var = 1, true_var_X = 1, true_icc = .2,
        true_reg_coef = .3, supposed_total_var = 1, supposed_var_X = 1,
        supposed_icc = .2, cluster_size = 10, width = .4, G = 2, seed = 1)),
      new = c("supposed_total_var", "supposed_var_X", "supposed_icc",
              "supposed_icc_X", "true_reg_coef", "mean_reg_coef",
              "median_reg_coef", "sd_reg_coef"),
      old = c("estimated_total_var", "estimated_var_X", "estimated_icc",
              "estimated_icc_X", "true_beta", "mean_beta", "median_beta",
              "sd_beta"))
  )
  for (member in names(echoes)) {
    e <- echoes[[member]]
    # The planners' messages and warnings are not the property under
    # test.
    res <- suppressMessages(suppressWarnings(eval(e$call)))
    expect_true(all(e$new %in% res$term), info = member)
    expect_false(any(e$old %in% res$term), info = member)
    expect_false(any(startsWith(res$term, "estimated_")), info = member)
  }
})
