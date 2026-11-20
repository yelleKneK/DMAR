# One notation, no silent arguments, and classed refusals: the checks that
# closed the 1.1.0 notation work. Each block states the rule it holds.

.term <- function(x, term) x$value[x$term == term]
.caught <- function(expr) tryCatch(expr, error = function(e) e)

# ---- No silent `...` --------------------------------------------------------

test_that("the planners that used `...` for nothing no longer take it", {
  for (f in list(ci_cv, ss_aipe_c, ss_aipe_cv, ss_power_R2, ss_aipe_R2,
                 ss_aipe_sc, ss_aipe_sc_ancova, ss_aipe_sm,
                 ss_aipe_sc_ancova_sensitivity)) {
    expect_false("..." %in% names(formals(f)))
  }
})

test_that("a name a planner does not have is an unused argument, not a different answer", {
  # ss_power_r() took N in 1.0.x, so N = is the likely slip here; it planned
  # for the default power and returned 36.
  expect_error(ss_power_R2(true_R2 = .3, p = 3, N = 50), "unused argument")
  expect_error(ss_power_R2(true_R2 = .2, population_R2 = .6, p = 4),
               "unused argument")
  expect_error(suppressWarnings(
    ss_aipe_R2(true_R2 = .3, p = 3, width = .2, conf.level = .99)),
    "unused argument")
  expect_error(ss_aipe_sm(true_sm = .5, width = .3, N = 20), "unused argument")
  expect_error(ss_aipe_c(true_var_error = 1, c_weights = c(1, -1),
                         width = .5, n = 10), "unused argument")
  expect_error(ss_aipe_sc(true_std_contrast = .5, c_weights = c(-1, 0, 1),
                          width = .8, desired_width = .2), "unused argument")
  expect_error(ss_aipe_sc_ancova(true_std_contrast = .8, psi = 4,
                                 true_sd_anova = 5, true_cor = .4, width = .5,
                                 c_weights = c(.5, .5, 0, -1)),
               "unused argument")
  expect_error(ci_cv(cv = .2, N = 30, n = 300), "unused argument")
})

test_that("ss_aipe_R2() takes the two simulation counts of verify_ss as arguments", {
  expect_true(all(c("G_initial", "G") %in% names(formals(ss_aipe_R2))))
  expect_identical(formals(ss_aipe_R2)$G_initial, 500)
  expect_identical(formals(ss_aipe_R2)$G, 10000)
})

# ---- One name for one quantity ---------------------------------------------

test_that("the null correlation of ss_power_r() is null_cor, as an argument and a row", {
  expect_identical(names(formals(ss_power_r))[1:2], c("true_cor", "null_cor"))
  res <- ss_power_r(true_cor = .40, null_cor = .20, desired_power = .80,
                    directional = TRUE)
  expect_equal(.term(res, "null_cor"), .20)
  expect_equal(.term(res, "true_cor"), .40)
  expect_false(any(c("rho_0", "null_rho", "true_rho") %in% res$term))
  # The 1.0.x planning value is no longer taken for the null value.
  expect_error(ss_power_r(.3, rho = .1, desired_power = .8), "unused argument")
  expect_error(ss_power_r(rho = .3, desired_power = .8), "unused argument")
  expect_error(ss_power_r(true_cor = .3, rho_0 = .1), "unused argument")
  # Nor are the names the first builds of 1.1.0 used.
  expect_error(ss_power_r(true_rho = .3, desired_power = .8),
               "unused argument")
  expect_error(ss_power_r(true_cor = .3, null_rho = .1, desired_power = .8),
               "unused argument")
})

test_that("the response and covariate correlation is true_cor in every ANCOVA planner", {
  for (f in list(ss_aipe_c_ancova, ss_aipe_c_ancova_sensitivity,
                 ss_aipe_sc_ancova, ss_aipe_sc_ancova_sensitivity,
                 ss_power_c_ancova, ss_power_composite_ancova_2group)) {
    expect_true("true_cor" %in% names(formals(f)))
    expect_false(any(c("rho", "true_rho") %in% names(formals(f))))
  }
  expect_error(ss_power_c_ancova(true_contrast = .5, c_weights = c(.5, -.5),
                                 true_sd = 1, rho = .5, desired_power = .8),
               "unused argument")
})

test_that("the desired width is the row `width` and the expected width `expected_width`", {
  plans <- list(
    ss_aipe_r(true_cor = .3, width = .2),
    ss_aipe_partial_r(true_cor = .3, p = 3, width = .2),
    suppressWarnings(ss_aipe_semipartial_r(true_r_sp = .3, p = 3, width = .2)),
    ss_aipe_icc(true_icc = .6, n_raters = 3, width = .2),
    ss_aipe_cliff_delta(true_delta = .3, width = .3),
    ss_aipe_mean(true_sd = 1, width = .3)
  )
  for (res in plans) {
    expect_true(all(c("width", "expected_width") %in% res$term))
    expect_false(any(c("width_target", "ci_width_expected") %in% res$term))
  }
  expect_equal(.term(plans[[1]], "width"), .2)
  mixed <- ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1,
                                 true_icc = .1, width = .2)
  expect_true(all(c("width", "expected_width") %in% mixed$term))
  expect_false("ci_width_expected" %in% mixed$term)
  dc <- design_consequences(true_effect = .2, se = .1, width = .3)
  expect_equal(.term(dc, "width"), .3)
  expect_false("width_target" %in% dc$term)
})

test_that("the coefficient and the mean a plan is made for are true_ values", {
  expect_true("true_reg_coef_j" %in% names(formals(ss_aipe_reg_coef)))
  expect_true("true_reg_coef_j" %in% names(formals(ss_aipe_rc)))
  expect_true("true_std_reg_coef_j" %in% names(formals(ss_aipe_src)))
  expect_true("true_reg_coef_j" %in% names(formals(ss_power_reg_coef)))
  expect_true("true_reg_coef_j" %in% names(formals(ss_power_rc)))
  for (f in list(ss_aipe_reg_coef, ss_aipe_rc, ss_aipe_src, ss_power_reg_coef,
                 ss_power_rc)) {
    expect_false(any(c("b_j", "beta_j", "true_b_j", "true_beta_j") %in%
                       names(formals(f))))
  }
  expect_true(all(c("true_mean", "true_sd") %in% names(formals(ss_aipe_cv))))
  expect_false(any(c("mu", "sigma") %in% names(formals(ss_aipe_cv))))
  expect_true("true_mean" %in% names(formals(ss_aipe_cv_sensitivity)))
  expect_false("mean" %in% names(formals(ss_aipe_cv_sensitivity)))
  # The two ways of stating the coefficient of variation agree.
  expect_identical(ss_aipe_cv(true_mean = 10, true_sd = 1, width = .1)$value,
                   ss_aipe_cv(true_cv = .1, width = .1)$value)
})

test_that("rows follow their arguments in the power planners and design_effect()", {
  pcm <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003,
                      true_var_error = .0262, frequency = 1, duration = 4,
                      desired_power = .8)
  expect_equal(.term(pcm, "frequency"), 1)
  expect_false("freq" %in% pcm$term)

  de <- design_effect(cluster_sizes = c(10, 12, 9), icc = .1)
  expect_equal(.term(de, "N"), 31)
  expect_false("n_total" %in% de$term)

  # The one-way planner keeps n_per_group, a name tidy() does not read as a
  # size, so that tidy() reports the total in the planning branch as it
  # does in the specified branch, which has no per-group row.
  one_way <- ss_power_one_way_anova(a = 3, true_f = .25, desired_power = .8)
  expect_equal(.term(one_way, "n_per_group") * 3, .term(one_way, "necessary_N"))
  expect_equal(generics::tidy(one_way)$estimate, .term(one_way, "necessary_N"))

  # Cohen's f2 at the planning value has one row name in the two planners
  # for a squared multiple correlation.
  omnibus <- ss_power_R2(true_R2 = .2, p = 4, random_predictors = FALSE)
  change <- ss_power_R2_change(true_R2_full = .2, true_R2_restricted = 0,
                               p_full = 4, p_restricted = 0,
                               random_predictors = FALSE)
  expect_equal(.term(omnibus, "true_f2"), .25)
  expect_equal(.term(omnibus, "true_f2"), .term(change, "true_f2"))
  expect_false("effect_size" %in% omnibus$term)
})

# ---- Refusals ---------------------------------------------------------------

test_that("desired_power outside (0, 1) is refused when a sample size is given too", {
  calls <- list(
    function(dp) ss_power_R2(true_R2 = .3, p = 3, specified_N = 60, desired_power = dp),
    function(dp) ss_power_R2_change(true_R2_full = .3, true_R2_restricted = .2,
                                    p_full = 5, p_restricted = 3,
                                    specified_N = 60, desired_power = dp),
    function(dp) ss_power_r(true_cor = .3, specified_N = 60, desired_power = dp),
    function(dp) ss_power_smd(true_smd = .5, n_1 = 30, n_2 = 30, desired_power = dp),
    function(dp) ss_power_c(true_contrast = 1, c_weights = c(1, -1), true_sd = 2,
                            specified_n_per_group = 20, desired_power = dp),
    function(dp) ss_power_sc(true_std_contrast = .5, c_weights = c(1, -1),
                             specified_n_per_group = 20, desired_power = dp),
    function(dp) ss_power_c_ancova(true_contrast = 1, c_weights = c(1, -1),
                                   true_sd = 2, true_cor = .4, specified_n_per_group = 20,
                                   desired_power = dp),
    function(dp) ss_power_one_way_anova(a = 3, true_f = .25, specified_N = 60,
                                        desired_power = dp),
    function(dp) ss_power_rm_anova(a = 3, true_f = .25, true_cor = .5, specified_N = 30,
                                   desired_power = dp),
    function(dp) ss_power_factorial_anova(factor_levels = c(2, 3),
                                          effect_indices = 2, true_f = .25,
                                          specified_n_per_cell = 10, desired_power = dp),
    function(dp) ss_power_mixed_effects(true_smd = .4, cluster_size = 20, true_icc = .1, specified_n_clusters_per_arm = 10,
                                        desired_power = dp),
    function(dp) ss_power_composite_ancova_2group(true_smd = .5, true_cor = .5,
                                                  specified_n_per_group = 100, desired_power = dp)
  )
  for (call in calls) {
    # The call runs with an admissible value ...
    expect_s3_class(call(.80), "data.frame")
    # ... and an inadmissible one is refused by class, though it is not used.
    for (bad in list(1, 0, 7, c(.8, .9))) {
      expect_s3_class(.caught(call(bad)), "dmar_desired_power_invalid")
    }
  }
})

test_that("ss_aipe_sc_ancova_sensitivity() asks for width when it plans a sample size", {
  err <- .caught(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, supposed_std_contrast = .5,
    c_weights = c(1, -1), true_cor = .3, G = 10, seed = 1,
    print_iter = FALSE))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_s3_class(err, "dmar_refusal")
  expect_identical(err$arg, "width")
  expect_match(conditionMessage(err), "'width' must be a single positive number")
  # The 1.0.x name is an unused argument, not "argument is of length zero".
  expect_error(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, supposed_std_contrast = .5,
    c_weights = c(1, -1), desired_width = .5, true_cor = .3, G = 10,
    seed = 1), "unused argument")
})

test_that("ss_aipe_R2_sensitivity() refuses supposed_R2 together with specified_N", {
  err <- .caught(ss_aipe_R2_sensitivity(
    true_R2 = .02, supposed_R2 = .02, width = .3, p = 3, specified_N = 6,
    G = 50, print_iter = FALSE, seed = 2))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, c("supposed_R2", "specified_N"))
  expect_match(conditionMessage(err),
               "You must specify 'supposed_R2' or 'specified_N', but not both.",
               fixed = TRUE)
  # The siblings give the same refusal by the same class.
  expect_s3_class(.caught(ss_aipe_cv_sensitivity(
    true_cv = .5, supposed_cv = .5, width = .5, specified_N = 4, G = 30)),
    "dmar_argument_invalid")
  expect_s3_class(.caught(ss_aipe_cliff_delta_sensitivity(
    true_delta = .3, supposed_delta = .3, width = .3, specified_N = 20,
    G = 30)), "dmar_argument_invalid")
})

test_that("the consistency refusals of ss_power_R2_change() have a class", {
  nested <- .caught(ss_power_R2_change(true_R2_full = .2,
                                       true_R2_restricted = .3,
                                       p_full = 5, p_restricted = 3))
  expect_s3_class(nested, "dmar_argument_invalid")
  expect_identical(nested$arg, c("true_R2_restricted", "true_R2_full"))
  expect_match(conditionMessage(nested), "exceeds 'true_R2_full'")

  counts <- .caught(ss_power_R2_change(true_R2_full = .3,
                                       true_R2_restricted = .2,
                                       p_full = 3, p_restricted = 3))
  expect_s3_class(counts, "dmar_argument_invalid")
  expect_identical(counts$arg, c("p_restricted", "p_full"))

  two_ways <- .caught(ss_power_R2_change(true_R2_full = .3,
                                         true_R2_restricted = .2,
                                         true_f2 = .1, p_full = 5,
                                         p_restricted = 3))
  expect_s3_class(two_ways, "dmar_argument_invalid")
  expect_match(conditionMessage(two_ways), "exactly one way")
})

test_that("the count of predictors in the partial family is refused by class below 2", {
  err <- .caught(var_partial_r(.3, 100, 1))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "p")
  expect_identical(err$value, 1)
  expect_s3_class(.caught(ss_aipe_partial_r(true_cor = .3, p = 1, width = .2)),
                  "dmar_argument_invalid")
  expect_s3_class(.caught(expected_partial_r(.3, 100, c(2, 1))),
                  "dmar_argument_invalid")
})

test_that("smd() from summary statistics refuses a group of one with its own standard deviation", {
  for (sizes in list(c(1, 2), c(2, 1), c(1, 1))) {
    err <- .caught(smd(mean_1 = 10, mean_2 = 8, s_1 = 2, s_2 = 3,
                       n_1 = sizes[1], n_2 = sizes[2]))
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, c("n_1", "n_2"))
    expect_match(conditionMessage(err), "no standard deviation")
  }
  # Two per group is the smallest admissible case and is unchanged.
  ok <- smd(mean_1 = 10, mean_2 = 8, s_1 = 2, s_2 = 3, n_1 = 2, n_2 = 2)
  expect_equal(ok$value, 2 / sqrt((4 + 9) / 2), tolerance = 1e-12)
  # The raw-data path refused already.
  expect_error(smd(group_1 = c(1, 2), group_2 = 3), "1 individual")
})

test_that("unbiased = TRUE without the sample sizes is refused by class", {
  err <- .caught(smd_c(mean_T = 10, mean_C = 8, s_C = 2, unbiased = TRUE))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "n_C")
  err <- .caught(smd(mean_1 = 10, mean_2 = 8, s = 2, unbiased = TRUE))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, c("n_1", "n_2"))
})

test_that("the Hedges refusal gives the reason that holds at its degrees of freedom", {
  one <- .caught(smd_c(mean_T = 1, mean_C = 0, s_C = 1, n_C = 2,
                       unbiased = TRUE))
  expect_s3_class(one, "dmar_hedges_correction_undefined")
  expect_identical(one$df, 1)
  expect_match(conditionMessage(one), "is 0 at df = 1", fixed = TRUE)
  expect_match(conditionMessage(one), "on 1 degree of freedom", fixed = TRUE)

  none <- .caught(smd_c(mean_T = 1, mean_C = 0, s_C = 1, n_C = 1,
                        unbiased = TRUE))
  expect_s3_class(none, "dmar_hedges_correction_undefined")
  expect_identical(none$df, 0)
  expect_match(conditionMessage(none), "is not defined at df = 0", fixed = TRUE)
  expect_no_match(conditionMessage(none), "is 0 at df = 1", fixed = TRUE)

  pooled <- .caught(smd(mean_1 = 1, mean_2 = 0, s = 1, n_1 = 1, n_2 = 1,
                        unbiased = TRUE))
  expect_s3_class(pooled, "dmar_hedges_correction_undefined")
  expect_identical(pooled$df, 0)
})

# ---- Warnings that joined a class -------------------------------------------

test_that("the search warnings of ss_aipe_sm() and ss_aipe_R2() have a class and a count", {
  # ss_aipe_sm() counted the warnings ci_nc_t() gave past a noncentrality
  # of 37.62 and gave one; the limits are now exact there and no warning
  # is given. The plan is unchanged.
  seen <- list()
  res <- withCallingHandlers(
    ss_aipe_sm(true_sm = 10, width = 2.5, assurance = .9),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(seen, 0L)
  expect_equal(res$value, 150)

  seen <- list()
  res <- withCallingHandlers(
    ss_aipe_R2(true_R2 = .3, p = 3, width = .2),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "dmar_nc_F_clamp_search")
  expect_gt(seen[[1]]$n_evaluations, 0)
  expect_equal(res$value, 226)
})

test_that("mediation_mbco() stitches its combined Monte Carlo warning as sentences about draws", {
  settle <- function(B, used, term) {
    .dmar_boot_settle(
      B, used, tail = 0.025, conf_level = 0.95, fun = "mediation_mbco",
      quantity = paste0("the effect '", term, "'"),
      dropped_because = paste0("gave a value of '", term,
                               "' that is not a finite number"),
      replications = "Monte Carlo draws", signal = FALSE)
  }
  w <- tryCatch(
    .mbco_mc_settle(list(ab = settle(2000, 1990, "ab"),
                         total = settle(2000, 900, "total"))),
    warning = function(w) w)
  expect_s3_class(w, "dmar_bootstrap_insufficient")
  msg <- conditionMessage(w)
  expect_match(msg, "follows. No interval is reported.", fixed = TRUE)
  expect_no_match(msg, "follows. no interval", fixed = TRUE)
  # Draws of the estimates are not resamples of the data.
  expect_no_match(msg, "a minority of the resamples", fixed = TRUE)
  expect_no_match(msg, "too sparse", fixed = TRUE)
  expect_match(msg, "a minority of the draws", fixed = TRUE)
  expect_match(msg, "not resamples of the data", fixed = TRUE)
})

test_that("the message of mediation_mbco() for one effect says draws throughout", {
  settle <- function(B, used) {
    .dmar_boot_settle(
      B, used, tail = 0.025, conf_level = 0.95, fun = "mediation_mbco",
      quantity = "the effect 'ab'",
      dropped_because = "gave a value of 'ab' that is not a finite number",
      replications = "Monte Carlo draws", signal = FALSE)
  }
  # Draws dropped, with the caution (5 percent of the draws, past the 2.5
  # percent beyond each limit).
  w <- tryCatch(.mbco_mc_settle(list(ab = settle(2000, 1900))),
                warning = function(w) w)
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  msg <- conditionMessage(w)
  expect_match(msg, "^mediation_mbco\\(\\): 100 of 2000 Monte Carlo draws")
  expect_match(msg, "describes draws in which the effect 'ab' is defined",
               fixed = TRUE)
  expect_match(msg, paste0("Caution: the dropped draws are at least the ",
                           "2.5% of the distribution of the draws"),
               fixed = TRUE)
  expect_match(msg, "where those draws would have fallen", fixed = TRUE)
  expect_no_match(msg, "resample|replication|bootstrap distribution")
  # A run below the required number of draws, marked as exploratory.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  e <- tryCatch(.mbco_mc_settle(list(ab = settle(200, 200))),
                warning = function(w) w)
  expect_s3_class(e, "dmar_bootstrap_exploratory")
  expect_match(conditionMessage(e), "B = 200 is below the 400 draws a 95%",
               fixed = TRUE)
  expect_match(conditionMessage(e), "With 200 draws each limit rests",
               fixed = TRUE)
  expect_no_match(conditionMessage(e), "resample|replication")
})

test_that("a bootstrap that loses most of its resamples still speaks of resamples and sparse data", {
  verdict <- .dmar_boot_settle(
    400, 150, tail = 0.025, conf_level = 0.95, fun = "cohen_kappa",
    quantity = "kappa", dropped_because = "gave an undefined kappa",
    signal = FALSE)
  msg <- .dmar_boot_message(verdict)
  expect_match(msg, "a minority of the resamples", fixed = TRUE)
  expect_match(msg, "too sparse for this bootstrap", fixed = TRUE)
})

test_that("the counts, the generators, and the coefficient planners take their one name and refuse the 1.0.x names", {
  unused <- "unused argument"
  # Counts are n_<things>.
  expect_error(cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3), unused)
  expect_error(cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6,
                                 covariates = 1), unused)
  expect_error(qbryant_paulson(0.95, num_covariates = 1, num_groups = 6,
                               df = 14), unused)
  expect_error(pbryant_paulson(4.64, num_covariates = 1, num_groups = 6,
                               df = 14), unused)
  expect_error(dbryant_paulson(4.64, num_covariates = 1, num_groups = 6,
                               df = 14), unused)
  adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
  # A 1.0.x name is not dropped through '...': with the per-group size given
  # by position, num_covariates = 2 would have run for one covariate.
  expect_error(ci_c_ancova_bp(adj, sqrt(0.01326), NULL, 4, num_covariates = 2,
                              df = 14), unused)
  # suppressWarnings: under warnPartialMatchArgs (tests/testthat/setup-options.R)
  # R warns of the partial match of the 1.0.x name before it stops; the error
  # is the property under test.
  expect_error(suppressWarnings(
    ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326), n = 4, df = 14)),
    "matches multiple formal arguments")
  expect_identical(names(formals(ci_c_ancova_bp))[4:5],
                   c("n_per_group", "n_covariates"))
  expect_error(sem_fit_measures(chisq = 115.85, df = 48, N = 301, groups = 2),
               unused)
  expect_error(equivalence_rmsea(chisq = 115.85, df = 48, N = 301,
                                 rmsea_bound = .08, groups = 2), unused)
  expect_error(ci_eigenvalue(diag(3), N = 50, k = 1), unused)
  expect_error(procrustes_phi(c(.8, .7, .6), c(.7, .7, .5), n_perm = 10),
               unused)
  # ecvi(): the divisor is the multiplier of the fit function. A complete
  # 1.0.x call is an unused argument; without 'npar' the old 'n' reaches
  # 'npar' and the call stops for the multiplier it lacks.
  expect_error(ecvi(chisq = 24.361, df = 8, npar = 13, n = 301), unused)
  expect_error(suppressWarnings(ecvi(chisq = 24.361, df = 8, n = 301)),
               "'multiplier'")
  # The generators: units per population and occasions per unit.
  expect_error(suppressWarnings(simulate_longitudinal_polynomial(
    n = 5, target_times = 0:3, fixed_coefficients = c(10, 1),
    random_variances = c(1, 0.1), error_variance = 1)),
    "matches multiple formal arguments")
  for (fn in list(simulate_longitudinal_polynomial,
                  simulate_longitudinal_logistic,
                  simulate_longitudinal_gompertz,
                  simulate_longitudinal_negative_exponential,
                  simulate_longitudinal_richards)) {
    expect_true(all(c("n_per_group", "n_occasions") %in% names(formals(fn))))
    expect_false(any(c("n", "occasions") %in% names(formals(fn))))
  }
  # Every population value of a coefficient planner carries true_.
  nuisance <- c("true_R2_Y_X", "true_R2_j_X_without_j", "true_sd_Y",
                "true_sd_X", "true_cor_XX", "true_cor_YX")
  for (fn in list(ss_aipe_reg_coef, ss_aipe_rc, ss_aipe_src, ss_power_reg_coef,
                  ss_power_rc)) {
    expect_true(all(nuisance %in% names(formals(fn))))
  }
  expect_error(ss_power_reg_coef(rho2_Y_X = .5, rho2_Y_X_without_j = .4, p = 5),
               unused)
  expect_error(ss_power_rc(cohen_f2 = .2, p = 5), unused)
  expect_error(ss_power_R2(cohen_f2 = .2, p = 5), unused)
  expect_error(ss_aipe_reg_coef(rho2_Y_X = .3, rho2_j_X_without_j = .2, p = 3,
                                true_reg_coef_j = .4, width = .3), unused)
  expect_error(ss_aipe_R2(true_R2 = .5, p = 5, width = .2, g = 10), unused)
  expect_true(all(c("true_cor_YX", "true_cor_XX") %in%
                    names(formals(ss_aipe_R2_sensitivity))))
  expect_true(all(c("true_cor_YX", "true_cor_XX") %in%
                    names(formals(ss_power_R2_sensitivity))))
})

test_that("rows follow the one name: true_f, n_groups, n_studies, n_occasions, n_resamples", {
  coef <- ss_power_reg_coef(true_f2 = .2, p = 5, specified_N = 60)
  expect_identical(coef$term, c("specified_N", "actual_power",
                                "noncentral_t_parm", "true_f"))
  expect_equal(.term(coef, "true_f"), sqrt(.2))
  expect_identical(ss_power_rc(true_f2 = .2, p = 5, specified_N = 60)$term,
                   coef$term)
  fit <- sem_fit_measures(chisq = 115.85, df = 48, N = 301, n_groups = 2)
  expect_equal(.term(fit, "n_groups"), 2)
  expect_false("groups" %in% fit$term)
  eq <- equivalence_rmsea(chisq = 115.85, df = 48, N = 301, rmsea_bound = .08,
                          n_groups = 2)
  expect_equal(.term(eq, "n_groups"), 2)
  d <- c(.2, .5, .3, .4); v <- rep(.04, 4)
  expect_equal(.term(meta_es(d, v), "n_studies"), 4)
  expect_equal(.term(meta_contrast(d, v, weights = c(1, -1, 1, -1)),
                     "n_studies"), 4)
  expect_false("k" %in% meta_es(d, v)$term)
  pcm <- ss_power_pcm(true_std_reg_coef = -.4, true_var_trend = .003,
                      true_var_error = .0262, frequency = 1,
                      duration = 4, desired_power = .8)
  expect_equal(.term(pcm, "n_occasions"), 5)
  phi <- procrustes_phi(c(.8, .7, .6, .5), c(.7, .7, .5, .6),
                        n_resamples = 20, seed = 113)
  expect_equal(.term(phi, "n_resamples"), 20)
})

test_that("the analysis of variance, contrast, and equivalence functions take their one name and refuse the 1.0.x names", {
  unused <- "unused argument"
  # Cohen's f, the repeated measures correlation, and the sphericity value
  # are planning values.
  expect_error(ss_power_one_way_anova(a = 3, f = .25), unused)
  expect_error(ss_power_rm_anova(a = 4, f = .25, true_cor = .5), unused)
  expect_error(ss_power_rm_anova(a = 4, true_f = .25, rho = .5), unused)
  expect_error(ss_power_rm_anova(a = 4, true_f = .25, epsilon = .75), unused)
  expect_error(ss_power_split_plot_anova(a = 2, b = 4, effect = "within",
                                         true_f = .25, rho = .5), unused)
  expect_error(ss_power_factorial_anova(factor_levels = c(2, 3),
                                        effect_indices = 1, f = .25), unused)
  expect_error(ss_power_factorial_ancova(factor_levels = c(2, 3),
                                         effect_indices = 1, f = .25), unused)
  expect_error(ss_power_factorial_ancova(factor_levels = c(2, 3),
                                         effect_indices = 1, true_f = .25,
                                         covariate_R2 = .2, n_covariates = 1),
               unused)
  # With the levels given by position, R matches the 1.0.x 'f' to
  # 'factor_levels'; the value is not a vector of levels and the call stops.
  # suppressWarnings: under warnPartialMatchArgs (tests/testthat/setup-options.R)
  # R warns of the partial match of the 1.0.x name before it stops; the error
  # is the property under test.
  expect_error(suppressWarnings(ss_power_factorial_anova(c(2, 3), 1, f = .25)),
               "'factor_levels'")
  expect_error(suppressWarnings(ss_power_factorial_ancova(c(2, 3), 1, f = .25)),
               "'factor_levels'")
  # The contrast planners.
  expect_error(ss_power_c(true_contrast = .5, c_weights = c(1, -1), sigma = 1),
               unused)
  expect_error(ss_power_contrast(c_weights = c(1, -1), mu = c(1, 2),
                                 true_var_error = 1), unused)
  expect_error(ss_power_contrast(c_weights = c(1, -1), true_means = c(1, 2),
                                 sigma_squared = 1), unused)
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), sigma = 1,
                                      contrast_upper = .5), unused)
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = 20,
                                   sigma = 1, contrast_upper = .5), unused)
  expect_error(ss_aipe_sc_ancova(true_std_contrast = .8, ratio = .6,
                                 width = .5, c_weights = c(1, -1)), unused)
  expect_error(ss_aipe_sc_ancova(true_contrast = 1, sigma_anova = 2,
                                 sigma_ancova = 1.5, divisor = "s_anova",
                                 width = .5, c_weights = c(1, -1)), unused)
  # ss_aipe_c(): the error variance is given one way.
  expect_identical(names(formals(ss_aipe_c)),
                   c("true_var_error", "c_weights", "width", "conf_level",
                     "assurance", "true_kurtosis"))
  expect_error(ss_aipe_c(MSwithin = 40, c_weights = c(1, -1), width = 3),
               unused)
  expect_error(ss_aipe_c(SD = 6, c_weights = c(1, -1), width = 3), unused)
  expect_error(ss_aipe_c(true_var_error = 40, c_weights = c(1, -1),
                         width = 3, assurance = .9, kurtosis = 1), unused)
  expect_error(ss_aipe_c(c_weights = c(1, -1), width = 3),
               class = "dmar_argument_invalid")
  expect_error(ss_aipe_c(true_var_error = NULL, c_weights = c(1, -1),
                         width = 3), "its square")
  expect_error(ss_aipe_c(true_var_error = -1, c_weights = c(1, -1),
                         width = 3), class = "dmar_argument_invalid")
  # A per-group size is n_per_group, a per-cell size n_per_cell, and the
  # size of the two one-sided tests functions n_per_treatment; no function
  # of the group keeps an argument 'n'.
  per_group <- list(ci_c, ci_sc, ci_c_ancova, ci_sc_ancova, equivalence_c,
                    cohen_f, is_orthogonal_set, variance_components_mls,
                    simulate_anova_data, simulate_ancova_data, ss_seq_c,
                    mr_smd, power_equivalence_c)
  for (fn in per_group) {
    expect_true("n_per_group" %in% names(formals(fn)))
    expect_false("n" %in% names(formals(fn)))
  }
  # A per-group size given so that the power at that size is returned is
  # specified_n_per_group.
  for (fn in list(ss_power_c, ss_power_sc, ss_power_contrast,
                  ss_power_split_plot_anova, ss_power_c_ancova,
                  ss_power_composite_ancova_2group)) {
    expect_true("specified_n_per_group" %in% names(formals(fn)))
    expect_false(any(c("n", "n_per_group") %in% names(formals(fn))))
  }
  expect_true("n_per_cell" %in% names(formals(simulate_ancova_factorial_data)))
  for (fn in list(power_equivalence_md, power_density_equivalence_md,
                  power_equivalence_md_plot)) {
    expect_true("n_per_treatment" %in% names(formals(fn)))
    expect_false("n" %in% names(formals(fn)))
  }
  # Where the 1.0.x 'n' matches two arguments, the call stops.
  multiple <- "matches multiple formal arguments"
  # suppressWarnings: under warnPartialMatchArgs (tests/testthat/setup-options.R)
  # R warns of the partial match of the 1.0.x name before it stops; the error
  # is the property under test.
  expect_error(suppressWarnings(
    ci_sc(means = c(1, 2, 3), s_anova = 1, c_weights = c(1, -1, 0), n = 10)),
    multiple)
  # In the equivalence power functions the degrees of freedom are df_error,
  # so 'n' matches n_per_treatment alone, the quantity it always was: a 1.0.x
  # call that gave the degrees of freedom by position runs as before.
  expect_identical(
    suppressWarnings(power_equivalence_md(.05, FALSE, -.2, .2, .05, .20,
                                          n = 24, 22))$value,
    power_equivalence_md(.05, FALSE, -.2, .2, .05, .20,
                         n_per_treatment = 24, df_error = 22)$value)
  expect_error(suppressWarnings(simulate_ancova_factorial_data(
    a = 2, b = 2, true_mean_Y = c(10, 12, 11, 14), true_mean_X = 5, true_sd_Y = 2,
    true_sd_X = 1, true_cor_YX = .4, n = 10)), multiple)
})

test_that("power_equivalence_c() reads n_per_group for independent groups and N for a paired design", {
  ind <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 24,
                             true_sd = .20, contrast_upper = .2,
                             true_contrast = .05, df_error = 22)
  expect_equal(ind$value, 0.8029678, tolerance = 1e-6)
  pai <- power_equivalence_c(N = 50, true_sd = 1, contrast_lower = .5,
                             side = "noninferiority", design = "paired")
  expect_equal(pai$value,
               pt(qt(.95, 49), 49, ncp = .5 * sqrt(50), lower.tail = FALSE),
               tolerance = 1e-12)
  # The size of the other design is refused, with a class and the names.
  e <- tryCatch(power_equivalence_c(c_weights = c(1, -1), N = 50,
                                    true_sd = 1, contrast_upper = .5),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, c("N", "design"))
  expect_match(conditionMessage(e), "'n_per_group'", fixed = TRUE)
  e <- tryCatch(power_equivalence_c(n_per_group = 50, true_sd = 1,
                                    contrast_lower = .5, side = "noninferiority",
                                    design = "paired"),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, c("n_per_group", "design"))
  expect_match(conditionMessage(e), "'N'", fixed = TRUE)
  # A missing size is refused for the design that needs it.
  expect_error(power_equivalence_c(c_weights = c(1, -1), true_sd = 1,
                                   contrast_upper = .5),
               "'n_per_group' must be given", class = "dmar_argument_invalid")
  expect_error(power_equivalence_c(true_sd = 1, contrast_lower = .5,
                                   side = "noninferiority", design = "paired"),
               "'N' must be the number of units",
               class = "dmar_argument_invalid")
  expect_error(power_equivalence_c(N = 1, true_sd = 1, contrast_lower = .5,
                                   side = "noninferiority", design = "paired"),
               "'N' must be the number of units")
  expect_error(power_equivalence_c(N = c(20, 20), true_sd = 1,
                                   contrast_lower = .5, side = "noninferiority",
                                   design = "paired"),
               "'N' must be the number of units")
})

test_that("rows follow the one name: true_f, noncentral_f_parm, df_effect, df_error, necessary_N, pilot sizes", {
  one <- ss_power_one_way_anova(a = 3, true_f = .25, specified_N = 60)
  expect_identical(one$term, c("specified_N", "a", "noncentral_f_parm",
                               "actual_power"))
  expect_equal(.term(one, "noncentral_f_parm"), 60 * .25^2)
  rm_plan <- ss_power_rm_anova(a = 4, true_f = .25, true_cor = .5,
                               specified_N = 20)
  expect_identical(rm_plan$term, c("specified_N", "a", "df_effect", "df_error",
                                   "noncentral_f_parm", "actual_power"))
  sp <- ss_power_split_plot_anova(a = 2, b = 4, effect = "within",
                                  true_f = .25, true_cor = .5, specified_n_per_group = 20)
  expect_identical(sp$term, c("specified_n_per_group", "total_N", "df_effect",
                              "df_error", "noncentral_f_parm", "actual_power"))
  fa <- ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 2,
                                 true_f = .25, specified_n_per_cell = 20)
  expect_identical(fa$term, c("specified_n_per_cell", "total_N", "df_effect",
                              "df_error", "noncentral_f_parm", "actual_power"))
  fc <- ss_power_factorial_ancova(factor_levels = c(2, 3), effect_indices = 2,
                                  true_f = .25, true_covariate_R2 = .25,
                                  n_covariates = 1, specified_n_per_cell = 20)
  expect_identical(fc$term, c("specified_n_per_cell", "total_N", "actual_power",
                              "df_effect", "df_error", "true_f",
                              "true_f_adjusted", "true_covariate_R2",
                              "n_covariates", "noncentral_f_parm",
                              "alpha_level"))
  expect_equal(.term(fc, "true_f_adjusted"), .25 / sqrt(.75))
  expect_equal(.term(fc, "noncentral_f_parm"),
               .term(fc, "total_N") * .term(fc, "true_f_adjusted")^2)
  con <- ss_power_contrast(c_weights = c(1, -1), true_means = c(10, 12),
                           true_var_error = 16, specified_n_per_group = 20)
  expect_identical(con$term, c("specified_n_per_group", "total_N",
                               "actual_power", "noncentral_t_parm", "true_f"))
  expect_equal(.term(con, "true_f"),
               abs(.term(con, "noncentral_t_parm")) / sqrt(40))
  expect_identical(
    mr_smd(structural_cost = 0.0001, epsilon = 0.05, d = 0.5, n_per_group = 30,
           sampling_cost = 0.001)$term, c("risk", "n_1", "n_2", "d"))
  expect_identical(mr_smd(pilot = TRUE, A = 10000, sampling_cost = 2.4)$term,
                   "pilot_n_per_group")
  expect_identical(mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75)$term,
                   "pilot_N")
  set.seed(113)
  cv <- mr_cv(data = rlnorm(60, 1, .3), A = 1, sampling_cost = .001,
              gamma = .49, m0 = 60)
  expect_identical(cv$term, c("risk", "N", "cv"))
  expect_equal(.term(cv, "N"), 60)
  # The four contrast sensitivity functions report pct_ci_less_width.
  s <- ss_aipe_c_sensitivity(true_var_error = 4,
                             supposed_var_error = 4,
                             c_weights = c(1, -1), width = 2, G = 20,
                             seed = 113)
  expect_true("pct_ci_less_width" %in% s$term)
  expect_false("pct_ci_less_w" %in% s$term)
})

test_that("the multilevel, intraclass correlation, and reliability functions take their one name and refuse the 1.0.x names", {
  has <- function(f, new, old) {
    fm <- names(formals(f))
    expect_true(all(new %in% fm))
    expect_false(any(old %in% fm))
  }
  has(ss_aipe_icc, c("true_icc", "n_raters"), c("rho", "true_rho", "k"))
  has(ss_aipe_icc_sensitivity, c("true_icc", "supposed_icc", "n_raters"),
      c("true_rho", "estimated_rho", "k"))
  has(var_icc, c("icc", "N", "n_raters"), c("rho", "n", "k"))
  has(var_alpha, c("alpha", "N", "n_items"), c("n", "p_items"))
  has(ss_power_mixed_effects,
      c("cluster_size", "true_icc", "specified_n_clusters_per_arm"), c("n", "rho", "J"))
  has(ss_aipe_reliability,
      c("true_cov_matrix", "n_items", "true_item_covariance", "true_factor_loading",
        "true_residual_var", "G_initial", "G"),
      c("S", "i", "cor_est", "lambda", "psi_square", "initial_iter",
        "final_iter"))
  has(ss_aipe_reliability_sensitivity, "n_items", c("i", "items"))
  # The estimators keep S, a sample covariance matrix.
  for (f in list(reliability, reliability_alpha, reliability_omega)) {
    expect_true("S" %in% names(formals(f)))
  }
  # Each 1.0.x name is an argument R no longer finds. The calls are built
  # with do.call() so that the old name is data here and not a call a
  # rename sweep would change.
  unused <- function(f, args) {
    expect_error(do.call(f, args), "unused argument")
  }
  unused(ss_aipe_icc, list(true_rho = .7, n_raters = 3, width = .2))
  unused(ss_aipe_icc, list(true_icc = .7, k = 3, width = .2))
  unused(ss_aipe_icc_sensitivity,
         list(true_icc = .7, estimated_rho = .7, n_raters = 3, width = .4,
              G = 5))
  unused(var_icc, list(rho = .6, N = 30, n_raters = 4))
  unused(var_icc, list(icc = .6, N = 30, k = 4))
  unused(var_alpha, list(alpha = .8, N = 100, p_items = 10))
  for (old in c("S", "i", "cor_est", "lambda", "psi_square", "initial_iter",
                "final_iter")) {
    unused(ss_aipe_reliability,
           c(list(model = "Parallel", type = "Normal Theory", width = .1),
             stats::setNames(list(1), old)))
  }
  unused(ss_aipe_reliability_sensitivity,
         list(true_reliability = .8, specified_N = 80, i = 4, width = .3,
              G = 4))
})

test_that("rows follow the one name: n_items, n_raters, true_icc, assurance_attained, G", {
  icc_plan <- ss_aipe_icc(true_icc = .7, n_raters = 3, width = .2)
  expect_identical(icc_plan$term, c("necessary_N", "expected_width", "true_icc",
                                    "n_raters", "width", "conf_level"))
  expect_identical(.term(icc_plan, "n_raters"), 3)
  icc_check <- ss_aipe_icc_sensitivity(true_icc = .5, specified_N = 30,
                                       n_raters = 3, width = .4, G = 10,
                                       seed = 113)
  expect_true(all(c("pct_ci_less_width", "total_N", "n_raters", "true_icc",
                    "supposed_icc") %in% icc_check$term))
  expect_false(any(c("pct_ci_less_w", "k", "true_rho", "estimated_rho") %in%
                     icc_check$term))
  # The variance of the intraclass correlation and of alpha, by the
  # formulas of their help pages (var_alpha(): van Zyl et al.'s N - 1
  # form, then the N - 2 form that follows from Bonett).
  expect_equal(var_icc(icc = .6, N = 30, n_raters = 4)$value,
               2 * (1 - .6)^2 * (1 + 3 * .6)^2 / (30 * 4 * 3))
  expect_equal(var_alpha(alpha = .8, N = 100, n_items = 10)$value,
               2 * 10 * (1 - .8)^2 / (9 * c(99, 98)))
  S <- matrix(.5, 4, 4); diag(S) <- 1
  for (r in list(reliability_alpha(S = S, N = 100, ci_method = "none"),
                 reliability(S = S, N = 100, type = "alpha",
                             ci_method = "none"))) {
    expect_identical(r$term[length(r$term)], "n_items")
    expect_identical(.term(r, "n_items"), 4)
    expect_false("J" %in% r$term)
    expect_identical(generics::glance(r)$n_items, 4)
  }
  plan <- ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                              width = .1, n_items = 6,
                              true_item_covariance = .6, true_residual_var = .4,
                              assurance = .85, G_initial = 20, G = 40,
                              seed = 113)
  expect_identical(plan$term, c("necessary_N", "width", "assurance",
                                "assurance_attained", "G"))
  expect_identical(.term(plan, "assurance"), .85)
  expect_identical(.term(plan, "G"), 40)
  expect_gte(.term(plan, "assurance_attained"), .85)
  for (bad in list(2.5, 0, NA_real_, c(5, 6))) {
    expect_error(ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                     width = .1, n_items = 6,
                                     true_item_covariance = .6,
                                     true_residual_var = .4, G = bad),
                 class = "dmar_count_invalid")
    expect_error(ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                     width = .1, n_items = 6,
                                     true_item_covariance = .6,
                                     true_residual_var = .4, G_initial = bad),
                 class = "dmar_count_invalid")
  }
  rel_check <- suppressWarnings(ss_aipe_reliability_sensitivity(
    true_reliability = .8, specified_N = 60, n_items = 4, width = .3, G = 10,
    seed = 113))
  expect_true(all(c("pct_ci_less_width", "n_items") %in% rel_check$term))
  expect_false(any(c("pct_ci_less_w", "items") %in% rel_check$term))
})

test_that("the composite planners, the mean and two-phase planners, and ss_power_c_ancova() take their one name and refuse the 1.0.x names", {
  has <- function(f, new, old) {
    fm <- names(formals(f))
    expect_true(all(new %in% fm))
    expect_false(any(old %in% fm))
  }
  has(ss_power_c_ancova, c("true_sd", "specified_n_per_group"), c("sigma", "n"))
  has(ss_power_composite_ancova_2group,
      c("true_sd", "true_sd_X", "specified_n_per_group"),
      c("sigma", "sd_cov", "n"))
  for (f in list(ss_power_composite_anova,
                 ss_power_composite_factorial_anova)) {
    has(f, c("true_means", "true_sd"), c("means", "sigma"))
  }
  has(ss_power_composite_factorial_ancova,
      c("true_means", "true_sd", "true_covariate_R2"),
      c("means", "sigma", "covariate_R2"))
  has(ss_power_composite_factorial_ancova_het,
      c("true_means", "true_cor", "true_sd", "true_sd_X"),
      c("means", "correlations", "sigma", "sd_cov"))
  has(ss_power_composite_ancova,
      c("true_means", "true_sd", "true_covariate_R2", "true_cor",
        "true_sd_X"),
      c("means", "sigma", "covariate_R2", "correlations", "sd_cov"))
  for (f in list(ss_aipe_composite_sem, ss_power_composite_sem)) {
    has(f, c("true_cov_matrix", "true_model", "true_mean"),
        c("Sigma", "pop_model", "mu"))
  }
  has(ss_aipe_rmsea_sensitivity, "true_cov_matrix", "Sigma")
  for (f in list(ss_aipe_mean, ss_aipe_mean_sensitivity, ss_aipe_two_phase,
                 ss_aipe_two_phase_sensitivity)) {
    has(f, "true_kurtosis", "kurtosis")
  }
  # Each 1.0.x name is an argument R no longer finds. The calls are built
  # with do.call() so that the old name is data here and not a call a
  # rename sweep would change.
  unused <- function(f, args) {
    expect_error(do.call(f, args), "unused argument")
  }
  eff <- list(list(factors = 1, true_f = .25))
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  unused(ss_power_c_ancova,
         list(true_contrast = .5, c_weights = c(1, -1), sigma = 1, true_cor = .3,
              n_per_group = 30))
  unused(ss_power_composite_ancova_2group,
         list(true_smd = .5, true_cor = .3, sigma = 1, n_per_group = 30))
  unused(ss_power_composite_ancova_2group,
         list(true_smd = .5, true_cor = .3, sd_cov = 2, n_per_group = 30))
  for (f in list(ss_power_composite_anova, ss_power_composite_factorial_anova,
                 ss_power_composite_factorial_ancova,
                 ss_power_composite_ancova)) {
    unused(f, list(factor_levels = c(2, 3), effects = eff, means = M,
                   true_sd = 4, specified_n_per_cell = 20))
    unused(f, list(factor_levels = c(2, 3), effects = eff, true_means = M,
                   sigma = 4, specified_n_per_cell = 20))
  }
  unused(ss_power_composite_factorial_ancova,
         list(factor_levels = c(2, 3), effects = eff, covariate_R2 = .2,
              n_covariates = 1, specified_n_per_cell = 20))
  unused(ss_power_composite_factorial_ancova_het,
         list(factor_levels = 2, effects = list(list(type = "covariate")),
              correlations = c(.3, .3), true_sd = 1, specified_n_per_cell = 20))
  unused(ss_power_composite_factorial_ancova_het,
         list(factor_levels = 2, effects = list(list(type = "covariate")),
              true_cor = c(.3, .3), true_sd = 1, sd_cov = 2,
              specified_n_per_cell = 20))
  unused(ss_aipe_mean, list(true_sd = 1, width = .4, assurance = .9,
                            kurtosis = 1))
  unused(ss_aipe_two_phase, list(true_sd_Y = 2, true_cor = .7, width = .6,
                                 assurance = .9, kurtosis = 1))
  # The 1.0.x 'n' of the per-group planners is no longer matched by R to
  # any argument (the size is 'specified_n_per_group'): the call stops.
  unused(ss_power_c_ancova,
         list(true_contrast = .5, c_weights = c(1, -1), true_sd = 1,
              true_cor = .3, n = 30))
  # ss_power_c_ancova() refuses a per-group size that is not one finite
  # number of at least 2, by class.
  for (bad in list(1, Inf, NA_real_, c(20, 30), "30")) {
    expect_error(ss_power_c_ancova(true_contrast = .5, c_weights = c(1, -1),
                                   true_sd = 1, true_cor = .3,
                                   specified_n_per_group = bad),
                 class = "dmar_argument_invalid")
  }
})

test_that("a field of an effect that the composite planners do not read is refused", {
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  # The 1.0.x field names, with effect sizes and with population means
  # (where an unread field was dropped without a word).
  old_f <- list(list(factors = 1, f = .25))
  old_pes <- list(list(factors = 1, partial_eta_squared = .06))
  for (f in list(ss_power_composite_anova, ss_power_composite_factorial_anova,
                 ss_power_composite_factorial_ancova,
                 ss_power_composite_ancova)) {
    e <- expect_error(f(factor_levels = c(2, 3), effects = old_f,
                        specified_n_per_cell = 20),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, "effects")
    expect_identical(e$fields, "f")
    expect_match(conditionMessage(e), "'true_f'", fixed = TRUE)
    expect_error(f(factor_levels = c(2, 3), effects = old_pes,
                   specified_n_per_cell = 20),
                 class = "dmar_argument_invalid")
    expect_error(f(factor_levels = c(2, 3), effects = old_f, true_means = M,
                   true_sd = 4, specified_n_per_cell = 20),
                 class = "dmar_argument_invalid")
  }
  e <- expect_error(ss_power_composite_factorial_ancova_het(
    factor_levels = 2,
    effects = list(list(type = "covariate", true_f = .3),
                   list(type = "slope", factors = 1, f = .2)),
    specified_n_per_cell = 20), class = "dmar_argument_invalid")
  expect_identical(e$effect, 2L)
  # An element with no name is refused too.
  expect_error(ss_power_composite_anova(
    factor_levels = c(2, 3), effects = list(list(1, true_f = .25)),
    specified_n_per_cell = 20), class = "dmar_argument_invalid")
  # The fields that are read are accepted, 'label' and 'type' among them.
  expect_no_error(ss_power_composite_anova(
    factor_levels = c(2, 3),
    effects = list(list(factors = 1, true_f = .25, label = "A"),
                   list(factors = 2, true_eta_squared_partial = .05)),
    specified_n_per_cell = 20))
  expect_no_error(ss_power_composite_factorial_ancova_het(
    factor_levels = 2,
    effects = list(list(type = "mean", factors = 1, true_f = .25),
                   list(type = "covariate", true_f = .30, label = "X")),
    specified_n_per_cell = 20))
})

test_that("rows follow the one name: true_f_<label>, noncentral_f_parm_<label>, df_error, total_N, true_sd_X, true_kurtosis, pct_ci_less_width", {
  plan <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 3),
    effects = list(list(factors = 1, true_f = .25, label = "A"),
                   list(factors = 2, true_f = .20, label = "B")),
    true_covariate_R2 = .25, n_covariates = 1, desired_power = .80)
  expect_identical(plan$term,
                   c("necessary_n_per_cell", "total_N", "composite_power",
                     "df_error", "power_A", "power_B", "true_f_A", "true_f_B",
                     "df_A", "df_B", "noncentral_f_parm_A",
                     "noncentral_f_parm_B", "true_covariate_R2",
                     "n_covariates", "cells", "alpha_level", "desired_power"))
  expect_identical(.term(plan, "total_N"),
                   6 * .term(plan, "necessary_n_per_cell"))
  expect_identical(.term(plan, "df_error"), .term(plan, "total_N") - 6 - 1)
  expect_identical(generics::tidy(plan)$estimate,
                   .term(plan, "necessary_n_per_cell"))
  at_n <- ss_power_composite_anova(
    factor_levels = c(2, 3),
    effects = list(list(factors = 1, true_f = .25)), specified_n_per_cell = 20)
  expect_identical(at_n$term[1:2], c("specified_n_per_cell", "total_N"))
  expect_identical(.term(at_n, "total_N"), 120)
  expect_false(any(c("necessary_N", "specified_N", "residual_df", "f_1",
                     "noncentral_parm_1", "true_covariate_R2") %in%
                     at_n$term))
  # The population values travel as attributes under the arguments' names.
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  from_means <- ss_power_composite_anova(
    factor_levels = c(2, 3), true_means = M, true_sd = 4,
    effects = list(list(factors = 1, label = "A")), specified_n_per_cell = 30)
  expect_equal(attr(from_means, "true_means"), M, ignore_attr = TRUE)
  expect_identical(attr(from_means, "true_sd"), 4)
  expect_null(attr(from_means, "means"))
  two <- ss_power_composite_ancova_2group(true_smd = .5, true_cor = .3,
                                          true_sd = 2, true_sd_X = 3,
                                          specified_n_per_group = 40)
  expect_true(all(c("specified_n_per_group", "total_N", "df_error",
                    "true_sd", "true_sd_X") %in% two$term))
  expect_identical(.term(two, "total_N"), 80)
  expect_false(any(c("specified_N", "residual_df", "sigma", "sd_cov") %in%
                     two$term))
  # Unequal absolute correlations: the planned total is relabeled with the
  # planned per-group size, and a total beside a specified size is not.
  planned <- suppressWarnings(ss_power_composite_ancova_2group(
    true_smd = .5, true_cor = c(.1, .5),
    composite_terms = c("group", "group_by_covariate"), desired_power = .80))
  expect_identical(planned$term[1:2],
                   c("approximate_n_per_group", "approximate_total_N"))
  expect_identical(generics::tidy(planned)$estimate,
                   .term(planned, "approximate_n_per_group"))
  given <- suppressWarnings(ss_power_composite_ancova_2group(
    true_smd = .5, true_cor = c(.1, .5), specified_n_per_group = 40))
  expect_identical(given$term[1:2], c("specified_n_per_group", "total_N"))
  # The planning value of the kurtosis is echoed under its name.
  mean_plan <- ss_aipe_mean(true_sd = 1, width = .4, assurance = .9,
                            true_kurtosis = 1)
  expect_identical(.term(mean_plan, "true_kurtosis"), 1)
  expect_false("kurtosis" %in% mean_plan$term)
  two_phase <- ss_aipe_two_phase(true_sd_Y = 2, true_cor = .7, width = .6,
                                 assurance = .9, true_kurtosis = -1)
  expect_identical(.term(two_phase, "true_kurtosis"), -1)
  expect_false("kurtosis" %in% two_phase$term)
  # The proportion of intervals no wider than the width, in the registry
  # the sensitivity family is checked against and in design_consequences().
  expect_true("pct_ci_less_width" %in% DMAR:::.SS_AIPE_SENS_CORE_TERMS)
  dc <- design_consequences(true_effect = .2, se = .1, width = .3)
  expect_true("pct_ci_less_width" %in% dc$term)
  expect_false("pct_ci_less_w" %in% dc$term)
})

test_that("the smaller pairs of names are one name each", {
  # The number of processes is ncpus, as in the bootstrap functions.
  crd_es <- c("ss_aipe_crd_es_n_clusters_fixed_width", "ss_aipe_crd_es_n_individuals_fixed_width",
              "ss_aipe_crd_es_n_clusters_fixed_budget", "ss_aipe_crd_es_n_individuals_fixed_budget",
              "ss_aipe_crd_es_both_fixed_budget", "ss_aipe_crd_es_both_fixed_width")
  for (fn in crd_es) {
    f <- get(fn, asNamespace("DMAR"))
    expect_true("ncpus" %in% names(formals(f)), info = fn)
    expect_identical(formals(f)$ncpus, 1, info = fn)
    expect_false(any(c("multicore", "num_proc") %in% names(formals(f))), info = fn)
  }
  expect_error(
    ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20, cluster_cost = 0,
                                           individual_cost = 1, ncpus = 1.5),
    class = "dmar_count_invalid")
  expect_error(
    do.call(ss_aipe_crd_es_n_clusters_fixed_budget,
            list(budget = 1000, cluster_size = 20, cluster_cost = 0, individual_cost = 1,
                 multicore = TRUE)),
    "unused argument")
  # Progress printing in the planners is print_iter.
  for (f in list(ss_power_R2, ss_power_reg_coef, ss_power_rc, ss_power_contrast)) {
    expect_true("print_iter" %in% names(formals(f)))
    expect_false("print_progress" %in% names(formals(f)))
  }
  # The error variances of the reliability planner are true_residual_var.
  expect_true("true_residual_var" %in% names(formals(ss_aipe_reliability)))
  expect_false("true_psi_square" %in% names(formals(ss_aipe_reliability)))
  # Two proportions are p_1 and p_2.
  expect_identical(names(formals(cohen_h)), c("p_1", "p_2"))
  expect_error(do.call(cohen_h, list(p1 = .5, p2 = .4)), "unused argument")
  # The observed statistic of ci_nc_t() has one name.
  expect_false("t_value" %in% names(formals(ci_nc_t)))
  # The total variance of a two-level plan is true_total_var.
  expect_true("true_total_var" %in% names(formals(ss_aipe_mixed_effects)))
  expect_true(all(c("true_total_var", "supposed_total_var") %in%
                    names(formals(ss_aipe_mixed_effects_sensitivity))))
  expect_true("true_eta_squared_partial" %in% names(formals(ss_power_factorial_anova)))
})

test_that("the crd effect size planners give the same width for one process and for two", {
  skip_on_cran()
  skip_on_os("windows")
  one <- ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
    cluster_cost = 0, individual_cost = 1, true_es = 0.5, es_type = 1, true_icc = 0.25,
    pr_treat = 0.5, G = 4, seed = 113, ncpus = 1)
  two <- ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
    cluster_cost = 0, individual_cost = 1, true_es = 0.5, es_type = 1, true_icc = 0.25,
    pr_treat = 0.5, G = 4, seed = 113, ncpus = 2)
  expect_equal(one$value, two$value)
})
