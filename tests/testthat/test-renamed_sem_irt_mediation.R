# The structural equation, psychometric, item response, and mediation
# functions take English quantity names in 1.1.0: a population value a
# planner is given carries true_, a planning assumption of a sensitivity
# function carries supposed_, and a model symbol (Sigma, lambda, nu, theta,
# a, b) is named for the quantity it holds. The DMAR 1.0.x names, and the
# names of the first 1.1.0 builds, stop. The calls with an old name are
# built with do.call() so that the old name is data here, not a call a
# rename sweep would change.

.rn_unused <- function(f, args) {
  expect_error(do.call(f, args), "unused argument")
}
.rn_has <- function(f, new, old) {
  fm <- names(formals(f))
  expect_true(all(new %in% fm), info = paste(setdiff(new, fm), collapse = ", "))
  expect_false(any(old %in% fm), info = paste(intersect(old, fm), collapse = ", "))
}

test_that("the mediation planners and the variance of the indirect effect take regression coefficient names", {
  .rn_has(ss_aipe_indirect_effect,
          c("true_std_reg_coef_a", "true_std_reg_coef_b", "max_N"),
          c("a", "b", "n_max", "true_a", "true_b"))
  .rn_has(ss_aipe_indirect_effect_sensitivity,
          c("true_std_reg_coef_a", "true_std_reg_coef_b",
            "supposed_std_reg_coef_a", "supposed_std_reg_coef_b"),
          c("true_a", "true_b", "estimated_a", "estimated_b"))
  .rn_has(ss_power_indirect_effect,
          c("true_std_reg_coef_a", "true_std_reg_coef_b",
            "true_std_reg_coef_direct", "specified_N"),
          c("a", "b", "c_prime", "N", "true_a", "true_b", "true_c_prime"))
  .rn_has(var_indirect_effect,
          c("reg_coef_a", "reg_coef_b", "var_reg_coef_a", "var_reg_coef_b",
            "cov_reg_coef_a_b"),
          c("a", "b", "var_a", "var_b", "cov_ab"))

  .rn_unused(ss_aipe_indirect_effect, list(a = .4, b = .4, width = .2))
  .rn_unused(ss_aipe_indirect_effect, list(true_a = .4, true_b = .4, width = .2))
  .rn_unused(ss_aipe_indirect_effect_sensitivity,
             list(true_a = .4, true_b = .3, estimated_a = .4, estimated_b = .3,
                  width = .4, G = 10))
  .rn_unused(ss_aipe_indirect_effect_sensitivity,
             list(true_std_reg_coef_a = .4, true_std_reg_coef_b = .3,
                  estimated_a = .4, estimated_b = .3, width = .4, G = 10))
  # 1.0.x a = is matched by R to alpha_level (with the partial matching
  # warning the test run turns on), and the call stops for b.
  suppressWarnings(.rn_unused(ss_power_indirect_effect,
                              list(a = .39, b = .39, desired_power = .8)))
  .rn_unused(ss_power_indirect_effect,
             list(true_a = .39, true_b = .39, desired_power = .8))
  .rn_unused(ss_power_indirect_effect,
             list(true_std_reg_coef_a = .39, true_std_reg_coef_b = .39,
                  true_c_prime = .1, desired_power = .8))
  .rn_unused(var_indirect_effect,
             list(a = .4, b = .4, var_a = .02, var_b = .02))
  .rn_unused(var_indirect_effect,
             list(reg_coef_a = .4, reg_coef_b = .4, var_reg_coef_a = .02,
                  var_reg_coef_b = .02, cov_ab = .005))
})

test_that("the mediation planners echo their paths under the new names", {
  plan <- ss_aipe_indirect_effect(true_std_reg_coef_a = .4,
                                  true_std_reg_coef_b = .4, width = .2)
  expect_identical(plan$term,
                   c("necessary_N", "expected_width", "true_std_reg_coef_a",
                     "true_std_reg_coef_b", "true_std_indirect_effect",
                     "width", "conf_level"))
  expect_equal(plan$value[plan$term == "true_std_indirect_effect"], .16)

  pw <- ss_power_indirect_effect(true_std_reg_coef_a = .39,
                                 true_std_reg_coef_b = .39,
                                 true_std_reg_coef_direct = .1,
                                 specified_N = 100)
  expect_identical(pw$term,
                   c("specified_N", "power_ab", "power_a", "power_b",
                     "true_std_reg_coef_a", "true_std_reg_coef_b",
                     "true_std_reg_coef_direct", "true_std_indirect_effect",
                     "alpha_level"))
  expect_equal(pw$value[pw$term == "true_std_reg_coef_direct"], .1)

  sens <- ss_aipe_indirect_effect_sensitivity(
    true_std_reg_coef_a = .4, true_std_reg_coef_b = .3,
    supposed_std_reg_coef_a = .4, supposed_std_reg_coef_b = .3,
    width = .4, G = 10, seed = 113)
  expect_true(all(c("mean_std_indirect_effect", "median_std_indirect_effect",
                    "sd_std_indirect_effect", "true_std_reg_coef_a",
                    "true_std_reg_coef_b", "true_std_indirect_effect",
                    "supposed_std_reg_coef_a", "supposed_std_reg_coef_b")
                  %in% sens$term))
  expect_false(any(c("mean_ab", "median_ab", "sd_ab", "true_a", "true_b",
                     "true_ab", "estimated_a", "estimated_b") %in% sens$term))
  expect_equal(sens$value[sens$term == "supposed_std_reg_coef_b"], .3)
  # With specified_N the supposed paths echo as NA.
  at_n <- ss_aipe_indirect_effect_sensitivity(
    true_std_reg_coef_a = .4, true_std_reg_coef_b = .3, specified_N = 60,
    width = .4, G = 10, seed = 113)
  expect_true(is.na(at_n$value[at_n$term == "supposed_std_reg_coef_a"]))
  # The per-replication file names the sample indirect effect the same way.
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  ss_aipe_indirect_effect_sensitivity(
    true_std_reg_coef_a = .4, true_std_reg_coef_b = .3, specified_N = 60,
    width = .4, G = 10, seed = 113, filename = f)
  expect_identical(names(utils::read.csv(f))[1], "std_indirect_effect")
})

test_that("the structural equation planners take true_cov_matrix, true_mean, and supposed_cov_matrix", {
  for (f in list(ss_aipe_composite_sem, ss_power_composite_sem)) {
    .rn_has(f, c("true_cov_matrix", "true_model", "true_mean"),
            c("Sigma", "pop_model", "mu", "true_Sigma", "true_mu"))
  }
  .rn_has(ss_aipe_sem_path, c("true_cov_matrix", "width"),
          c("Sigma", "desired_width", "true_Sigma"))
  .rn_has(ss_aipe_sem_path_sensitivity,
          c("supposed_cov_matrix", "true_cov_matrix", "width", "specified_N"),
          c("est_Sigma", "estimated_Sigma", "true_Sigma", "desired_width", "N"))
  .rn_has(ss_aipe_rmsea_sensitivity, c("true_cov_matrix", "specified_N"),
          c("Sigma", "true_Sigma", "N"))
  expect_identical(eval(formals(ss_aipe_sem_path_sensitivity)$true_cov_matrix,
                        list(supposed_cov_matrix = "same")), "same")
})

test_that("an old structural equation name in ... is refused by name before anything is fitted", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  model <- "f =~ x1 + x2 + x3 + x4"
  S <- matrix(0.5, 4, 4, dimnames = list(paste0("x", 1:4), paste0("x", 1:4)))
  diag(S) <- 1
  refused <- function(f, args, bad) {
    err <- tryCatch(do.call(f, args), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, bad)
    expect_match(conditionMessage(err), "lavaan::lavOptions", fixed = TRUE)
  }
  refused(ss_aipe_composite_sem,
          list(model = model, true_Sigma = S, width = .3), "true_Sigma")
  refused(ss_aipe_composite_sem,
          list(model = model, true_cov_matrix = S, width = .3,
               true_mu = rep(0, 4)), "true_mu")
  refused(ss_power_composite_sem,
          list(model = model, true_Sigma = S), "true_Sigma")
  refused(ss_power_composite_sem,
          list(model = model, true_cov_matrix = S, true_mu = rep(0, 4)),
          "true_mu")
  refused(ss_aipe_sem_path,
          list(model = model, true_Sigma = S, width = .3, which_path = "f=~x2"),
          "true_Sigma")
  refused(ss_aipe_sem_path_sensitivity,
          list(model = model, estimated_Sigma = S, which_path = "f=~x2",
               width = .3, G = 2), "estimated_Sigma")
  refused(ss_aipe_sem_path_sensitivity,
          list(model = model, supposed_cov_matrix = S, true_Sigma = S,
               which_path = "f=~x2", width = .3, G = 2), "true_Sigma")
  refused(ss_aipe_rmsea_sensitivity,
          list(width = .05, model = model, true_Sigma = S, G = 2),
          "true_Sigma")
})

test_that("the reliability planners take true_cov_matrix, true_factor_loading, true_residual_var, and supposed_reliability", {
  .rn_has(ss_aipe_reliability,
          c("true_cov_matrix", "true_factor_loading", "true_residual_var"),
          c("S", "lambda", "psi_square", "true_Sigma", "true_lambda",
            "true_psi_squared"))
  .rn_has(ss_aipe_reliability_sensitivity, "supposed_reliability",
          "estimated_reliability")
  base <- list(model = "Congeneric", type = "Factor Analytic", width = .15,
               n_items = 4)
  .rn_unused(ss_aipe_reliability,
             c(base, list(true_lambda = c(.8, .7, .7, .8),
                          true_psi_squared = c(.4, .5, .5, .4))))
  .rn_unused(ss_aipe_reliability,
             c(base, list(lambda = c(.8, .7, .7, .8),
                          psi_square = c(.4, .5, .5, .4))))
  .rn_unused(ss_aipe_reliability,
             list(model = "True Score", type = "Normal Theory", width = .15,
                  true_Sigma = diag(4)))
  .rn_unused(ss_aipe_reliability_sensitivity,
             list(true_reliability = .8, estimated_reliability = .8,
                  n_items = 4, width = .15, G = 2))
  new <- ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                             width = .15, n_items = 4,
                             true_factor_loading = c(.8, .7, .7, .8),
                             true_residual_var = c(.4, .5, .5, .4))
  expect_identical(new$term, "necessary_N")
  sens <- ss_aipe_reliability_sensitivity(true_reliability = .8,
                                          supposed_reliability = .8,
                                          n_items = 4, width = .3, G = 3,
                                          seed = 113)
  expect_true("supposed_reliability" %in% sens$term)
  expect_false("estimated_reliability" %in% sens$term)
  expect_equal(sens$value[sens$term == "supposed_reliability"], .8)
})

test_that("dmacs() takes and reports factor loadings and intercepts", {
  .rn_has(dmacs,
          c("factor_loading_reference", "factor_loading_focal",
            "intercept_reference", "intercept_focal"),
          c("lambda_reference", "lambda_focal", "nu_reference", "nu_focal"))
  .rn_unused(dmacs, list(lambda_reference = c(.8, .7),
                         lambda_focal = c(.8, .7),
                         nu_reference = c(2, 2.3), nu_focal = c(2, 2),
                         sd_pooled = 1))
  res <- dmacs(factor_loading_reference = c(.8, .7),
               factor_loading_focal = c(.8, .7),
               intercept_reference = c(2, 2.3), intercept_focal = c(2, 2),
               sd_pooled = 1)
  expect_identical(names(res),
                   c("item", "factor_loading_reference", "factor_loading_focal",
                     "intercept_reference", "intercept_focal", "sd_pooled",
                     "dmacs"))
  expect_identical(names(generics::tidy(res)),
                   c("term", "estimate", "factor_loading_reference",
                     "factor_loading_focal", "intercept_reference",
                     "intercept_focal", "sd_pooled"))
})

test_that("irt_information() takes discrimination, threshold, and trait_value and reports under them", {
  .rn_has(irt_information,
          c("discrimination", "threshold", "item", "trait_value", "grm"),
          c("a", "b", "theta"))
  .rn_unused(irt_information, list(a = 1.2, b = 0))
  .rn_unused(irt_information, list(discrimination = 1.2, threshold = 0,
                                   theta = c(-1, 0, 1)))
  res <- irt_information(discrimination = c(x1 = 1.2, x2 = .8),
                         threshold = c(-1, .5, 0), item = c("x1", "x1", "x2"),
                         trait_value = c(-1, 0, 1))
  expect_identical(names(res), c("trait_value", "test_information", "se"))
  expect_identical(attr(res, "discrimination"), c(x1 = 1.2, x2 = .8))
  expect_identical(unname(attr(res, "threshold")), c(-1, .5, 0))
  expect_true(attr(res, "trait_value_max_information") %in% c(-1, 0, 1))
  expect_null(attr(res, "a"))
  expect_null(attr(res, "theta_max_information"))
  # Positional calls are unaffected.
  expect_equal(irt_information(1.2, 0)$test_information,
               irt_information(discrimination = 1.2,
                               threshold = 0)$test_information)
})

test_that("simple_structure(), sem_information_criteria(), var_ete(), df_missing(), and df_barnard_rubin() take the English names", {
  L <- rbind(c(.8, .1), c(.7, .0), c(.1, .8), c(.0, .7))
  expect_identical(names(formals(simple_structure))[1], "factor_loading_matrix")
  .rn_unused(simple_structure, list(Lambda = L))
  expect_identical(simple_structure(factor_loading_matrix = L)$value,
                   simple_structure(L)$value)

  .rn_has(sem_information_criteria, "parameter_estimate", "theta")
  .rn_unused(sem_information_criteria,
             list(log_likelihood = -120, npar = 2, N = 200,
                  measures = "SPBIC", theta = c(.4, 1), prior_mean = c(0, 1),
                  information_observed = diag(c(100, 50))))
  ic <- sem_information_criteria(log_likelihood = -120, npar = 2, N = 200,
                                 measures = "SPBIC",
                                 parameter_estimate = c(.4, 1),
                                 prior_mean = c(0, 1),
                                 information_observed = diag(c(100, 50)))
  info <- attr(ic, "sem_information")
  expect_identical(info$parameter_estimate, c(.4, 1))
  expect_null(info$theta)

  .rn_has(var_ete, c("error_variance", "var_Z", "reg_coef_1", "reg_coef_2",
                     "mean_Z"),
          c("sigma2", "sigma2_Z", "beta_1", "beta_2", "mu_Z"))
  .rn_unused(var_ete, list(error_variance = 4, var_Z = 1, n_1 = 20, n_2 = 20,
                           beta_1 = .5, beta_2 = .3))
  .rn_unused(var_ete, list(error_variance = 4, var_Z = 1, n_1 = 20, n_2 = 20,
                           reg_coef_1 = .5, reg_coef_2 = .3, mu_Z = 0))

  .rn_has(df_missing, "df_complete", "nu_com")
  .rn_has(df_barnard_rubin, "df_complete", "nu_com")
  .rn_unused(df_missing, list(fmi = .5, nu_com = 23, method = "complete_df"))
  .rn_unused(df_barnard_rubin, list(between = .04, within = .1, m = 5,
                                    nu_com = 23))
  dm <- df_missing(.5, df_complete = 23, method = "complete_df")
  expect_identical(dm$df_complete, 23)
  expect_false("nu_com" %in% names(dm))
  br <- df_barnard_rubin(between = .04, within = .1, m = 5, df_complete = 23)
  expect_identical(attr(br, "df_complete"), 23)
  expect_null(attr(br, "nu_com"))
})
