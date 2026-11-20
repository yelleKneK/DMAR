## Monte Carlo sensitivity-analysis tests. These are inherently slow; gated by
## skip_on_cran() to keep the CRAN test budget in check. The closed-form
## planners (test-ss_aipe_R2.R, test-ss_aipe_rmsea.R, etc.) cover the same
## functional surface that runs on CRAN.

## ---------------------------------------------------------------------------
## The family-wide return-schema contract.
##
## Every ss_aipe_*_sensitivity() member reports the core rows named in the
## registry vector DMAR:::.SS_AIPE_SENS_CORE_TERMS (the realized-width
## summaries, the width-attainment proportion, and the tail-specific and
## overall empirical non-coverage, all proportions on the 0 to 1 scale),
## plus input echoes named for their unit: a sample-size echo (total_N,
## n_per_group, n_1, n_clusters, or n_units for a crossed design of units
## and raters), a true_* echo of the data generating
## value, width, conf_level, and, when an assurance was supplied, an
## assurance echo (absent otherwise, matching the planner convention). The
## calls below use a tiny G: the contract is about the table's shape and
## scale, not its Monte Carlo precision.
## ---------------------------------------------------------------------------

# One entry per member. None of these calls supplies an assurance, so the
# schema check also asserts the assurance echo row is absent (it appears
# only when one is supplied, matching the planner convention). Members that
# need lavaan are tested separately so the loop does not skip everything
# when lavaan is absent.
.sens_family_calls <- function() {
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  list(
    ss_aipe_R2_sensitivity = quote(ss_aipe_R2_sensitivity(
        true_R2 = 0.30, supposed_R2 = 0.30, width = 0.30, p = 3,
        G = 5, print_iter = FALSE)),
    ss_aipe_c_ancova_sensitivity = quote(ss_aipe_c_ancova_sensitivity(
        true_var_error_ancova = 0.8, supposed_var_error_ancova = 0.8,
        true_var_error_anova = 1, supposed_var_error_anova = 1,
        true_cor = 0.4, supposed_cor = 0.4, G = 5,
        true_mean_Y = c(10, 11, 12), true_sd_Y = 1, true_mean_X = 5, true_sd_X = 1,
        c_weights = c(-1, 0, 1), width = 0.40)),
    ss_aipe_c_sensitivity = quote(ss_aipe_c_sensitivity(
        true_var_error = 1, supposed_var_error = 1,
        c_weights = c(-1, 1), true_contrast = 0.5, width = 1.2,
        G = 5, print_iter = FALSE)),
    ss_aipe_cliff_delta_sensitivity = quote(ss_aipe_cliff_delta_sensitivity(
        true_delta = 0.30, supposed_delta = 0.30, width = 0.50,
        G = 5, print_iter = FALSE)),
    ss_aipe_crossed_sensitivity = quote(ss_aipe_crossed_sensitivity(
        true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
        supposed_var_unit = 0.5, supposed_var_rater = 0.1,
        supposed_var_residual = 0.4, width = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_cv_sensitivity = quote(ss_aipe_cv_sensitivity(
        true_cv = 0.25, supposed_cv = 0.25, width = 0.20,
        G = 5, print_iter = FALSE)),
    ss_aipe_icc_sensitivity = quote(ss_aipe_icc_sensitivity(
        true_icc = 0.8, supposed_icc = 0.8, n_raters = 4, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_indirect_effect_sensitivity = quote(ss_aipe_indirect_effect_sensitivity(
        true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.4, supposed_std_reg_coef_a = 0.4, supposed_std_reg_coef_b = 0.4,
        width = 0.40, G = 5, print_iter = FALSE)),
    ss_aipe_mean_sensitivity = quote(ss_aipe_mean_sensitivity(
        true_sd = 1, supposed_sd = 1, width = 0.50,
        G = 5, print_iter = FALSE)),
    ss_aipe_mixed_effects_sensitivity = quote(ss_aipe_mixed_effects_sensitivity(
        true_total_var = 1, true_var_X = 1, true_icc = 0.2,
        supposed_total_var = 1, supposed_var_X = 1,
        supposed_icc = 0.2, cluster_size = 10, width = 0.20,
        G = 3, print_iter = FALSE)),
    ss_aipe_omega_squared_sensitivity = quote(ss_aipe_omega_squared_sensitivity(
        true_omega_squared = 0.2, supposed_omega_squared = 0.2,
        df_effect = 2, width = 0.30, G = 5, print_iter = FALSE)),
    ss_aipe_partial_r_sensitivity = quote(ss_aipe_partial_r_sensitivity(
        true_cor = 0.4, supposed_cor = 0.4, p = 2, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_pcm_sensitivity = quote(ss_aipe_pcm_sensitivity(
        true_var_trend = 0.1, true_var_error = 1,
        supposed_var_trend = 0.1, supposed_var_error = 1,
        duration = 4, frequency = 1, width = 0.50,
        G = 3, print_iter = FALSE)),
    ss_aipe_r_sensitivity = quote(ss_aipe_r_sensitivity(
        true_cor = 0.4, supposed_cor = 0.4, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_rc_sensitivity = bquote(ss_aipe_rc_sensitivity(
        true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
        supposed_var_Y = 1, supposed_cov_YX = .(cov_YX),
        supposed_cov_XX = .(Sigma_X), which_predictor = 1, width = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_reg_coef_sensitivity = bquote(ss_aipe_reg_coef_sensitivity(
        true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
        supposed_var_Y = 1, supposed_cov_YX = .(cov_YX),
        supposed_cov_XX = .(Sigma_X), which_predictor = 1, width = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_reliability_sensitivity = quote(ss_aipe_reliability_sensitivity(
        true_reliability = 0.8, supposed_reliability = 0.8, n_items = 5,
        width = 0.20, G = 3, print_iter = FALSE)),
    ss_aipe_sc_ancova_sensitivity = quote(ss_aipe_sc_ancova_sensitivity(
        true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
        width = 0.80, true_cor = 0.4, G = 5, print_iter = FALSE)),
    ss_aipe_sc_sensitivity = quote(ss_aipe_sc_sensitivity(
        true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
        width = 0.80, G = 5, print_iter = FALSE)),
    ss_aipe_semipartial_r_sensitivity = quote(ss_aipe_semipartial_r_sensitivity(
        true_r_sp = 0.4, supposed_r_sp = 0.4, p = 2, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_sm_sensitivity = quote(ss_aipe_sm_sensitivity(
        true_sm = 5, supposed_sm = 5, width = 1,
        G = 5, print_iter = FALSE)),
    ss_aipe_smd_sensitivity = quote(ss_aipe_smd_sensitivity(
        true_smd = 0.5, supposed_smd = 0.5, width = 0.80,
        G = 5, print_iter = FALSE)),
    ss_aipe_src_sensitivity = bquote(ss_aipe_src_sensitivity(
        true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
        supposed_var_Y = 1, supposed_cov_YX = .(cov_YX),
        supposed_cov_XX = .(Sigma_X), which_predictor = 1, width = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_equivalence_smd_sensitivity = quote(ss_aipe_equivalence_smd_sensitivity(
        true_smd = 0.1, supposed_smd = 0.1,
        smd_upper = 0.5, width = 0.40,
        G = 3, print_iter = FALSE)),
    ss_aipe_equivalence_r_sensitivity = quote(ss_aipe_equivalence_r_sensitivity(
        true_cor = 0.1, supposed_cor = 0.1,
        cor_upper = 0.5, width = 0.40,
        G = 3, print_iter = FALSE)),
    ss_aipe_two_phase_sensitivity = quote(ss_aipe_two_phase_sensitivity(
        true_sd_Y = 1, true_cor = 0.6, supposed_sd_Y = 1,
        supposed_cor = 0.6, n_pool = 2000, width = 0.40,
        G = 5, print_iter = FALSE))
  )
}

# The shared assertions, used for the closed-form members below and for the
# lavaan-backed members in their own test.
.expect_sens_schema <- function(res, member) {
  core <- DMAR:::.SS_AIPE_SENS_CORE_TERMS
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_false(any(duplicated(res$term)),
               info = paste(member, "has duplicated terms"))
  missing_core <- setdiff(core, res$term)
  if (length(missing_core)) {
    fail(paste(member, "is missing core rows:",
               paste(missing_core, collapse = ", ")))
    return(invisible(NULL))
  }

  v <- stats::setNames(res$value, res$term)

  # Scale contract: the attainment and non-coverage rows are proportions on
  # the 0 to 1 scale, and the total Type I error is the sum of its tails.
  props <- v[c("pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high",
               "total_type_I_error")]
  expect_false(anyNA(props), info = member)
  expect_true(all(props >= 0 & props <= 1), info = member)
  expect_equal(unname(v["total_type_I_error"]),
               unname(v["pct_ci_miss_low"] + v["pct_ci_miss_high"]),
               tolerance = 1e-12, info = member)

  # Width summaries are nonnegative and ordered sensibly.
  expect_true(v["mean_ci_width"] > 0, info = member)
  expect_true(v["sd_ci_width"] >= 0, info = member)

  # Echoes named for their unit: a sample-size echo, a true_* echo of the
  # data generating value, the width target, and the confidence level.
  expect_true(any(c("total_N", "n_per_group", "n_1", "n_clusters", "n_units")
                  %in% res$term), info = member)
  expect_true(any(startsWith(res$term, "true_")), info = member)
  expect_true("width" %in% res$term, info = member)
  expect_true("conf_level" %in% res$term, info = member)
  expect_false(is.na(v["width"]), info = member)
  expect_false(is.na(v["conf_level"]), info = member)

  # Members that accept an assurance argument echo it only when one was
  # supplied, matching the planner convention (see test-ss_aipe_smd.R).
  # None of the calls in this file's family loop supply one.
  expect_false("assurance" %in% res$term, info = member)

  .expect_sens_record(res, member)
}

# The record of the run every sensitivity function attaches as
# attr(, "simulation") (DMAR:::.ss_sens_simulation()): the three counts,
# the one denominator, and a width entry that is the table's own
# pct_ci_less_width row. `G` is the number of replications the call asked for;
# `width_target` is FALSE for a function with no width target.
.expect_sens_record <- function(res, member, G = NULL, width_target = TRUE) {
  s <- attr(res, "simulation")
  expect_type(s, "list")
  expect_true(all(c("attempted", "returned", "failed", "denominator",
                    "width_target_probability", "width_target_mc_ci")
                  %in% names(s)), info = member)
  if (!is.null(G)) expect_identical(s$attempted, as.integer(G), info = member)
  expect_type(s$attempted, "integer")
  expect_identical(s$returned + s$failed, s$attempted, info = member)
  expect_identical(s$denominator, s$returned, info = member)
  if (width_target) {
    pct <- res$value[res$term == "pct_ci_less_width"]
    expect_equal(s$width_target_probability, pct, tolerance = 1e-12,
                 info = member)
    expect_equal(as.numeric(s$width_target_mc_ci),
                 as.numeric(stats::binom.test(
                   round(pct * s$returned), s$returned,
                   conf.level = attr(res, "conf_level"))$conf.int),
                 tolerance = 1e-12, info = member)
  } else {
    expect_identical(s$width_target_probability, NA_real_, info = member)
    expect_identical(s$width_target_mc_ci, c(NA_real_, NA_real_),
                     info = member)
  }
  expect_identical(dmar_metadata(res)$simulation, s, info = member)
  expect_match(DMAR:::.dmar_tbl_footer(res)[["simulation"]],
               paste0("^Replications: ", s$attempted), info = member)
}

test_that("every closed-form ss_aipe_*_sensitivity() member honors the family return schema", {
  skip_on_cran()
  members <- .sens_family_calls()
  for (member in names(members)) {
    set.seed(113)
    res <- suppressWarnings(suppressMessages(eval(members[[member]])))
    .expect_sens_schema(res, member)
    expect_identical(attr(res, "simulation")$attempted,
                     as.integer(members[[member]]$G), info = member)
  }
})

test_that("every exported sensitivity function is in a test of the record", {
  exported <- grep("_sensitivity$", getNamespaceExports("DMAR"), value = TRUE)
  # The members of the loop above, the lavaan-backed members, and the
  # functions with no width target, each tested in this file.
  tested <- c(names(.sens_family_calls()),
              "ss_aipe_rmsea_sensitivity", "ss_aipe_sem_path_sensitivity",
              "ss_power_R2_sensitivity", "ss_power_reg_coef_sensitivity",
              "ss_seq_c_sensitivity")
  expect_setequal(exported, tested)
})

test_that("the sensitivity functions with no width target record the counts alone", {
  skip_on_cran()
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  res <- suppressMessages(ss_power_R2_sensitivity(
    true_R2 = 0.20, specified_N = 60, p = 3, G = 6, print_iter = FALSE,
    seed = 113))
  .expect_sens_record(res, "ss_power_R2_sensitivity", G = 6,
                      width_target = FALSE)
  res <- suppressMessages(ss_power_reg_coef_sensitivity(
    true_var_Y = 1, true_cov_YX = c(0.4, 0.3), true_cov_XX = Sigma_X,
    which_predictor = 1, specified_N = 60, G = 6, print_iter = FALSE,
    seed = 113))
  .expect_sens_record(res, "ss_power_reg_coef_sensitivity", G = 6,
                      width_target = FALSE)
  res <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                              true_sd = 15.67, G = 10, seed = 113)
  .expect_sens_record(res, "ss_seq_c_sensitivity", G = 10,
                      width_target = FALSE)
  expect_identical(attr(res, "simulation")$failed, 0L)
})

test_that("a power sensitivity table with a replication without a test says so in its footer", {
  x <- data.frame(term = "empirical_power", value = .8)
  class(x) <- c("dmar_ss_power_sensitivity", "dmar_tbl", "data.frame")
  x <- DMAR:::.ss_sens_record(x, 10, rep(c(TRUE, FALSE), c(9L, 1L)))
  expect_identical(DMAR:::.dmar_tbl_footer(x)[["simulation"]],
                   "Replications: 10 attempted, 9 with a test, 1 without")
})

test_that("the lavaan-backed sensitivity members honor the family return schema", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  # ss_aipe_rmsea_sensitivity: the model is deliberately misspecified so the
  # population RMSEA is positive.
  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
  dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  set.seed(113)
  res <- suppressWarnings(suppressMessages(
    ss_aipe_rmsea_sensitivity(width = 0.05,
                              model = "g =~ x1 + x2 + x3 + x4 + x5 + x6",
                              true_cov_matrix = Sigma, specified_N = 200, G = 5)
  ))
  .expect_sens_schema(res, "ss_aipe_rmsea_sensitivity")
  expect_identical(attr(res, "simulation")$attempted, 5L)
  expect_lte(attr(res, "simulation")$returned, 5L)
  expect_true(all(c("df", "true_rmsea") %in% res$term))
  expect_false("suc_rep" %in% res$term)

  # ss_aipe_sem_path_sensitivity: a two-factor model with a labeled
  # structural path.
  pop_model <- "
    f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
    f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
    f2 ~ 0.5*f1
    f1 ~~ 1*f1
    f2 ~~ 0.75*f2
    y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
    y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
  "
  Sigma_sem <- cov_sem(pop_model)$sigma_theta
  analysis_model <- "
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f2 ~ b*f1
  "
  set.seed(113)
  out <- suppressWarnings(suppressMessages(
    ss_aipe_sem_path_sensitivity(model = analysis_model,
                                 supposed_cov_matrix = Sigma_sem,
                                 true_cov_matrix = Sigma_sem,
                                 which_path = "b", width = 0.30,
                                 specified_N = 150, G = 5)
  ))
  .expect_sens_schema(out, "ss_aipe_sem_path_sensitivity")
  expect_identical(attr(out, "simulation")$attempted, 5L)
  expect_lte(attr(out, "simulation")$returned, 5L)
  # Data sets whose fit did not converge are drawn again and counted.
  expect_type(attr(out, "simulation")$n_redrawn, "integer")
  expect_gte(attr(out, "simulation")$n_redrawn, 0L)
  expect_true("true_path" %in% out$term)
  expect_false("suc_rep" %in% out$term)

  # ss_aipe_reliability_sensitivity under missing data (FIML in every
  # replication) and with the omega estimator: the same schema, with the
  # missing data inputs carried as attributes rather than rows.
  for (variant in list(list(missing_prop = 0.2, fmi = 0.18),
                       list(estimator = "omega"))) {
    set.seed(113)
    rel <- suppressWarnings(suppressMessages(do.call(
      ss_aipe_reliability_sensitivity,
      c(list(true_reliability = 0.8, supposed_reliability = 0.8, n_items = 5,
             width = 0.20, G = 3), variant))))
    .expect_sens_schema(rel, paste0("ss_aipe_reliability_sensitivity/",
                                    names(variant)[1]))
  }
})

test_that("ss_seq_c_sensitivity() echoes its inputs under their unit names", {
  skip_on_cran()
  res <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                              true_sd = 15.67, true_means = c(70, 65),
                              G = 10, seed = 113)
  expect_named(res, c("term", "value"))
  expect_true(all(c("n_star", "mean_N", "median_N", "sd_N",
                    "ratio_mean_N_n_star", "coverage", "se_coverage",
                    "half_width", "true_contrast", "true_sd",
                    "alpha_level", "m0") %in% res$term))
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v["half_width"]), 2.5)
  expect_equal(unname(v["true_contrast"]), 5)
  expect_equal(unname(v["true_sd"]), 15.67)
  expect_equal(unname(v["alpha_level"]), 0.05)
  expect_equal(unname(v["m0"]), 10)
})

## ---------------------------------------------------------------------------
## Member-specific behavior beyond the shared schema.
## ---------------------------------------------------------------------------

test_that("ss_aipe_R2_sensitivity() reports the R2 summaries, limit summaries, and echoes", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_R2_sensitivity(true_R2 = 0.30, supposed_R2 = 0.30, width = 0.20,
                           p = 5, G = 20, print_iter = FALSE)
  )
  expect_true(all(c("mean_R2", "median_R2", "sd_R2",
                    "mean_ci_lower", "mean_ci_upper",
                    "mean_ci_width_lower", "mean_ci_width_upper",
                    "total_N", "p", "true_R2", "supposed_R2") %in% res$term))
  # The count of replications without an interval is in the simulation
  # record, not in a row.
  expect_false("num_probs_with_cis" %in% res$term)
  expect_identical(attr(res, "simulation")$attempted, 20L)
  expect_identical(attr(res, "simulation")$failed,
                   20L - attr(res, "simulation")$returned)
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v["true_R2"]), 0.30)
  expect_equal(unname(v["supposed_R2"]), 0.30)
  expect_equal(unname(v["p"]), 5)
  expect_equal(unname(v["width"]), 0.20)
  # The one-sided widths add up to the full width in the mean.
  expect_equal(unname(v["mean_ci_width_lower"] + v["mean_ci_width_upper"]),
               unname(v["mean_ci_width"]), tolerance = 1e-12)
})

test_that("ss_aipe_smd_sensitivity() reports the realized SMD beside the widths", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_smd_sensitivity(true_smd = 0.5, supposed_smd = 0.5,
                            width = 0.20, G = 20, print_iter = FALSE)
  )
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_smd", "median_smd", "sd_smd") %in% res$term))
  expect_lt(abs(v["mean_smd"] - 0.5), 0.15)
  expect_equal(unname(v["total_N"]), unname(2 * v["n_per_group"]))
  expect_equal(unname(v["true_smd"]), 0.5)
})

test_that("ss_aipe_sm_sensitivity() reports the realized standardized mean beside the widths", {
  skip_on_cran()
  set.seed(113)
  # At true_sm = 5 the noncentrality parameter passes 37.62, where
  # ci_nc_t() warned in DMAR 1.0.x; its limits are now exact there, and
  # the call is silent.
  res <- expect_no_warning(suppressMessages(
    ss_aipe_sm_sensitivity(true_sm = 5, supposed_sm = 5,
                           width = 0.5, G = 20, print_iter = FALSE)
  ))
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_sm", "median_sm", "sd_sm") %in% res$term))
  expect_lt(abs(v["mean_sm"] - 5), 1)
  expect_equal(unname(v["true_sm"]), 5)
})

test_that("ss_aipe_mean_sensitivity() reports the realized mean beside the widths, for either population", {
  skip_on_cran()
  res <- ss_aipe_mean_sensitivity(true_sd = 2, supposed_sd = 2, width = 1,
                                  true_mean = 10, G = 20, seed = 113)
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_mean", "median_mean", "sd_mean", "total_N",
                    "true_sd", "true_mean", "supposed_sd") %in% res$term))
  expect_lt(abs(v["mean_mean"] - 10), 1)
  expect_equal(unname(v["true_mean"]), 10)
  # A resampled population and a fixed size: the same schema, with the
  # planning value echoed as NA.
  res_p <- ss_aipe_mean_sensitivity(specified_N = 30, width = 8,
                                    population = pygmalion$iq_gain,
                                    G = 5, seed = 113)
  .expect_sens_schema(res_p, "ss_aipe_mean_sensitivity/population")
  expect_true(is.na(res_p$value[res_p$term == "supposed_sd"]))
})

test_that("ss_aipe_sc_sensitivity() reports the realized contrast, both one-sided widths, and echoes", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_sc_sensitivity(true_std_contrast = 0.5, supposed_std_contrast = 0.5,
                           c_weights = c(-1, 0, 1),
                           width = 0.40, conf_level = 0.95,
                           G = 50, print_iter = FALSE)
  )
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_std_contrast", "median_std_contrast", "sd_std_contrast",
                    "mean_ci_width_lower", "mean_ci_width_upper") %in% res$term))
  expect_lt(abs(v["mean_std_contrast"] - 0.5), 0.15)
  expect_equal(unname(v["total_N"]), unname(3 * v["n_per_group"]))
  # A well-specified plan attains the target width often; on the percentage
  # scale this row would exceed 1, so the bound is a real scale check.
  expect_gt(v["pct_ci_less_width"], 0)
  expect_true(all(v[c("pct_ci_miss_low", "pct_ci_miss_high")] >= 0 &
                    v[c("pct_ci_miss_low", "pct_ci_miss_high")] <= 1))
})

test_that("ss_aipe_sc_ancova_sensitivity() honors the schema on both divisor branches", {
  skip_on_cran()
  for (divisor in c("s_ancova", "s_anova")) {
    set.seed(113)
    res <- suppressMessages(
      ss_aipe_sc_ancova_sensitivity(
        true_std_contrast = 0.5, supposed_std_contrast = 0.5,
        c_weights = c(-1, 0, 1),
        width = 0.40, true_cor = 0.4, divisor = divisor,
        conf_level = 0.95, G = 20, print_iter = FALSE
      )
    )
    .expect_sens_schema(res, paste0("ss_aipe_sc_ancova_sensitivity/", divisor))
    v <- stats::setNames(res$value, res$term)
    expect_equal(unname(v["true_cor"]), 0.4)
    expect_equal(unname(v["true_std_contrast"]), 0.5)
  }
})

test_that("ss_aipe_c_ancova_sensitivity() reports proportions and the resolved planning inputs", {
  skip_on_cran()
  # true_sd_Y must equal sqrt(true_var_error_anova) per the function's own
  # input-consistency check. Use true_sd_Y = 1 with anova error var = 1.
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_c_ancova_sensitivity(
      true_var_error_ancova = 0.8, supposed_var_error_ancova = 0.8,
      true_var_error_anova = 1, supposed_var_error_anova = 1,
      true_cor = 0.4, supposed_cor = 0.4, G = 50,
      true_mean_Y = c(10, 11, 12), true_sd_Y = 1, true_mean_X = 5, true_sd_X = 1,
      c_weights = c(-1, 0, 1), width = 0.40, conf_level = 0.95
    )
  )
  v <- stats::setNames(res$value, res$term)
  expect_true("mean_se_ratio" %in% res$term)
  expect_equal(unname(v["supposed_var_error_ancova"]), 0.8)
  expect_equal(unname(v["true_cor"]), 0.4)
  # The population contrast implied by true_mean_Y and the weights.
  expect_equal(unname(v["true_contrast"]), sum(c(-1, 0, 1) * c(10, 11, 12)))
  expect_equal(unname(v["total_N"]), unname(3 * v["n_per_group"]))
  # Proportions on the 0 to 1 scale (this member once reported percentages).
  expect_true(all(v[c("pct_ci_less_width", "pct_ci_miss_low",
                      "pct_ci_miss_high", "total_type_I_error")] <= 1))
})

test_that("ss_aipe_reg_coef_sensitivity() echoes the implied coefficients", {
  skip_on_cran()
  set.seed(113)
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  res <- suppressMessages(
    ss_aipe_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      supposed_var_Y = 1, supposed_cov_YX = cov_YX, supposed_cov_XX = Sigma_X,
      which_predictor = 1, width = 0.40, conf_level = 0.95,
      G = 50, print_iter = FALSE
    )
  )
  v <- stats::setNames(res$value, res$term)
  b <- (solve(Sigma_X) %*% cov_YX)[1]
  expect_equal(unname(v["true_reg_coef_j"]), b, tolerance = 1e-10)
  expect_equal(unname(v["supposed_reg_coef_j"]), b, tolerance = 1e-10)
  expect_equal(unname(v["p"]), 2)
  expect_equal(unname(v["which_predictor"]), 1)
  expect_true(all(c("mean_R2", "median_R2", "sd_R2") %in% res$term))
})

test_that("ss_aipe_rc_sensitivity() and ss_aipe_src_sensitivity() share the engine of ss_aipe_reg_coef_sensitivity()", {
  skip_on_cran()
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  for (fn in list(ss_aipe_rc_sensitivity, ss_aipe_src_sensitivity)) {
    set.seed(113)
    res <- suppressMessages(
      fn(true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
         supposed_var_Y = 1, supposed_cov_YX = cov_YX,
         supposed_cov_XX = Sigma_X,
         which_predictor = 1, width = 0.40, conf_level = 0.95,
         G = 50, print_iter = FALSE)
    )
    expect_s3_class(res, "data.frame")
    expect_true(all(DMAR:::.SS_AIPE_SENS_CORE_TERMS %in% res$term))
  }
})

test_that("ss_aipe_src_sensitivity() names its coefficient rows std_reg_coef_j, the standardized coefficient", {
  skip_on_cran()
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  run <- function(fn) {
    set.seed(113)
    suppressMessages(
      fn(true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
         supposed_var_Y = 1, supposed_cov_YX = cov_YX,
         supposed_cov_XX = Sigma_X,
         which_predictor = 1, width = 0.40, conf_level = 0.95,
         G = 50, print_iter = FALSE))
  }
  src <- run(ss_aipe_src_sensitivity)
  rc  <- run(ss_aipe_rc_sensitivity)
  beta_rows <- c("mean_std_reg_coef_j", "median_std_reg_coef_j",
                 "sd_std_reg_coef_j", "true_std_reg_coef_j",
                 "supposed_std_reg_coef_j")
  b_rows <- sub("std_reg_coef_j", "reg_coef_j", beta_rows)
  expect_true(all(beta_rows %in% src$term))
  expect_false(any(b_rows %in% src$term))
  expect_true(all(b_rows %in% rc$term))
  expect_false(any(beta_rows %in% rc$term))
  expect_identical(attr(src, "simulation"), attr(src, "simulation"))
  expect_false(is.null(attr(src, "simulation")))
})

test_that("a supplied assurance value is echoed under its own name", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_smd_sensitivity(true_smd = 0.5, supposed_smd = 0.5,
                            width = 0.30, assurance = 0.90,
                            G = 5, print_iter = FALSE)
  )
  expect_equal(res$value[res$term == "assurance"], 0.90)
})

test_that("evaluating at a fixed size echoes NA for the planning value", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_smd_sensitivity(true_smd = 0.5, specified_n_per_group = 30,
                            width = 0.30,
                            G = 5, print_iter = FALSE)
  )
  v <- stats::setNames(res$value, res$term)
  expect_true(is.na(v["supposed_smd"]))
  expect_equal(unname(v["n_per_group"]), 30)
  expect_equal(unname(v["total_N"]), 60)
})

test_that("the contrast, SMD, and CV members write per-replication results only to a named file", {
  skip_on_cran()
  # Nothing is written unless the caller names a file; the old `save`
  # switch and its default path in the working directory are gone, and a
  # malformed filename is refused before any replication runs. Each entry
  # is a function of the filename so the same call runs three times: once
  # refused, once creating the file, once appending to it.
  members <- list(
    ss_aipe_c_ancova_sensitivity = function(f) ss_aipe_c_ancova_sensitivity(
      true_var_error_ancova = 0.8, supposed_var_error_ancova = 0.8, true_cor = 0.4,
      G = 3, true_mean_Y = c(10, 11, 12), true_sd_Y = 1, true_mean_X = 5, true_sd_X = 1,
      c_weights = c(-1, 0, 1), width = 0.40, filename = f),
    ss_aipe_sc_sensitivity = function(f) ss_aipe_sc_sensitivity(
      true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
      width = 0.80, G = 3, print_iter = FALSE, filename = f),
    ss_aipe_smd_sensitivity = function(f) ss_aipe_smd_sensitivity(
      true_smd = 0.5, supposed_smd = 0.5, width = 0.80,
      G = 3, print_iter = FALSE, filename = f),
    ss_aipe_cv_sensitivity = function(f) ss_aipe_cv_sensitivity(
      true_cv = 0.25, supposed_cv = 0.25, width = 0.20, G = 3,
      print_iter = FALSE, filename = f)
  )
  columns <- list(
    ss_aipe_c_ancova_sensitivity = c("contrast", "se_contrast", "se_contrast_restricted",
                                     "ratio", "width_obs", "type_I_error",
                                     "type_I_error_upper", "type_I_error_lower"),
    ss_aipe_sc_sensitivity = c("std_contrast", "full_width", "width_lower",
                               "width_upper", "type_I_error_upper",
                               "type_I_error_lower", "type_I_error",
                               "ci_lower", "ci_upper"),
    ss_aipe_smd_sensitivity = c("d", "full_width", "width_lower",
                                "width_upper", "type_I_error_upper",
                                "type_I_error_lower", "type_I_error",
                                "ci_lower", "ci_upper"),
    ss_aipe_cv_sensitivity = c("ci_lower", "ci_upper", "cv", "int_ok",
                               "Width")
  )
  for (member in names(members)) {
    fn <- get(member)
    expect_null(eval(formals(fn)$filename), info = member)
    expect_false("save" %in% names(formals(fn)), info = member)
    expect_error(members[[member]](42),
                 "'filename' must be NULL or a single character string",
                 fixed = TRUE, info = member)

    out_file <- tempfile(fileext = ".csv")
    set.seed(113)
    suppressWarnings(members[[member]](out_file))
    written <- utils::read.csv(out_file)
    expect_named(written, columns[[member]], info = member)
    expect_equal(nrow(written), 3L, info = member)

    # A second run appends to the existing file and says so.
    set.seed(113)
    expect_message(suppressWarnings(members[[member]](out_file)),
                   "already exists", fixed = TRUE, info = member)
    expect_equal(nrow(utils::read.csv(out_file)), 6L, info = member)
    unlink(out_file)
  }
})

test_that("the correlation, reliability, and regression coefficient members write per-replication results only to a named file", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  # Nothing is written unless the caller names a file; the old `save`
  # switch and its default path in the working directory are gone, and a
  # malformed filename is refused before any replication runs. The two
  # regression coefficient functions forward `filename` alone to the engine
  # they share with ss_aipe_reg_coef_sensitivity(), which does their writing, so for them
  # the check is that the forwarded path receives one row per replication
  # and is appended to on a second run; the column names are that
  # function's to test.
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  members <- list(
    ss_aipe_partial_r_sensitivity = function(f) ss_aipe_partial_r_sensitivity(
      true_cor = 0.4, supposed_cor = 0.4, p = 3, width = 0.40,
      G = 3, print_iter = FALSE, filename = f),
    ss_aipe_semipartial_r_sensitivity = function(f) ss_aipe_semipartial_r_sensitivity(
      true_r_sp = 0.4, supposed_r_sp = 0.4, p = 3, width = 0.40,
      G = 3, print_iter = FALSE, filename = f),
    ss_aipe_reliability_sensitivity = function(f) ss_aipe_reliability_sensitivity(
      true_reliability = 0.8, supposed_reliability = 0.8, n_items = 5,
      width = 0.20, G = 3, print_iter = FALSE, filename = f),
    ss_aipe_rc_sensitivity = function(f) ss_aipe_rc_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      supposed_var_Y = 1, supposed_cov_YX = cov_YX,
      supposed_cov_XX = Sigma_X, which_predictor = 1, width = 0.60,
      G = 3, print_iter = FALSE, filename = f),
    ss_aipe_src_sensitivity = function(f) ss_aipe_src_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      supposed_var_Y = 1, supposed_cov_YX = cov_YX,
      supposed_cov_XX = Sigma_X, which_predictor = 1, width = 0.60,
      G = 3, print_iter = FALSE, filename = f)
  )
  columns <- list(
    ss_aipe_partial_r_sensitivity = c("partial_r", "ci_lower", "ci_upper",
                                      "ci_width", "type_I_lower",
                                      "type_I_upper"),
    ss_aipe_semipartial_r_sensitivity = c("r_sp", "ci_lower", "ci_upper",
                                          "ci_width", "type_I_lower",
                                          "type_I_upper"),
    ss_aipe_reliability_sensitivity = c("rel_hat", "ci_lower", "ci_upper",
                                        "ci_width", "type_I_lower",
                                        "type_I_upper")
  )
  for (member in names(members)) {
    fn <- get(member)
    expect_null(eval(formals(fn)$filename), info = member)
    expect_false("save" %in% names(formals(fn)), info = member)
    expect_error(members[[member]](42), "filename", info = member)

    out_file <- tempfile(fileext = ".csv")
    set.seed(113)
    invisible(utils::capture.output(
      suppressMessages(suppressWarnings(members[[member]](out_file)))
    ))
    written <- utils::read.csv(out_file)
    if (!is.null(columns[[member]]))
      expect_named(written, columns[[member]], info = member)
    expect_equal(nrow(written), 3L, info = member)

    # A second run appends to the existing file without repeating the header.
    set.seed(113)
    invisible(utils::capture.output(
      suppressMessages(suppressWarnings(members[[member]](out_file)))
    ))
    expect_equal(nrow(utils::read.csv(out_file)), 6L, info = member)
    unlink(out_file)
  }
})
