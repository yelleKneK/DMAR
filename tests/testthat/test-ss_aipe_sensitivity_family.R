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
## n_per_group, n_1, or n_clusters), a true_* echo of the data generating
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
        true_R2 = 0.30, estimated_R2 = 0.30, w = 0.30, p = 3,
        G = 5, print_iter = FALSE)),
    ss_aipe_c_ancova_sensitivity = quote(ss_aipe_c_ancova_sensitivity(
        true_error_var_ancova = 0.8, est_error_var_ancova = 0.8,
        true_error_var_anova = 1, est_error_var_anova = 1,
        rho = 0.4, est_rho = 0.4, G = 5,
        mu_y = c(10, 11, 12), sigma_y = 1, mu_x = 5, sigma_x = 1,
        c_weights = c(-1, 0, 1), width = 0.40)),
    ss_aipe_c_sensitivity = quote(ss_aipe_c_sensitivity(
        true_error_variance = 1, estimated_error_variance = 1,
        c_weights = c(-1, 1), true_psi = 0.5, width = 1.2,
        G = 5, print_iter = FALSE)),
    ss_aipe_cliff_delta_sensitivity = quote(ss_aipe_cliff_delta_sensitivity(
        true_delta = 0.30, estimated_delta = 0.30, width = 0.50,
        G = 5, print_iter = FALSE)),
    ss_aipe_cv_sensitivity = quote(ss_aipe_cv_sensitivity(
        true_cv = 0.25, estimated_cv = 0.25, width = 0.20,
        G = 5, print_iter = FALSE)),
    ss_aipe_icc_sensitivity = quote(ss_aipe_icc_sensitivity(
        true_rho = 0.8, estimated_rho = 0.8, k = 4, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_indirect_effect_sensitivity = quote(ss_aipe_indirect_effect_sensitivity(
        true_a = 0.4, true_b = 0.4, estimated_a = 0.4, estimated_b = 0.4,
        width = 0.40, G = 5, print_iter = FALSE)),
    ss_aipe_mixed_effects_sensitivity = quote(ss_aipe_mixed_effects_sensitivity(
        true_sigma2_y = 1, true_sigma2_x = 1, true_icc = 0.2,
        estimated_sigma2_y = 1, estimated_sigma2_x = 1,
        estimated_icc = 0.2, cluster_size = 10, width = 0.20,
        G = 3, print_iter = FALSE)),
    ss_aipe_omega_squared_sensitivity = quote(ss_aipe_omega_squared_sensitivity(
        true_omega_squared = 0.2, estimated_omega_squared = 0.2,
        df_effect = 2, width = 0.30, G = 5, print_iter = FALSE)),
    ss_aipe_partial_r_sensitivity = quote(ss_aipe_partial_r_sensitivity(
        true_rho = 0.4, estimated_rho = 0.4, J = 1, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_pcm_sensitivity = quote(ss_aipe_pcm_sensitivity(
        true_variance_trend = 0.1, true_error_variance = 1,
        estimated_variance_trend = 0.1, estimated_error_variance = 1,
        duration = 4, frequency = 1, width = 0.50,
        G = 3, print_iter = FALSE)),
    ss_aipe_r_sensitivity = quote(ss_aipe_r_sensitivity(
        true_rho = 0.4, estimated_rho = 0.4, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_rc_sensitivity = bquote(ss_aipe_rc_sensitivity(
        true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
        estimated_var_Y = 1, estimated_cov_YX = .(cov_YX),
        estimated_cov_XX = .(Sigma_X), which_predictor = 1, w = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_reg_coef_sensitivity = bquote(ss_aipe_reg_coef_sensitivity(
        true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
        estimated_var_Y = 1, estimated_cov_YX = .(cov_YX),
        estimated_cov_XX = .(Sigma_X), which_predictor = 1, w = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_reliability_sensitivity = quote(ss_aipe_reliability_sensitivity(
        true_reliability = 0.8, estimated_reliability = 0.8, i = 5,
        width = 0.20, G = 3, print_iter = FALSE)),
    ss_aipe_sc_ancova_sensitivity = quote(ss_aipe_sc_ancova_sensitivity(
        true_psi = 0.5, estimated_psi = 0.5, c_weights = c(-1, 0, 1),
        desired_width = 0.80, rho = 0.4, G = 5, print_iter = FALSE)),
    ss_aipe_sc_sensitivity = quote(ss_aipe_sc_sensitivity(
        true_psi = 0.5, estimated_psi = 0.5, c_weights = c(-1, 0, 1),
        desired_width = 0.80, G = 5, print_iter = FALSE)),
    ss_aipe_semipartial_r_sensitivity = quote(ss_aipe_semipartial_r_sensitivity(
        true_r_sp = 0.4, estimated_r_sp = 0.4, J = 1, width = 0.40,
        G = 5, print_iter = FALSE)),
    ss_aipe_sm_sensitivity = quote(ss_aipe_sm_sensitivity(
        true_sm = 5, estimated_sm = 5, desired_width = 1,
        G = 5, print_iter = FALSE)),
    ss_aipe_smd_sensitivity = quote(ss_aipe_smd_sensitivity(
        true_delta = 0.5, estimated_delta = 0.5, desired_width = 0.80,
        G = 5, print_iter = FALSE)),
    ss_aipe_src_sensitivity = bquote(ss_aipe_src_sensitivity(
        true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
        estimated_var_Y = 1, estimated_cov_YX = .(cov_YX),
        estimated_cov_XX = .(Sigma_X), which_predictor = 1, w = 0.60,
        G = 5, print_iter = FALSE)),
    ss_aipe_equivalence_smd_sensitivity = quote(ss_aipe_equivalence_smd_sensitivity(
        true_smd = 0.1, estimated_smd = 0.1,
        delta_upper = 0.5, width = 0.40,
        G = 3, print_iter = FALSE)),
    ss_aipe_equivalence_r_sensitivity = quote(ss_aipe_equivalence_r_sensitivity(
        true_r = 0.1, estimated_r = 0.1,
        rho_upper = 0.5, width = 0.40,
        G = 3, print_iter = FALSE))
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
  props <- v[c("pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high",
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
  expect_true(any(c("total_N", "n_per_group", "n_1", "n_clusters")
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
}

test_that("every closed-form ss_aipe_*_sensitivity() member honors the family return schema", {
  skip_on_cran()
  members <- .sens_family_calls()
  for (member in names(members)) {
    set.seed(113)
    res <- suppressWarnings(suppressMessages(eval(members[[member]])))
    .expect_sens_schema(res, member)
  }
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
                              Sigma = Sigma, N = 200, G = 5)
  ))
  .expect_sens_schema(res, "ss_aipe_rmsea_sensitivity")
  expect_true(all(c("suc_rep", "df", "true_rmsea") %in% res$term))

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
                                 est_Sigma = Sigma_sem,
                                 true_Sigma = Sigma_sem,
                                 which_path = "b", desired_width = 0.30,
                                 N = 150, G = 5)
  ))
  .expect_sens_schema(out, "ss_aipe_sem_path_sensitivity")
  expect_true(all(c("suc_rep", "true_path") %in% out$term))
})

test_that("ss_seq_c_sensitivity() echoes its inputs under their unit names", {
  skip_on_cran()
  res <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                              true_sigma = 15.67, true_means = c(70, 65),
                              G = 10, seed = 113)
  expect_named(res, c("term", "value"))
  expect_true(all(c("n_star", "mean_N", "median_N", "sd_N",
                    "ratio_mean_N_n_star", "coverage", "se_coverage",
                    "half_width", "true_psi", "true_sigma",
                    "alpha_level", "m0") %in% res$term))
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v["half_width"]), 2.5)
  expect_equal(unname(v["true_psi"]), 5)
  expect_equal(unname(v["true_sigma"]), 15.67)
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
    ss_aipe_R2_sensitivity(true_R2 = 0.30, estimated_R2 = 0.30, w = 0.20,
                           p = 5, G = 20, print_iter = FALSE)
  )
  expect_true(all(c("mean_R2", "median_R2", "sd_R2",
                    "mean_lower_limit", "mean_upper_limit",
                    "mean_ci_width_lower", "mean_ci_width_upper",
                    "num_probs_with_cis",
                    "total_N", "p", "true_R2", "estimated_R2") %in% res$term))
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v["true_R2"]), 0.30)
  expect_equal(unname(v["estimated_R2"]), 0.30)
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
    ss_aipe_smd_sensitivity(true_delta = 0.5, estimated_delta = 0.5,
                            desired_width = 0.20, G = 20, print_iter = FALSE)
  )
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_smd", "median_smd", "sd_smd") %in% res$term))
  expect_lt(abs(v["mean_smd"] - 0.5), 0.15)
  expect_equal(unname(v["total_N"]), unname(2 * v["n_per_group"]))
  expect_equal(unname(v["true_delta"]), 0.5)
})

test_that("ss_aipe_sm_sensitivity() reports the realized standardized mean beside the widths", {
  skip_on_cran()
  set.seed(113)
  # At true_sm = 5 the noncentrality parameter is large enough that
  # conf_limits_nct() warns about the accurate range of R's noncentral t
  # functions; that warning is expected here and not under test.
  res <- suppressWarnings(suppressMessages(
    ss_aipe_sm_sensitivity(true_sm = 5, estimated_sm = 5,
                           desired_width = 0.5, G = 20, print_iter = FALSE)
  ))
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_sm", "median_sm", "sd_sm") %in% res$term))
  expect_lt(abs(v["mean_sm"] - 5), 1)
  expect_equal(unname(v["true_sm"]), 5)
})

test_that("ss_aipe_sc_sensitivity() reports the realized contrast, both one-sided widths, and echoes", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_sc_sensitivity(true_psi = 0.5, estimated_psi = 0.5,
                           c_weights = c(-1, 0, 1),
                           desired_width = 0.40, conf_level = 0.95,
                           G = 50, print_iter = FALSE)
  )
  v <- stats::setNames(res$value, res$term)
  expect_true(all(c("mean_psi", "median_psi", "sd_psi",
                    "mean_ci_width_lower", "mean_ci_width_upper") %in% res$term))
  expect_lt(abs(v["mean_psi"] - 0.5), 0.15)
  expect_equal(unname(v["total_N"]), unname(3 * v["n_per_group"]))
  # A well-specified plan attains the target width often; on the percentage
  # scale this row would exceed 1, so the bound is a real scale check.
  expect_gt(v["pct_ci_less_w"], 0)
  expect_true(all(v[c("pct_ci_miss_low", "pct_ci_miss_high")] >= 0 &
                    v[c("pct_ci_miss_low", "pct_ci_miss_high")] <= 1))
})

test_that("ss_aipe_sc_ancova_sensitivity() honors the schema on both divisor branches", {
  skip_on_cran()
  for (divisor in c("s_ancova", "s_anova")) {
    set.seed(113)
    res <- suppressMessages(
      ss_aipe_sc_ancova_sensitivity(
        true_psi = 0.5, estimated_psi = 0.5,
        c_weights = c(-1, 0, 1),
        desired_width = 0.40, rho = 0.4, divisor = divisor,
        conf_level = 0.95, G = 20, print_iter = FALSE
      )
    )
    .expect_sens_schema(res, paste0("ss_aipe_sc_ancova_sensitivity/", divisor))
    v <- stats::setNames(res$value, res$term)
    expect_equal(unname(v["rho"]), 0.4)
    expect_equal(unname(v["true_psi"]), 0.5)
  }
})

test_that("ss_aipe_c_ancova_sensitivity() reports proportions and the resolved planning inputs", {
  skip_on_cran()
  # sigma_y must equal sqrt(true_error_var_anova) per the function's own
  # input-consistency check. Use sigma_y = 1 with anova error var = 1.
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_c_ancova_sensitivity(
      true_error_var_ancova = 0.8, est_error_var_ancova = 0.8,
      true_error_var_anova = 1, est_error_var_anova = 1,
      rho = 0.4, est_rho = 0.4, G = 50,
      mu_y = c(10, 11, 12), sigma_y = 1, mu_x = 5, sigma_x = 1,
      c_weights = c(-1, 0, 1), width = 0.40, conf_level = 0.95
    )
  )
  v <- stats::setNames(res$value, res$term)
  expect_true("mean_se_ratio" %in% res$term)
  expect_equal(unname(v["est_error_var_ancova"]), 0.8)
  expect_equal(unname(v["rho"]), 0.4)
  # The population contrast implied by mu_y and the weights.
  expect_equal(unname(v["true_psi"]), sum(c(-1, 0, 1) * c(10, 11, 12)))
  expect_equal(unname(v["total_N"]), unname(3 * v["n_per_group"]))
  # Proportions on the 0 to 1 scale (this member once reported percentages).
  expect_true(all(v[c("pct_ci_less_w", "pct_ci_miss_low",
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
      estimated_var_Y = 1, estimated_cov_YX = cov_YX, estimated_cov_XX = Sigma_X,
      which_predictor = 1, w = 0.40, conf_level = 0.95,
      G = 50, print_iter = FALSE
    )
  )
  v <- stats::setNames(res$value, res$term)
  b <- (solve(Sigma_X) %*% cov_YX)[1]
  expect_equal(unname(v["true_b_j"]), b, tolerance = 1e-10)
  expect_equal(unname(v["estimated_b_j"]), b, tolerance = 1e-10)
  expect_equal(unname(v["p"]), 2)
  expect_equal(unname(v["which_predictor"]), 1)
  expect_true(all(c("mean_R2", "median_R2", "sd_R2") %in% res$term))
})

test_that("ss_aipe_rc_sensitivity() and ss_aipe_src_sensitivity() delegate to ss_aipe_reg_coef_sensitivity", {
  skip_on_cran()
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  for (fn in list(ss_aipe_rc_sensitivity, ss_aipe_src_sensitivity)) {
    set.seed(113)
    res <- suppressMessages(
      fn(true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
         estimated_var_Y = 1, estimated_cov_YX = cov_YX,
         estimated_cov_XX = Sigma_X,
         which_predictor = 1, w = 0.40, conf_level = 0.95,
         G = 50, print_iter = FALSE)
    )
    expect_s3_class(res, "data.frame")
    expect_true(all(DMAR:::.SS_AIPE_SENS_CORE_TERMS %in% res$term))
  }
})

test_that("a supplied assurance value is echoed under its own name", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_smd_sensitivity(true_delta = 0.5, estimated_delta = 0.5,
                            desired_width = 0.30, assurance = 0.90,
                            G = 5, print_iter = FALSE)
  )
  expect_equal(res$value[res$term == "assurance"], 0.90)
})

test_that("evaluating at a fixed size echoes NA for the planning value", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_aipe_smd_sensitivity(true_delta = 0.5, n_per_group = 30,
                            desired_width = 0.30,
                            G = 5, print_iter = FALSE)
  )
  v <- stats::setNames(res$value, res$term)
  expect_true(is.na(v["estimated_delta"]))
  expect_equal(unname(v["n_per_group"]), 30)
  expect_equal(unname(v["total_N"]), 60)
})
