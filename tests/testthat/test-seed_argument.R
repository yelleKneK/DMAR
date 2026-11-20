## The seed argument of the exports that draw random numbers and gained
## seed = NULL as their last formal in 1.1.0: the sensitivity functions,
## procrustes_phi(), ss_aipe_equivalence_r() (whose assurance step is a
## simulation), ss_aipe_reg_coef() and ss_aipe_src() (whose assurance step
## with assurance_method = "simulation" is a simulation), and the simulate_*
## data generators.
## ss_aipe_equivalence_smd() computes its assurance probability from the
## noncentral t distribution, draws nothing, and has no seed; the last test
## below checks that.
##
## Three properties, one loop over the calls below (tiny G, not on CRAN):
##   1. seed = NULL leaves the function on the caller's stream, so
##      set.seed(k); f() gives what it gave before the argument existed
##      (the old-against-new comparison is in the development record; here
##      the same call with and without seed = NULL is compared);
##   2. f(seed = k) is reproducible whatever the state of the caller's
##      stream, and equals set.seed(k); f();
##   3. f(seed = k) leaves the caller's .Random.seed as it found it.

# One call per function, without a seed. Functions that need lavaan are
# marked so the loop can skip them where it is absent.
.seed_calls <- function() {
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Sigma_cfa <- Lambda %*% matrix(c(1, 0.5, 0.5, 1), 2, 2) %*% t(Lambda) +
    diag(1 - 0.7^2, 6)
  dimnames(Sigma_cfa) <- list(paste0("x", 1:6), paste0("x", 1:6))
  pop_model <- "
    f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
    f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
    f2 ~ 0.5*f1
    f1 ~~ 1*f1
    f2 ~~ 0.75*f2
    y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
    y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
  "
  analysis_model <- "
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f2 ~ b*f1
  "
  list(
    ss_aipe_cliff_delta_sensitivity = quote(ss_aipe_cliff_delta_sensitivity(
      true_delta = 0.3, specified_N = 60, width = 0.5, G = 5,
      print_iter = FALSE)),
    ss_aipe_partial_r_sensitivity = quote(ss_aipe_partial_r_sensitivity(
      true_cor = 0.3, specified_N = 40, p = 2, width = 0.4, G = 5,
      print_iter = FALSE)),
    ss_aipe_semipartial_r_sensitivity = quote(ss_aipe_semipartial_r_sensitivity(
      true_r_sp = 0.3, specified_N = 40, p = 2, width = 0.4, G = 5,
      print_iter = FALSE)),
    ss_aipe_cv_sensitivity = quote(ss_aipe_cv_sensitivity(
      true_cv = 0.2, specified_N = 30, width = 0.2, G = 5,
      print_iter = FALSE)),
    ss_aipe_R2_sensitivity = quote(ss_aipe_R2_sensitivity(
      true_R2 = 0.4, specified_N = 60, width = 0.3, p = 3, G = 5,
      print_iter = FALSE)),
    ss_aipe_c_sensitivity = quote(ss_aipe_c_sensitivity(
      true_var_error = 1, supposed_var_error = 1,
      c_weights = c(-1, 1), true_contrast = 0.5, width = 1.2,
      G = 5, print_iter = FALSE)),
    ss_aipe_c_ancova_sensitivity = quote(ss_aipe_c_ancova_sensitivity(
      true_var_error_ancova = 0.8, supposed_var_error_ancova = 0.8,
      true_var_error_anova = 1, supposed_var_error_anova = 1,
      true_cor = 0.4, supposed_cor = 0.4, G = 5,
      true_mean_Y = c(10, 11, 12), true_sd_Y = 1, true_mean_X = 5, true_sd_X = 1,
      c_weights = c(-1, 0, 1), width = 0.40)),
    ss_aipe_equivalence_r_sensitivity = quote(ss_aipe_equivalence_r_sensitivity(
      true_cor = 0.1, supposed_cor = 0.1, cor_upper = 0.5, width = 0.40,
      G = 3, print_iter = FALSE)),
    ss_aipe_equivalence_smd_sensitivity = quote(ss_aipe_equivalence_smd_sensitivity(
      true_smd = 0.1, supposed_smd = 0.1, smd_upper = 0.5, width = 0.40,
      G = 3, print_iter = FALSE)),
    ss_aipe_icc_sensitivity = quote(ss_aipe_icc_sensitivity(
      true_icc = 0.8, supposed_icc = 0.8, n_raters = 4, width = 0.40,
      G = 5, print_iter = FALSE)),
    ss_aipe_indirect_effect_sensitivity = quote(ss_aipe_indirect_effect_sensitivity(
      true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.4, supposed_std_reg_coef_a = 0.4, supposed_std_reg_coef_b = 0.4,
      width = 0.40, G = 5, print_iter = FALSE)),
    ss_aipe_mixed_effects_sensitivity = quote(ss_aipe_mixed_effects_sensitivity(
      true_total_var = 1, true_var_X = 1, true_icc = 0.2,
      supposed_total_var = 1, supposed_var_X = 1,
      supposed_icc = 0.2, cluster_size = 10, width = 0.20,
      G = 3, print_iter = FALSE)),
    ss_aipe_omega_squared_sensitivity = quote(ss_aipe_omega_squared_sensitivity(
      true_omega_squared = 0.2, supposed_omega_squared = 0.2,
      df_effect = 2, width = 0.30, G = 5, print_iter = FALSE)),
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
    ss_aipe_rmsea_sensitivity = bquote(ss_aipe_rmsea_sensitivity(
      width = 0.10, model = "g =~ x1 + x2 + x3 + x4 + x5 + x6",
      true_cov_matrix = .(Sigma_cfa), specified_N = 150, G = 3)),
    ss_aipe_sc_ancova_sensitivity = quote(ss_aipe_sc_ancova_sensitivity(
      true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
      width = 0.80, true_cor = 0.4, G = 5, print_iter = FALSE)),
    ss_aipe_sc_sensitivity = quote(ss_aipe_sc_sensitivity(
      true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
      width = 0.80, G = 5, print_iter = FALSE)),
    ss_aipe_sem_path_sensitivity = bquote(ss_aipe_sem_path_sensitivity(
      model = .(analysis_model),
      supposed_cov_matrix = cov_sem(.(pop_model))$sigma_theta,
      which_path = "b", width = 0.30, specified_N = 150, G = 3)),
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
    ss_power_R2_sensitivity = quote(ss_power_R2_sensitivity(
      true_R2 = 0.30, specified_N = 50, p = 5, G = 20, print_iter = FALSE)),
    ss_power_reg_coef_sensitivity = bquote(ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
      which_predictor = 1, desired_power = 0.80, G = 20,
      print_iter = FALSE)),
    procrustes_phi = quote(procrustes_phi(
      loadings_1 = c(0.8, 0.7, 0.6, 0.5, 0.75),
      loadings_2 = c(0.75, 0.72, 0.55, 0.52, 0.7), n_resamples = 200)),
    ss_aipe_equivalence_r = quote(ss_aipe_equivalence_r(
      true_cor = 0.30, width = 0.30, assurance = 0.80)),
    # The regression coefficient planners draw only in the simulated
    # assurance step.
    ss_aipe_reg_coef = bquote(ss_aipe_reg_coef(
      true_cor_XX = .(Sigma_X), true_cor_YX = .(cov_YX), which_predictor = 1,
      width = 0.30, noncentral = TRUE, assurance = 0.80,
      assurance_method = "simulation", G = 1000)),
    ss_aipe_src = bquote(ss_aipe_src(
      true_cor_XX = .(Sigma_X), true_cor_YX = .(cov_YX), which_predictor = 1,
      width = 0.30, assurance = 0.80, assurance_method = "simulation",
      G = 1000)),
    simulate_ancova_data = quote(simulate_ancova_data(
      true_mean_Y = c(3, 5), true_mean_X = 10, true_sd_Y = 1, true_sd_X = 2, true_cor = 0.8,
      a = 2, n_per_group = 6)),
    simulate_ancova_factorial_data = quote(simulate_ancova_factorial_data(
      a = 2, b = 2, true_mean_Y = c(50, 55, 52, 57),
      true_mean_X = matrix(10, nrow = 4, ncol = 1), true_sd_Y = 8, true_sd_X = 3,
      true_cor_YX = 0.40, n_per_cell = 5)),
    simulate_longitudinal_polynomial = quote(simulate_longitudinal_polynomial(
      n_per_group = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(4, 0.25), error_variance = 1)),
    simulate_longitudinal_gompertz = quote(simulate_longitudinal_gompertz(
      n_per_group = 6, target_times = 0:6,
      fixed_parameters = c(alpha = 0.75, beta = 2, gamma = 1, zeta = 0.25),
      random_variances = c(0.01, 0.04, 0, 0), error_variance = 0.01)),
    simulate_longitudinal_logistic = quote(simulate_longitudinal_logistic(
      n_per_group = 6, target_times = 0:12,
      fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9, zeta = 10),
      random_variances = c(25, 0.25, 0, 0), error_variance = 4)),
    simulate_longitudinal_negative_exponential =
      quote(simulate_longitudinal_negative_exponential(
        n_per_group = 6, target_times = 0:10,
        fixed_parameters = c(alpha = 1, zeta = -1, gamma = 0.4),
        random_variances = c(0.01, 0.01, 0), error_variance = 0.01)),
    simulate_longitudinal_richards = quote(simulate_longitudinal_richards(
      n_per_group = 6, target_times = 0:12,
      fixed_parameters = c(alpha = 80, beta = 6, gamma = 0.9, delta = 1,
                           zeta = 10),
      random_variances = c(25, 0.25, 0, 0, 0), error_variance = 4))
  )
}
.seed_needs_lavaan <- c("ss_aipe_rmsea_sensitivity",
                        "ss_aipe_sem_path_sensitivity")

# Evaluate one call, with a seed argument appended when one is given.
# Messages and warnings that the tiny simulations give are not the subject
# here.
.seed_run <- function(call, seed) {
  if (!missing(seed)) call$seed <- seed
  suppressWarnings(suppressMessages(eval(call)))
}

test_that("each function that gained seed has it as its last formal, defaulting to NULL", {
  for (nm in names(.seed_calls())) {
    f <- formals(get(nm, envir = asNamespace("DMAR")))
    expect_identical(names(f)[length(f)], "seed", info = nm)
    expect_null(f$seed, info = nm)
  }
})

test_that("seed = NULL stays on the caller's stream, a seed is reproducible, and the caller's stream is left alone", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  calls <- .seed_calls()
  for (nm in names(calls)) {
    if (nm %in% .seed_needs_lavaan && !requireNamespace("lavaan", quietly = TRUE))
      next
    cl <- calls[[nm]]
    # 1. The default is the caller's stream.
    set.seed(4107); a0 <- .seed_run(cl)
    set.seed(4107); a1 <- .seed_run(cl, NULL)
    expect_identical(a1, a0, info = paste(nm, "seed = NULL"))
    # 2. A seed is reproducible from two different states of the caller's
    #    stream, and is the stream set.seed() gives.
    set.seed(1); b1 <- .seed_run(cl, 4107)
    set.seed(2); invisible(stats::runif(7)); b2 <- .seed_run(cl, 4107)
    expect_identical(b2, b1, info = paste(nm, "same seed twice"))
    expect_identical(b1, a0, info = paste(nm, "seed against set.seed()"))
    # 3. The caller's generator state is as it was.
    set.seed(99); before <- .Random.seed
    invisible(.seed_run(cl, 4107))
    expect_identical(.Random.seed, before, info = paste(nm, "caller's state"))
    # A different seed gives a different result, where the result is not
    # a whole number or a permutation p-value that two streams can share.
    if (!nm %in% c("procrustes_phi", "ss_aipe_equivalence_r")) {
      d <- .seed_run(cl, 4108)
      expect_false(identical(d, b1), info = paste(nm, "different seed"))
    }
  }
})

test_that("a seed that is not a single number is refused", {
  expect_error(procrustes_phi(c(0.8, 0.7, 0.6), c(0.7, 0.7, 0.5), n_resamples = 10,
                              seed = c(1, 2)),
               "'seed' must be a single number or NULL")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10, true_sd_Y = 1,
                                    true_sd_X = 2, true_cor = 0.8, a = 2, n_per_group = 6,
                                    seed = "a"),
               "'seed' must be a single number or NULL")
  expect_error(ss_aipe_r_sensitivity(true_cor = 0.4, specified_N = 50,
                                     width = 0.4, G = 2, seed = NA),
               "'seed' must be a single number or NULL")
})

test_that("the equivalence planners draw nothing without an assurance, and the SMD planner draws nothing at all", {
  set.seed(8); before <- .Random.seed
  a <- ss_aipe_equivalence_r(true_cor = 0.3, width = 0.3)
  b <- ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.5)
  b_assured <- ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.5,
                                       assurance = 0.8)
  expect_identical(.Random.seed, before)
  expect_identical(ss_aipe_equivalence_r(true_cor = 0.3, width = 0.3,
                                         seed = 5), a)
  # Its assurance probability is computed, not simulated, so it has no
  # seed (it had none in DMAR 1.0.0 either).
  expect_false("seed" %in% names(formals(ss_aipe_equivalence_smd)))
  set.seed(9)
  expect_identical(ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.5,
                                           assurance = 0.8), b_assured)
})

test_that("every sensitivity function and every data generator has a seed argument", {
  ex <- getNamespaceExports("DMAR")
  fam <- grep("_sensitivity$|^simulate_", ex, value = TRUE)
  expect_gt(length(fam), 30)
  for (nm in fam) {
    f <- formals(get(nm, envir = asNamespace("DMAR")))
    expect_true("seed" %in% names(f), info = nm)
    expect_null(f$seed, info = nm)
  }
})
