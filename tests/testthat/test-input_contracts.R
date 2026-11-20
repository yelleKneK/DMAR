# Inputs that once gave a wrong number, an NA table, or an error from deep
# inside base R are refused with a message that names the argument, and the
# admissible calls beside them keep their values.

val <- function(tab, t) tab$value[tab$term == t]

test_that("ci_proportion() refuses counts that are not finite", {
  expect_error(ci_proportion(1, Inf), "'N'")
  expect_error(ci_proportion(Inf, Inf), "'successes'")
  expect_equal(val(ci_proportion(0, 20), "ci_upper"), 0.1611, tolerance = 5e-4)
  expect_identical(val(ci_proportion(0, 20), "ci_lower"), 0)
})

test_that("sd_unbiased() refuses inputs with no unbiased estimate behind them", {
  expect_error(sd_unbiased(s = -1, N = 20), "'s'")
  expect_error(sd_unbiased(s = 2), "'s' and 'N'")
  expect_error(sd_unbiased(N = 20), "'s' and 'N'")
  expect_error(sd_unbiased(X = 5), "'X'")
  expect_error(sd_unbiased(X = c(1, 2, Inf)), "'X'")
  expect_error(sd_unbiased(X = c("a", "b")), "'X'")
  expect_error(sd_unbiased(s = 1, N = 1), "'N'")
  expect_error(sd_unbiased(s = 1, N = NA_real_), "'N'")
  expect_error(sd_unbiased(s = NA_real_, N = 10), "'s'")
  expect_error(sd_unbiased(s = Inf, N = 20), "'s'")
  expect_error(sd_unbiased(X = c(TRUE, FALSE, TRUE)), "'X'")
  # A missing element of a vector gives a missing row, as it did; an
  # inadmissible element refuses the call.
  expect_equal(sd_unbiased(s = c(1, NA), N = 10)$value,
               c(sd_unbiased(s = 1, N = 10)$value, NA))
  expect_equal(sd_unbiased(s = 2, N = c(10, NA))$value,
               c(sd_unbiased(s = 2, N = 10)$value, NA))
  expect_error(sd_unbiased(s = 2, N = c(10, 1)), "'N'")
  expect_error(sd_unbiased(s = c(1, -1), N = 10), "'s'")
  expect_error(sd_unbiased(s = c(1, 2, 3), N = c(10, 20)), "same length")
})

test_that("sd_unbiased() keeps its admissible calls", {
  expect_equal(sd_unbiased(s = 1, N = 10.5)$value, 1.026616, tolerance = 1e-6)
  two <- sd_unbiased(s = c(1, 2), N = c(10, 20))
  expect_identical(nrow(two), 2L)
  expect_equal(two$value, c(sd_unbiased(s = 1, N = 10)$value,
                            sd_unbiased(s = 2, N = 20)$value),
               tolerance = 1e-14)
  expect_identical(nrow(sd_unbiased(s = 1, N = c(10, 20))), 2L)
  expect_identical(sd_unbiased(s = 1, N = Inf)$value, 1)
  expect_identical(sd_unbiased(s = 0, N = 5)$value, 0)
})

test_that("cv() refuses inputs that define no coefficient of variation", {
  expect_error(cv(mean = 2, sd = -1), "'sd'")
  expect_error(cv(mean = 0, sd = 1), "'mean'")
  expect_error(cv(mean = Inf, sd = 1), "'mean'")
  expect_error(cv(mean = c(1, 2), sd = 1), "'mean'")
  expect_error(cv(cv = .2, mean = 999), "do not specify")
  expect_error(cv(cv = .2, sd = 3), "do not specify")
  expect_error(cv(mean = 100, sd = 15, cv = 0.15), "do not specify")
  expect_error(cv(mean = 100), "either 'cv' or both 'mean' and 'sd'")
  expect_error(cv(cv = "a"), "'cv'")
  expect_error(cv(cv = c(.1, .2)), "'cv'")
  expect_error(cv(cv = NA_real_), "'cv'")
  expect_error(cv(cv = .2, unbiased = TRUE), "'N'")
  expect_error(cv(cv = .2, N = 0, unbiased = TRUE), "'N'")
  expect_error(cv(cv = .2, N = 1, unbiased = TRUE), "'N'")
  expect_error(cv(cv = .2, N = Inf, unbiased = TRUE), "'N'")
  expect_error(cv(cv = .2, N = 10, unbiased = "yes"), "'unbiased'")
  expect_error(cv(cv = .2, N = 10, unbiased = 1), "'unbiased'")
  expect_error(cv(cv = .2, N = 10, unbiased = NA), "'unbiased'")
})

test_that("cv() keeps its admissible calls", {
  expect_equal(cv(mean = -10, sd = 2)$value, -0.2, tolerance = 1e-15)
  expect_equal(cv(cv = -0.2)$value, -0.2, tolerance = 1e-15)
  # N is not examined when the sample value is asked for.
  expect_equal(cv(mean = 100, sd = 15, N = 10)$value, 0.15, tolerance = 1e-15)
  expect_equal(cv(cv = .2, N = 0)$value, 0.2, tolerance = 1e-15)
  expect_equal(cv(cv = .2, N = 10, unbiased = TRUE)$value, 0.205,
               tolerance = 1e-15)
  expect_equal(cv(cv = .2, N = 10.5, unbiased = TRUE)$value,
               0.2 * (1 + 1 / 42), tolerance = 1e-15)
  expect_identical(cv(mean = 5, sd = 0)$value, 0)
})

test_that("ci_cv() refuses a coefficient, a sample, or a sample size with no interval", {
  no_interval <- "requires a positive mean and a positive standard deviation"
  invalid <- "dmar_argument_invalid"
  expect_error(ci_cv(cv = -.2, N = 20), no_interval, class = invalid)
  expect_error(ci_cv(cv = 0, N = 20), no_interval, class = invalid)
  expect_error(ci_cv(cv = c(.2, .3), N = 20), no_interval, class = invalid)
  expect_error(ci_cv(cv = NA_real_, N = 20), no_interval, class = invalid)
  expect_error(ci_cv(mean = 10, sd = -2, N = 20), no_interval, class = invalid)
  expect_error(ci_cv(mean = 10, sd = 0, N = 20), no_interval, class = invalid)
  expect_error(ci_cv(data = c(-1, -2, -3, -4)), no_interval, class = invalid)
  expect_error(ci_cv(data = c(5, 5, 5, 5)), no_interval, class = invalid)
  expect_error(ci_cv(data = c(1, 2, NA, 4)), "'data'", class = invalid)
  expect_error(ci_cv(data = 3), "'data'", class = invalid)
  expect_error(ci_cv(data = c("1", "2", "3")), "'data'", class = invalid)
  expect_error(ci_cv(cv = .2, N = Inf), "'N'", class = invalid)
  expect_error(ci_cv(cv = .2, N = 1), "'N'", class = invalid)
  expect_error(ci_cv(cv = .2, N = NA_real_), "'N'", class = invalid)
  expect_error(ci_cv(mean = 10, sd = 2, N = NA_real_), "'N'", class = invalid)
  # A mean that is not positive has one message, whether it was given as
  # 'mean' or computed from the scores, and the message says what the
  # interval requires.
  not_positive <- "The sample mean is -?[0-9.]+, and it must be positive"
  expect_error(ci_cv(mean = -10, sd = 2, N = 20), not_positive, class = invalid)
  expect_error(ci_cv(mean = 0, sd = 2, N = 20), not_positive, class = invalid)
  expect_error(ci_cv(data = c(-1, -2, -3, -4)), not_positive, class = invalid)
  expect_error(ci_cv(data = c(-2, -1, 1, 2)), not_positive, class = invalid)
  expect_error(ci_cv(mean = -10, sd = 2, N = 20), no_interval)
  # A mean that is not one finite number is a malformed argument.
  expect_error(ci_cv(mean = NA_real_, sd = 2, N = 20),
               "'mean' must be a single finite number", class = invalid)
  expect_error(ci_cv(mean = 5, sd = "a", N = 30),
               "'sd' must be a single finite number", class = invalid)
})

test_that("ci_cv() keeps its admissible calls, a fractional n included", {
  for (n in c(1.5, 20, 20.5)) {
    # At n = 1.5 the upper limit is infinite, which is part of the exact
    # interval: nothing is signaled.
    out <- expect_no_condition(ci_cv(cv = .2, N = n))
    expect_identical(is.infinite(val(out, "ci_upper")), n < 2)
    expect_identical(out$term, c("ci_lower", "cv", "cv_unbiased",
                                 "ci_upper"))
    expect_identical(val(out, "cv"), .2)
    expect_lt(val(out, "ci_lower"), .2)
    expect_gt(val(out, "ci_upper"), .2)
  }
  set.seed(113)
  X <- rnorm(15, 5, 1)
  a <- ci_cv(data = X)
  b <- ci_cv(mean = mean(X), sd = sd(X), N = 15)
  d <- ci_cv(cv = sd(X) / mean(X), N = 15)
  expect_equal(a$value, b$value, tolerance = 1e-12)
  expect_equal(a$value, d$value, tolerance = 1e-12)
})

test_that("power_equivalence_c() refuses missing and infinite inputs", {
  call_with <- function(...) {
    args <- utils::modifyList(
      list(c_weights = c(1, -1), n_per_group = 20, true_sd = 1, contrast_upper = .5),
      list(...))
    do.call(power_equivalence_c, args)
  }
  expect_error(call_with(true_contrast = NA_real_), "'true_contrast'")
  expect_error(call_with(true_contrast = Inf), "'true_contrast'")
  expect_error(call_with(n_per_group = Inf), "'n_per_group'")
  expect_error(call_with(n_per_group = c(20, NA)), "'n_per_group'")
  expect_error(call_with(contrast_upper = Inf), "'contrast_upper'")
  expect_error(call_with(contrast_lower = Inf), "'contrast_lower'")
  expect_error(call_with(true_sd = NA_real_), "'true_sd'")
  expect_error(call_with(true_sd = Inf), "'true_sd'")
  expect_error(call_with(alpha_level = NA_real_), "'alpha_level'")
  expect_error(call_with(df_error = Inf), "degrees of freedom")
  expect_error(call_with(df_error = NA_real_), "degrees of freedom")
  expect_error(power_equivalence_c(N = Inf, true_sd = 1, contrast_lower = .5,
                                   side = "noninferiority", design = "paired"),
               "'N'")
  # A fractional n_per_group and a fractional df_error stay admissible.
  expect_equal(call_with(n_per_group = 20.5)$value, 0.03502764, tolerance = 1e-7)
  expect_equal(call_with(df_error = 17.3)$value, 0.04811175, tolerance = 1e-7)
})

test_that("every exported function with a G argument checks it as a count", {
  ns <- asNamespace("DMAR")
  with_G <- Filter(function(f) {
    obj <- get(f, envir = ns)
    is.function(obj) && "G" %in% names(formals(obj))
  }, sort(getNamespaceExports("DMAR")))
  expect_gte(length(with_G), 43L)
  checked <- vapply(with_G, function(f) {
    any(grepl(".dmar_check_count(G, \"G\"", deparse(body(get(f, envir = ns))),
              fixed = TRUE))
  }, logical(1))
  expect_identical(names(checked)[!checked], character(0))
})

test_that("a G that is not a count of replications is refused before any simulation", {
  bad <- list(0, 2.5, -3, NA_real_, Inf, c(5, 6), "5", NULL, TRUE)
  for (G in bad) {
    expect_error(ss_aipe_smd_sensitivity(true_smd = .5, specified_n_per_group = 20,
                                         G = G),
                 class = "dmar_count_invalid")
    expect_error(ss_aipe_r_sensitivity(true_cor = .3, width = .2,
                                       specified_N = 50, G = G),
                 class = "dmar_count_invalid")
    expect_error(ss_aipe_cv_sensitivity(true_cv = .2, specified_N = 20, G = G),
                 class = "dmar_count_invalid")
  }
  expect_error(ss_aipe_smd_sensitivity(true_smd = .5, specified_n_per_group = 20, G = 0),
               "'G' must be a single whole number of at least 1")
  expect_error(ss_aipe_indirect_effect(true_std_reg_coef_a = .4, true_std_reg_coef_b = .4, width = .2,
                                       method = "monte_carlo", G = 9),
               "'G' must be a single whole number of at least 10")
  expect_error(ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                                    true_sd = 15.67, true_means = c(70, 65),
                                    G = 1),
               "'G' must be a single whole number of at least 2")
  expect_error(ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                                    true_sd = 15.67, true_means = c(70, 65),
                                    G = 2.5),
               class = "dmar_count_invalid")
  # The two wrappers leave a missing G to the function they call.
  expect_error(ss_aipe_rc_sensitivity(G = 2.5), class = "dmar_count_invalid")
  expect_error(ss_aipe_src_sensitivity(G = 2.5), class = "dmar_count_invalid")
})

test_that("ss_aipe_rmsea_sensitivity() refuses a G that is not a count", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
  dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  for (G in list(0, 2.5, NA_real_, Inf, c(5, 6), "5")) {
    expect_error(
      ss_aipe_rmsea_sensitivity(width = 0.05,
                                model = "g =~ x1 + x2 + x3 + x4 + x5 + x6",
                                true_cov_matrix = Sigma, specified_N = 200, G = G),
      class = "dmar_count_invalid")
  }
})

test_that("a whole G given as a double runs as before", {
  set.seed(113)
  a <- ss_aipe_cv_sensitivity(true_cv = .2, specified_N = 20, G = 5)
  set.seed(113)
  b <- ss_aipe_cv_sensitivity(true_cv = .2, specified_N = 20, G = 5L)
  expect_identical(a$value, b$value)
  expect_null(attr(a, "n_without_interval"))
  expect_identical(attr(a, "simulation")$failed, 0L)
  expect_identical(val(a, "total_type_I_error"),
                   val(a, "pct_ci_miss_low") + val(a, "pct_ci_miss_high"))
})

test_that("ss_aipe_cv_sensitivity() leaves out replications whose sample mean is not positive", {
  skip_on_cran()
  # P(sample mean <= 0) is pnorm(-sqrt(10)), about 7.8e-4, a replication; in
  # 3000 replications at least one occurs with probability 0.90, and this
  # seed has some.
  set.seed(7)
  seen <- list()
  # With a coefficient of variation this large many of the intervals that
  # exist have an infinite upper limit, which the function reports in a
  # warning of its own; that warning is not what is tested here.
  out <- suppressWarnings(withCallingHandlers(
    ss_aipe_cv_sensitivity(true_cv = 1, specified_N = 10, G = 3000),
    dmar_replications_without_interval = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }))
  expect_length(seen, 1L)
  expect_match(conditionMessage(seen[[1L]]), "not positive")
  n_out <- seen[[1L]]$n_without_interval
  expect_gt(n_out, 0L)
  expect_identical(seen[[1L]]$G, 3000L)
  expect_identical(attr(out, "simulation")$failed, n_out)
  expect_s3_class(out, "dmar_tbl")
  width_terms <- c("mean_ci_width", "median_ci_width", "sd_ci_width",
                   "pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high",
                   "total_type_I_error")
  # The intervals that exist can still have an infinite upper limit. The
  # mean and SD of the width are over the finite widths, so they are
  # finite; the median is over every width.
  expect_true(all(!is.na(out$value[out$term %in%
    c("median_ci_width", "pct_ci_miss_low", "pct_ci_miss_high",
      "total_type_I_error")])))
  expect_true(all(is.finite(out$value[out$term %in%
    c("mean_ci_width", "sd_ci_width")])))
  expect_gt(attr(out, "simulation")$n_infinite_width, 0L)
  expect_identical(attr(out, "simulation")$failed, n_out)
  expect_identical(attr(out, "simulation")$attempted, 3000L)
  expect_true(all(width_terms %in% out$term))
  expect_equal(val(out, "total_type_I_error"),
               val(out, "pct_ci_miss_low") + val(out, "pct_ci_miss_high"),
               tolerance = 1e-12)
  expect_true(all(is.finite(out$value[out$term %in%
    c("mean_cv", "median_cv", "sd_cv")])))
})

test_that("each of the 45 functions with a G argument refuses a fractional G", {
  # One call of each function with G = 2.5, in place of a search of the
  # source for the check. The call is refused before anything is simulated,
  # so each takes no time; the other arguments are those of a call that
  # runs.
  cov_YX  <- c(0.3, 0.3)
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), 2)
  analysis_model <- "m ~ a * x\ny ~ b * m\nab := a * b"
  pop_model <- "m ~ 0.4 * x\ny ~ 0.4 * m\nx ~~ 1 * x\nm ~~ 1 * m\ny ~~ 1 * y"
  one_factor <- "f =~ x1 + x2 + x3 + x4"
  Sigma <- matrix(0.4, 4, 4, dimnames = list(paste0("x", 1:4),
                                             paste0("x", 1:4)))
  diag(Sigma) <- 1
  path_Sigma <- matrix(c(1, 0.4, 0.3, 0.4, 1, 0.5, 0.3, 0.5, 1), 3,
                       dimnames = list(c("x", "m", "y"), c("x", "m", "y")))
  G <- 2.5
  calls <- list(
    ss_aipe_c_ancova_sensitivity = quote(ss_aipe_c_ancova_sensitivity(
      true_var_error_ancova = 30, supposed_var_error_ancova = 30, true_cor = 0.2,
      true_mean_Y = c(10, 12, 15, 13), true_mean_X = 2, G = G, true_sd_X = 1.3, true_sd_Y = 2,
      c_weights = c(1, 0, -1, 0), width = 3)),
    ss_aipe_c_sensitivity = quote(ss_aipe_c_sensitivity(
      true_var_error = 4, supposed_var_error = 4,
      c_weights = c(-1, 0, 1), width = 1, G = G)),
    ss_aipe_cliff_delta_sensitivity = quote(ss_aipe_cliff_delta_sensitivity(
      true_delta = 0.3, supposed_delta = 0.3, width = 0.3, G = G)),
    ss_aipe_composite_sem = quote(ss_aipe_composite_sem(
      model = analysis_model, true_model = pop_model, parameters = "ab",
      width = c(ab = 0.25), specified_N = 200, G = G, seed = 113)),
    ss_aipe_crossed_sensitivity = quote(ss_aipe_crossed_sensitivity(
      true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
      supposed_var_unit = 0.5, supposed_var_rater = 0.1,
      supposed_var_residual = 0.4, width = 0.4, G = G, seed = 113)),
    ss_aipe_cv_sensitivity = quote(ss_aipe_cv_sensitivity(
      true_cv = 0.25, supposed_cv = 0.25, width = 0.1, G = G)),
    ss_aipe_equivalence_r_sensitivity = quote(
      ss_aipe_equivalence_r_sensitivity(true_cor = 0, supposed_cor = 0,
                                        width = 0.3, cor_upper = 0.2,
                                        G = G)),
    ss_aipe_equivalence_smd_sensitivity = quote(
      ss_aipe_equivalence_smd_sensitivity(true_smd = 0, supposed_smd = 0,
                                          width = 0.3, smd_upper = 0.2,
                                          G = G)),
    ss_aipe_icc_sensitivity = quote(ss_aipe_icc_sensitivity(
      true_icc = 0.7, supposed_icc = 0.7, n_raters = 3, width = 0.4, G = G)),
    ss_aipe_indirect_effect = quote(ss_aipe_indirect_effect(
      true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.4, width = 0.2, method = "monte_carlo", G = G)),
    ss_aipe_indirect_effect_sensitivity = quote(
      ss_aipe_indirect_effect_sensitivity(true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.3,
                                          supposed_std_reg_coef_a = 0.4,
                                          supposed_std_reg_coef_b = 0.3, width = 0.4,
                                          method = "closed_form", G = G)),
    ss_aipe_mean_sensitivity = quote(ss_aipe_mean_sensitivity(
      true_sd = 10, supposed_sd = 10, width = 4, G = G)),
    ss_aipe_mixed_effects_sensitivity = quote(
      ss_aipe_mixed_effects_sensitivity(true_total_var = 1, true_var_X = 1,
                                        true_icc = 0.1, true_reg_coef = 0.3,
                                        specified_n_clusters = 25, cluster_size = 10,
                                        width = 0.4, G = G)),
    ss_aipe_omega_squared_sensitivity = quote(
      ss_aipe_omega_squared_sensitivity(true_omega_squared = 0.1,
                                        supposed_omega_squared = 0.1,
                                        df_effect = 2, width = 0.1, G = G)),
    ss_aipe_partial_r_sensitivity = quote(ss_aipe_partial_r_sensitivity(
      true_cor = 0.4, supposed_cor = 0.4, p = 4, width = 0.4, G = G)),
    ss_aipe_pcm_sensitivity = quote(ss_aipe_pcm_sensitivity(
      true_var_trend = 0.003, true_var_error = 0.0262,
      supposed_var_trend = 0.003, supposed_var_error = 0.0262,
      duration = 4, frequency = 1, width = 0.05, G = G)),
    ss_aipe_r_sensitivity = quote(ss_aipe_r_sensitivity(
      true_cor = 0.3, supposed_cor = 0.3, width = 0.4, G = G)),
    ss_aipe_R2 = quote(ss_aipe_R2(
      true_R2 = 0.5, width = 0.1, p = 5, verify_ss = TRUE, G = G)),
    ss_aipe_R2_sensitivity = quote(ss_aipe_R2_sensitivity(
      true_R2 = 0.5, supposed_R2 = 0.4, width = 0.1, p = 5, G = G)),
    ss_aipe_rc_sensitivity = quote(ss_aipe_rc_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      supposed_var_Y = 1, supposed_cov_YX = cov_YX,
      supposed_cov_XX = Sigma_X, which_predictor = 1, width = 0.2, G = G)),
    ss_aipe_reg_coef_sensitivity = quote(ss_aipe_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      supposed_var_Y = 1, supposed_cov_YX = cov_YX,
      supposed_cov_XX = Sigma_X, which_predictor = 1, width = 0.2, G = G)),
    ss_aipe_reliability_sensitivity = quote(ss_aipe_reliability_sensitivity(
      true_reliability = 0.8, supposed_reliability = 0.8, n_items = 4,
      width = 0.15, estimator = "alpha", G = G)),
    ss_aipe_rmsea_sensitivity = quote(ss_aipe_rmsea_sensitivity(
      width = 0.05, model = one_factor, true_cov_matrix = Sigma, specified_N = 200, G = G)),
    ss_aipe_sc_ancova_sensitivity = quote(ss_aipe_sc_ancova_sensitivity(
      true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
      width = 0.5, true_cor = 0.4, G = G)),
    # The planner draws only for its simulated assurance step and checks G
    # only then.
    ss_aipe_sc_ancova = quote(ss_aipe_sc_ancova(
      true_std_contrast = 0.5, c_weights = c(1, -1), width = 0.5,
      assurance = 0.8, assurance_method = "simulation", G = G)),
    ss_aipe_sc_sensitivity = quote(ss_aipe_sc_sensitivity(
      true_std_contrast = 0.5, supposed_std_contrast = 0.5, c_weights = c(-1, 0, 1),
      width = 0.4, G = G)),
    ss_aipe_sem_path_sensitivity = quote(ss_aipe_sem_path_sensitivity(
      model = analysis_model, supposed_cov_matrix = path_Sigma,
      true_cov_matrix = path_Sigma, which_path = "b", width = 0.3,
      G = G)),
    # The planner draws only for its assurance step and checks G only then.
    ss_aipe_semipartial_r = quote(ss_aipe_semipartial_r(
      true_r_sp = 0.3, p = 2, width = 0.2, assurance = 0.8, G = G)),
    # The regression coefficient planners draw only for the simulated
    # assurance step and check G only then.
    ss_aipe_reg_coef = quote(ss_aipe_reg_coef(
      true_cor_XX = Sigma_X, true_cor_YX = cov_YX, which_predictor = 1,
      width = 0.3, noncentral = TRUE, assurance = 0.8,
      assurance_method = "simulation", G = G)),
    ss_aipe_src = quote(ss_aipe_src(
      true_cor_XX = Sigma_X, true_cor_YX = cov_YX, which_predictor = 1,
      width = 0.3, assurance = 0.8, assurance_method = "simulation", G = G)),
    ss_aipe_semipartial_r_sensitivity = quote(
      ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3,
                                        supposed_r_sp = 0.3, p = 4,
                                        width = 0.2, G = G)),
    ss_aipe_sm_sensitivity = quote(ss_aipe_sm_sensitivity(
      true_sm = 1.5, supposed_sm = 1.5, width = 0.5, G = G)),
    ss_aipe_smd_sensitivity = quote(ss_aipe_smd_sensitivity(
      true_smd = 0.5, supposed_smd = 0.5, width = 0.3, G = G)),
    ss_aipe_src_sensitivity = quote(ss_aipe_src_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      supposed_var_Y = 1, supposed_cov_YX = cov_YX,
      supposed_cov_XX = Sigma_X, which_predictor = 1, width = 0.2, G = G)),
    ss_aipe_two_phase_sensitivity = quote(ss_aipe_two_phase_sensitivity(
      true_sd_Y = 0.43, true_cor = 0.45, supposed_sd_Y = 0.43,
      supposed_cor = 0.55, n_pool = 20000, width = 0.1, G = G,
      seed = 113)),
    ss_power_composite_sem = quote(ss_power_composite_sem(
      model = analysis_model, true_model = pop_model, specified_N = 200, G = G,
      seed = 113)),
    ss_power_R2_sensitivity = quote(ss_power_R2_sensitivity(
      true_R2 = 0.3, supposed_R2 = 0.3, desired_power = 0.8, p = 5, G = G)),
    ss_power_reg_coef_sensitivity = quote(ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 1, desired_power = 0.8, G = G)),
    ss_seq_c_sensitivity = quote(ss_seq_c_sensitivity(
      c_weights = c(1, -1), half_width = 2.5, true_sd = 15.67, G = G,
      seed = 113)),
    ss_aipe_reliability = quote(ss_aipe_reliability(
      model = "Parallel", type = "Normal Theory", width = 0.1, n_items = 6,
      true_item_covariance = 0.6, true_residual_var = 0.4, assurance = 0.85,
      G = G)),
    ss_aipe_crd_es_n_clusters_fixed_width = quote(
      ss_aipe_crd_es_n_clusters_fixed_width(
        width = 0.4, cluster_size = 20, true_es = 0.5, true_icc = 0.1,
        pr_treat = 0.5, G = G)),
    ss_aipe_crd_es_n_individuals_fixed_width = quote(
      ss_aipe_crd_es_n_individuals_fixed_width(
        width = 0.4, n_clusters = 100, true_es = 0.5, true_icc = 0.1,
        pr_treat = 0.5, G = G)),
    ss_aipe_crd_es_n_clusters_fixed_budget = quote(
      ss_aipe_crd_es_n_clusters_fixed_budget(
        budget = 1000, cluster_size = 20, cluster_cost = 0,
        individual_cost = 1, G = G)),
    ss_aipe_crd_es_n_individuals_fixed_budget = quote(
      ss_aipe_crd_es_n_individuals_fixed_budget(
        budget = 1000, n_clusters = 200, cluster_cost = 0,
        individual_cost = 1, G = G)),
    ss_aipe_crd_es_both_fixed_budget = quote(
      ss_aipe_crd_es_both_fixed_budget(
        budget = 1000, cluster_cost = 5, individual_cost = 1, true_es = 0.5,
        true_icc = 0.1, pr_treat = 0.5, G = G)),
    ss_aipe_crd_es_both_fixed_width = quote(
      ss_aipe_crd_es_both_fixed_width(
        width = 0.5, cluster_cost = 5, individual_cost = 1, true_es = 0.5,
        true_icc = 0.25, pr_treat = 0.5, G = G))
  )
  # The list is every exported function with the argument.
  ns <- asNamespace("DMAR")
  with_G <- Filter(function(f) {
    obj <- get(f, envir = ns)
    is.function(obj) && "G" %in% names(formals(obj))
  }, getNamespaceExports(ns))
  expect_setequal(names(calls), with_G)
  for (f in names(calls)) {
    expect_error(eval(calls[[f]]), class = "dmar_count_invalid", info = f)
  }
  # The smallest count: 2 for the two functions whose summaries need two
  # replications, 10 for the three planners.
  G <- 1
  expect_error(eval(calls$ss_aipe_crossed_sensitivity),
               class = "dmar_count_invalid")
  expect_error(eval(calls$ss_seq_c_sensitivity),
               class = "dmar_count_invalid")
  G <- 9
  for (f in c("ss_aipe_indirect_effect", "ss_aipe_composite_sem",
              "ss_power_composite_sem")) {
    expect_error(eval(calls[[f]]), class = "dmar_count_invalid", info = f)
  }
})

test_that("every summary row of ss_aipe_cv_sensitivity() describes the replications that have an interval", {
  skip_on_cran()
  # The file holds every replication. All ten statistic rows are over the
  # rows that have limits (a positive sample mean); the mean and SD of the
  # width are over the finite widths among them.
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  set.seed(7)
  out <- suppressWarnings(
    ss_aipe_cv_sensitivity(true_cv = 1, specified_N = 10, G = 3000,
                           width = 2, filename = f))
  raw <- utils::read.csv(f)
  expect_identical(nrow(raw), 3000L)
  has <- !is.na(raw$ci_lower)
  expect_identical(sum(!has), attr(out, "simulation")$failed)
  expect_gt(sum(!has), 0L)
  expect_equal(val(out, "mean_cv"), mean(raw$cv[has]), tolerance = 1e-12)
  expect_equal(val(out, "median_cv"), stats::median(raw$cv[has]),
               tolerance = 1e-12)
  expect_equal(val(out, "sd_cv"), stats::sd(raw$cv[has]), tolerance = 1e-12)
  expect_gt(val(out, "mean_cv"), 0)
  expect_equal(val(out, "pct_ci_miss_low"),
               mean(raw$ci_lower[has] >= 1), tolerance = 1e-12)
  expect_equal(val(out, "pct_ci_miss_high"),
               mean(raw$ci_upper[has] <= 1), tolerance = 1e-12)
  w <- raw$ci_upper[has] - raw$ci_lower[has]
  expect_gt(sum(is.infinite(w)), 0L)
  expect_identical(attr(out, "simulation")$n_infinite_width,
                   sum(is.infinite(w)))
  expect_equal(val(out, "median_ci_width"), stats::median(w),
               tolerance = 1e-12)
  expect_equal(val(out, "mean_ci_width"), mean(w[is.finite(w)]),
               tolerance = 1e-12)
  expect_equal(val(out, "sd_ci_width"), stats::sd(w[is.finite(w)]),
               tolerance = 1e-12)
  # An infinite width does not meet the target; the row is over the
  # replications with an interval and equals the record's entry, whose
  # denominator is the number of replications with an interval.
  expect_equal(val(out, "pct_ci_less_width"), mean(w <= 2), tolerance = 1e-12)
  sim <- attr(out, "simulation")
  expect_equal(sim$width_target_probability,
               val(out, "pct_ci_less_width"), tolerance = 1e-12)
  expect_identical(sim$returned, sum(has))
  expect_identical(sim$denominator, sum(has))
  # The estimate of a replication without an interval is the one the sample
  # gave, sign included, and it stays in the file.
  expect_true(all(raw$cv[!has] <= 0))
})

test_that("ss_aipe_cv_sensitivity() reports the replications with an infinite limit once, with their count", {
  # Seed 88, true_cv = 3, N = 5, G = 500: 114 replications have a mean
  # that is not positive and no interval, and 346 of the other 386 have an
  # interval whose upper limit is infinite. ci_cv() signals nothing for
  # such an interval; the function says once how many there were and how
  # they enter the rows.
  seen <- list()
  set.seed(88)
  out <- withCallingHandlers(
    ss_aipe_cv_sensitivity(true_cv = 3, specified_N = 5, width = .5, G = 500),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  classes <- vapply(seen, function(w) class(w)[1L], "")
  expect_identical(sort(classes), c("dmar_replications_with_infinite_limit",
                                    "dmar_replications_without_interval"))
  infinite <- seen[[which(classes == "dmar_replications_with_infinite_limit")]]
  expect_identical(infinite$code, "replications_with_infinite_limit")
  expect_identical(infinite$G, 500L)
  expect_identical(infinite$n_replications, 346L)
  expect_identical(conditionMessage(infinite), paste0(
    "In 346 of the 500 replications the upper limit of the confidence ",
    "interval was infinite; such an interval counts as not meeting the ",
    "width target and is left out of mean_ci_width and sd_ci_width ",
    "(median_ci_width is over every width)."))
  expect_identical(attr(out, "simulation")$failed, 114L)
  expect_gt(val(out, "mean_cv"), 0)
  expect_identical(attr(out, "simulation")$n_infinite_width, 346L)
  expect_identical(val(out, "median_ci_width"), Inf)
  expect_true(is.finite(val(out, "mean_ci_width")))
  expect_identical(val(out, "pct_ci_less_width"), 0)
  expect_output(print(out),
                "Replications: 500 attempted, 386 with an interval, 114 without",
                fixed = TRUE)
  expect_output(print(out),
                "Intervals with an infinite upper limit: 346 of 386",
                fixed = TRUE)
  # true_cv = 1, N = 10: every replication has an interval here, some of
  # them unbounded above, and the mean width is finite. Without a width
  # target the warning does not speak of one.
  set.seed(88)
  w <- tryCatch(ss_aipe_cv_sensitivity(true_cv = 1, specified_N = 10, G = 500),
                warning = function(w) w)
  expect_s3_class(w, "dmar_replications_with_infinite_limit")
  expect_identical(w$n_replications, 96L)
  expect_no_match(conditionMessage(w), "width target")
  set.seed(88)
  out <- suppressWarnings(
    ss_aipe_cv_sensitivity(true_cv = 1, specified_N = 10, width = .5, G = 500))
  expect_true(is.finite(val(out, "mean_ci_width")))
  expect_identical(attr(out, "simulation")$n_infinite_width, 96L)
  expect_identical(attr(out, "simulation")$failed, 0L)
  # No replication with an infinite limit, no warning.
  expect_no_warning(
    ss_aipe_cv_sensitivity(true_cv = .2, specified_N = 20, G = 20, seed = 88))
})

test_that("ss_aipe_cv_sensitivity() refuses inputs that define no simulation", {
  call_with <- function(...) {
    args <- utils::modifyList(list(true_cv = .25, specified_N = 30, G = 5),
                              list(...))
    do.call(ss_aipe_cv_sensitivity, args)
  }
  refused <- function(expr, arg, class = "dmar_argument_invalid") {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, class)
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$arg, arg)
  }
  # The mean the scores are generated with: one positive number.
  refused(call_with(true_mean = -100), "true_mean")
  refused(call_with(true_mean = 0), "true_mean")
  refused(call_with(true_mean = c(1, 2)), "true_mean")
  refused(call_with(true_mean = Inf), "true_mean")
  # The two coefficients of variation.
  refused(call_with(true_cv = c(.2, .3)), "true_cv")
  refused(call_with(true_cv = 0), "true_cv")
  refused(call_with(true_cv = NULL), "true_cv")
  refused(ss_aipe_cv_sensitivity(true_cv = .25, supposed_cv = c(.2, .3),
                                 width = .1, G = 5), "supposed_cv")
  refused(ss_aipe_cv_sensitivity(true_cv = .25, supposed_cv = -1, width = .1,
                                 G = 5), "supposed_cv")
  # The sample size is a count of at least 4, the smallest size the
  # planner returns; a character string does not turn the value column
  # into text, and a fraction is not echoed.
  refused(call_with(specified_N = "30"), "specified_N", "dmar_count_invalid")
  refused(call_with(specified_N = 50.5), "specified_N", "dmar_count_invalid")
  refused(call_with(specified_N = 3), "specified_N", "dmar_count_invalid")
  refused(call_with(specified_N = c(20, 30)), "specified_N", "dmar_count_invalid")
  # The width and the assurance are checked with a specified size too.
  refused(call_with(width = -1), "width")
  refused(call_with(width = c(.1, .2)), "width")
  refused(call_with(assurance = 5), "assurance", "dmar_assurance_invalid")
  refused(call_with(assurance = "high"), "assurance", "dmar_assurance_invalid")
  # The admissible calls beside them.
  out <- call_with(specified_N = 4, seed = 1)
  expect_true(is.numeric(out$value))
  expect_identical(val(out, "total_N"), 4)
  expect_true(is.nan(val(out, "pct_ci_less_width")))
  out <- call_with(width = .2, assurance = .8, seed = 1)
  expect_identical(val(out, "assurance"), .8)
  expect_identical(val(out, "width"), .2)
})
