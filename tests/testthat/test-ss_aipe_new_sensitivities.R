## Tests for the 10 new AIPE sensitivity functions added in this commit.
## All Monte Carlo tests behind skip_on_cran() to preserve CRAN's runtime
## budget; the planners they pair with are covered on the CRAN-side check
## by the closed-form planner tests.

# ---------- ss_aipe_c_sensitivity ----------

test_that("ss_aipe_c_sensitivity() returns the tidy schema and respects monotonicity in error_variance", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_c_sensitivity(
    true_error_variance = 4, estimated_error_variance = 4,
    c_weights = c(-1, 0, 1), width = 1,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("mean_psi", "mean_ci_width", "pct_ci_less_w",
                    "n_per_group") %in% res$term))

  # Higher true error variance should inflate the mean realized CI width.
  set.seed(113)
  res_hi <- ss_aipe_c_sensitivity(
    true_error_variance = 16, estimated_error_variance = 4,
    c_weights = c(-1, 0, 1), width = 1,
    G = 60, print_iter = FALSE
  )
  w_lo <- res$value[res$term == "mean_ci_width"]
  w_hi <- res_hi$value[res_hi$term == "mean_ci_width"]
  expect_gt(w_hi, w_lo)
})

# ---------- ss_aipe_omega_squared_sensitivity ----------

test_that("ss_aipe_omega_squared_sensitivity() returns the tidy schema and recovers the planning value", {
  skip_on_cran()
  set.seed(113)
  # suppressMessages: the planner's iterative search summarizes its
  # noncentral F lower-limit clamps in a message that is not under test.
  res <- suppressMessages(ss_aipe_omega_squared_sensitivity(
    true_omega_squared = 0.10, estimated_omega_squared = 0.10,
    df_effect = 2, width = 0.15, G = 60, print_iter = FALSE
  ))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("mean_omega_squared", "mean_ci_width", "total_N") %in% res$term))
  # The mean realized estimate should be roughly near the truth at moderate G.
  expect_lt(abs(res$value[res$term == "mean_omega_squared"] - 0.10), 0.10)
})

# ---------- ss_aipe_partial_r_sensitivity ----------

test_that("ss_aipe_partial_r_sensitivity() recovers the planning rho on average", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_partial_r_sensitivity(
    true_rho = 0.40, estimated_rho = 0.40, J = 3, width = 0.30,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_partial_r", "mean_ci_width", "total_N") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_partial_r"] - 0.40), 0.10)
})

# ---------- ss_aipe_semipartial_r_sensitivity ----------

test_that("ss_aipe_semipartial_r_sensitivity() recovers the planning value on average", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_semipartial_r_sensitivity(
    true_r_sp = 0.30, estimated_r_sp = 0.30, J = 3, width = 0.30,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_r_sp", "mean_ci_width") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_r_sp"] - 0.30), 0.10)
})

# ---------- ss_aipe_indirect_effect_sensitivity ----------

test_that("ss_aipe_indirect_effect_sensitivity() returns a tidy summary; closed_form and monte_carlo both work", {
  skip_on_cran()
  set.seed(113)
  res_s <- ss_aipe_indirect_effect_sensitivity(
    true_a = 0.4, true_b = 0.3, estimated_a = 0.4, estimated_b = 0.3,
    width = 0.30, method = "closed_form", G = 50, print_iter = FALSE
  )
  expect_s3_class(res_s, "data.frame")
  expect_true(all(c("mean_ab", "mean_ci_width", "true_ab") %in% res_s$term))
  expect_equal(res_s$value[res_s$term == "true_ab"], 0.12, tolerance = 1e-10)

  set.seed(113)
  res_m <- ss_aipe_indirect_effect_sensitivity(
    true_a = 0.4, true_b = 0.3, estimated_a = 0.4, estimated_b = 0.3,
    width = 0.30, method = "monte_carlo", B = 200, G = 30, print_iter = FALSE
  )
  expect_s3_class(res_m, "data.frame")
})

# ---------- ss_aipe_equivalence_smd_sensitivity ----------

test_that("ss_aipe_equivalence_smd_sensitivity() reports pct_equivalent inside the bounds", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_equivalence_smd_sensitivity(
    true_smd = 0.0, estimated_smd = 0.0, width = 0.40,
    delta_upper = 0.5,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_smd", "pct_equivalent",
                    "delta_lower", "delta_upper") %in% res$term))
  expect_equal(res$value[res$term == "delta_lower"], -0.5)
  expect_equal(res$value[res$term == "delta_upper"],  0.5)
})

# ---------- ss_aipe_cliff_delta_sensitivity ----------

test_that("ss_aipe_cliff_delta_sensitivity() recovers Cliff's delta on average", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_cliff_delta_sensitivity(
    true_delta = 0.30, estimated_delta = 0.30, width = 0.40,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_cliff_delta", "mean_ci_width", "n_1", "n_2") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_cliff_delta"] - 0.30), 0.15)
})

# ---------- ss_aipe_pcm_sensitivity ----------

test_that("ss_aipe_pcm_sensitivity() returns the tidy schema", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_pcm_sensitivity(
    true_variance_trend = 0.003, true_error_variance = 0.0262,
    estimated_variance_trend = 0.003, estimated_error_variance = 0.0262,
    duration = 4, frequency = 1, width = 0.025,
    G = 30, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("mean_slope_diff", "mean_ci_width", "n_per_group",
                    "n_timepoints") %in% res$term))
  expect_equal(res$value[res$term == "n_timepoints"], 5)  # f*D + 1 = 1*4 + 1
})

test_that("ss_aipe_pcm_sensitivity() evaluates the planner's two-group estimand", {
  skip_on_cran()
  # ss_aipe_pcm() sizes n PER GROUP so the expected width of the confidence
  # interval on the BETWEEN-group difference in change coefficients is at most
  # `width` (a two-sample pooled-variance t interval on per-subject estimated
  # slopes, with 2n - 2 degrees of freedom). The sensitivity simulator must
  # evaluate that same two-group estimand, so it has to (a) resolve the per-
  # group n the planner returns and (b) realize a mean CI width at or just
  # under the target with coverage near conf_level. A regression to a one-group
  # quantity (a single group, df = n - 1, standard error sd / sqrt(n)) would
  # land near half the target width and fail the width band below, which is the
  # load-bearing assertion here.

  # frequency = 1 (the Kelley & Rausch tables): f^(2p) factor is 1.
  n1 <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                    duration = 4, frequency = 1, width = 0.025)$value
  set.seed(113)
  r1 <- ss_aipe_pcm_sensitivity(
    true_variance_trend = 0.003, true_error_variance = 0.0262,
    estimated_variance_trend = 0.003, estimated_error_variance = 0.0262,
    duration = 4, frequency = 1, width = 0.025, G = 100, print_iter = FALSE
  )
  v1 <- function(t) r1$value[r1$term == t]
  expect_equal(v1("n_per_group"), n1)        # simulator uses the planned n
  expect_lt(abs(v1("mean_ci_width") - 0.025), 0.20 * 0.025)  # width on target
  expect_lt(v1("total_type_I_error"), 0.15)      # coverage ~ conf_level

  # frequency = 2 (f^(2p) = 4 at p = 1): the planner's V, and therefore n, must
  # scale with frequency^(2p). If that factor were dropped from the planner the
  # returned n would be too small and the simulator (which fits OLS on the true
  # time grid) would realize an interval wider than the target, failing the
  # width band. So this case guards the f^(2p) factor end to end. The pure
  # frequency effect is isolated by comparing the planner at f = 1 and f = 2 at
  # a common width: more frequent sampling inflates the per-unit-time slope
  # variance V by f^(2p), so the f = 2 design needs the larger n.
  n2 <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                    duration = 2, frequency = 2, width = 0.05)$value
  n1_w05 <- ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
                        duration = 4, frequency = 1, width = 0.05)$value
  expect_gt(n2, n1_w05)                          # V scales by f^(2p) = 4
  set.seed(113)
  r2 <- ss_aipe_pcm_sensitivity(
    true_variance_trend = 0.003, true_error_variance = 0.0262,
    estimated_variance_trend = 0.003, estimated_error_variance = 0.0262,
    duration = 2, frequency = 2, width = 0.05, G = 100, print_iter = FALSE
  )
  v2 <- function(t) r2$value[r2$term == t]
  expect_equal(v2("n_per_group"), n2)
  expect_lt(abs(v2("mean_ci_width") - 0.05), 0.20 * 0.05)
  expect_lt(v2("total_type_I_error"), 0.15)
})

# ---------- ss_aipe_mixed_effects_sensitivity ----------

test_that("ss_aipe_mixed_effects_sensitivity() returns the tidy schema", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  set.seed(113)
  # suppressMessages: some replications land on lme4's boundary (singular)
  # fit, which lme4 reports as a message.
  res <- suppressMessages(suppressWarnings(ss_aipe_mixed_effects_sensitivity(
    true_sigma2_y = 1, true_sigma2_x = 1, true_icc = 0.10, true_beta = 0.30,
    estimated_sigma2_y = 1, estimated_sigma2_x = 1, estimated_icc = 0.10,
    width = 0.30, cluster_size = 20,
    G = 25, print_iter = FALSE
  )))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_beta", "mean_ci_width", "n_clusters",
                    "cluster_size") %in% res$term))
})

# ---------- ss_aipe_reliability_sensitivity ----------

test_that("ss_aipe_reliability_sensitivity() recovers the planning reliability for alpha", {
  skip_on_cran()
  set.seed(113)
  res <- suppressWarnings(ss_aipe_reliability_sensitivity(
    true_reliability = 0.80, estimated_reliability = 0.80, i = 6,
    width = 0.15, estimator = "alpha", G = 30, print_iter = FALSE
  ))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_reliability", "mean_ci_width", "total_N",
                    "items") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_reliability"] - 0.80), 0.10)
})

# ---------- Input-validation sanity check (sample) ----------

test_that("the new sensitivity functions reject conflicting estimated_* and specified_N", {
  expect_error(
    ss_aipe_c_sensitivity(true_error_variance = 4, estimated_error_variance = 4,
                          n_per_group = 50, c_weights = c(-1, 0, 1), width = 1, G = 5),
    "but not both"
  )
  expect_error(
    ss_aipe_omega_squared_sensitivity(true_omega_squared = 0.1,
                                      estimated_omega_squared = 0.1,
                                      specified_N = 50, df_effect = 2, width = 0.1, G = 5),
    "but not both"
  )
  expect_error(
    ss_aipe_partial_r_sensitivity(true_rho = 0.4, estimated_rho = 0.4,
                                  specified_N = 30, J = 3, width = 0.2, G = 5),
    "but not both"
  )
})
