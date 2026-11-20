test_that("the expected width is the documented expression in each of the three cases", {
  v <- function(x, t) x$value[x$term == t]
  V <- 1; rho <- 0.10; vx <- 1; m <- 20
  tau2 <- rho * V; sig2 <- (1 - rho) * V
  # Individual-level predictor centered at its cluster mean: no design
  # effect; sigma / (sigma_W sqrt(K (m - 1) - 1)).
  r1 <- ss_aipe_mixed_effects(V, vx, true_icc = rho, width = 0.20,
                              cluster_size = m, true_icc_X = 0.25)
  K <- v(r1, "necessary_n_clusters"); nu <- K * (m - 1) - 1
  expect_equal(v(r1, "df"), nu)
  expect_equal(v(r1, "expected_width"),
               2 * qt(.975, nu) * sqrt(sig2 / (0.75 * vx)) / sqrt(nu),
               tolerance = 1e-12)
  # The planner is minimal: one fewer cluster misses the width.
  nu_less <- (K - 1) * (m - 1) - 1
  expect_gt(2 * qt(.975, nu_less) * sqrt(sig2 / 0.75) / sqrt(nu_less), 0.20)
  # Cluster-level predictor: the regression of K cluster means.
  r2 <- ss_aipe_mixed_effects(V, vx, true_icc = rho, width = 0.20,
                              cluster_size = m, predictor_level = 2)
  K2 <- v(r2, "necessary_n_clusters")
  expect_identical(attr(r2, "centering"), "grand_mean")
  expect_equal(v(r2, "df"), K2 - 2)
  expect_equal(v(r2, "expected_width"),
               2 * qt(.975, K2 - 2) * sqrt((tau2 + sig2 / m) / vx) /
                 sqrt(K2 - 2), tolerance = 1e-12)
  expect_gt(K2, K)
  # Individual-level predictor not centered: the pooled slope, the same
  # under either name, and never more clusters than the centered analysis.
  r3 <- ss_aipe_mixed_effects(V, vx, true_icc = rho, width = 0.20,
                              cluster_size = m, centering = "none",
                              true_icc_X = 0.25)
  r4 <- ss_aipe_mixed_effects(V, vx, true_icc = rho, width = 0.20,
                              cluster_size = m, centering = "grand_mean",
                              true_icc_X = 0.25)
  expect_equal(r3$value, r4$value)
  K3 <- v(r3, "necessary_n_clusters")
  expect_lte(K3, K)
  info <- (K3 * (m - 1) - 1) * 0.75 / sig2 +
    (K3 - 2) * (0.25 + 0.75 / m) / (tau2 + sig2 / m)
  expect_equal(v(r3, "expected_width"),
               2 * qt(.975, K3 * (m - 1) - 1) / sqrt(info), tolerance = 1e-12)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_mixed_effects() returns documented rows", {
  res <- ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
                                    width = 0.20, cluster_size = 20)
  expect_setequal(res$term,
                  c("necessary_n_clusters", "total_N", "width",
                    "expected_width", "df", "true_icc", "true_icc_X",
                    "cluster_size", "predictor_level"))
  expect_identical(attr(res, "centering"), "cluster_mean")
})

test_that("ss_aipe_mixed_effects() expected width matches target", {
  res <- ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
                                    width = 0.20)
  expect_lte(res$value[res$term == "expected_width"], 0.20 + 1e-8)
})

test_that("ss_aipe_mixed_effects() smaller width => more clusters", {
  r1 <- ss_aipe_mixed_effects(1, 1, true_icc = 0.10, width = 0.30)
  r2 <- ss_aipe_mixed_effects(1, 1, true_icc = 0.10, width = 0.10)
  expect_lt(r1$value[1], r2$value[1])
})

test_that("ss_aipe_mixed_effects() rejects bad inputs", {
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 1.2, width = 0.2),
               "in \\[0, 1\\)")
  expect_error(ss_aipe_mixed_effects(-1, 1, true_icc = 0.1, width = 0.2),
               "positive")
})

test_that("ss_aipe_mixed_effects() fires a type guard for every formal", {
  # Non-numeric or non-scalar values are rejected by name; no argument
  # is silently accepted.
  expect_error(ss_aipe_mixed_effects("1", 1, true_icc = 0.1, width = 0.2),
               "'true_total_var' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, c(1, 2), true_icc = 0.1, width = 0.2),
               "'true_var_X' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = "0.1", width = 0.2),
               "'true_icc' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = TRUE),
               "'width' must be a single positive number")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     cluster_size = "20"),
               "'cluster_size' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     conf_level = "0.95"),
               "'conf_level' must be a single number")
})

test_that("ss_aipe_mixed_effects() fires a range guard for every formal", {
  expect_error(ss_aipe_mixed_effects(-1, 1, true_icc = 0.1, width = 0.2),
               "'true_total_var' must be positive")
  expect_error(ss_aipe_mixed_effects(1, -1, true_icc = 0.1, width = 0.2),
               "'true_var_X' must be positive")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 1.2, width = 0.2),
               "'true_icc' must be in \\[0, 1\\)")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = -0.2),
               "'width' must be a single positive number")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     cluster_size = 1),
               "'cluster_size' must be at least 2")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     conf_level = 1.5),
               "'conf_level' must be a single number in \\(0, 1\\)")
})

test_that("ss_aipe_mixed_effects() no longer accepts a 'beta' argument", {
  # The closed-form variance of the slope does not involve the slope
  # value, so the never-read formal was removed; a supplied beta fails
  # loudly instead of being silently ignored.
  expect_error(ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.1,
                                     width = 0.2, beta = 0.3),
               "unused argument")
})

test_that("ss_aipe_mixed_effects() refuses a design that has no slope or is misdescribed", {
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     predictor_level = 2,
                                     centering = "cluster_mean"),
               class = "dmar_argument_invalid")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     predictor_level = 3),
               "'predictor_level' must be 1 or 2")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     predictor_level = 2, true_icc_X = 0.2),
               "'true_icc_X'")
  expect_error(ss_aipe_mixed_effects(1, 1, true_icc = 0.1, width = 0.2,
                                     true_icc_X = 1),
               "'true_icc_X'")
})

test_that("the sensitivity function simulates the design the planner plans", {
  skip_on_cran()
  v <- function(x, t) x$value[x$term == t]
  check <- function(...) {
    plan <- ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1,
                                  true_icc = 0.10, width = 0.40,
                                  cluster_size = 10, ...)
    sens <- ss_aipe_mixed_effects_sensitivity(
      true_total_var = 1, true_var_X = 1, true_icc = 0.10, true_reg_coef = 0.3,
      supposed_total_var = 1, supposed_var_X = 1, supposed_icc = 0.10,
      width = 0.40, cluster_size = 10, G = 400, seed = 113, ...)
    expect_equal(v(sens, "n_clusters"), v(plan, "necessary_n_clusters"))
    expect_equal(v(sens, "expected_width"), v(plan, "expected_width"))
    expect_identical(attr(sens, "centering"), attr(plan, "centering"))
    # The simulated mean width is the planned one, within Monte Carlo error
    # and the first-order error of the pooled case.
    expect_lt(abs(v(sens, "mean_ci_width") / v(plan, "expected_width") - 1),
              0.04)
    expect_lt(v(sens, "total_type_I_error"), 0.10)
    expect_lt(abs(v(sens, "mean_reg_coef") - 0.3), 0.05)
  }
  check()
  check(true_icc_X = 0.3)
  check(predictor_level = 2)
  check(predictor_level = 2, centering = "none")
  skip_if_not_installed("lme4")
  check(centering = "grand_mean", true_icc_X = 0.3)
})

test_that("the sensitivity function evaluates few clusters and refuses too few", {
  s <- ss_aipe_mixed_effects_sensitivity(
    true_total_var = 1, true_var_X = 1, true_icc = 0.1,
    specified_n_clusters = 3, cluster_size = 4, width = 1,
    predictor_level = 2, G = 20, seed = 1)
  expect_equal(s$value[s$term == "n_clusters"], 3)
  expect_error(ss_aipe_mixed_effects_sensitivity(
    true_total_var = 1, true_var_X = 1, true_icc = 0.1,
    specified_n_clusters = 2, cluster_size = 4, width = 1,
    predictor_level = 2, G = 20), class = "dmar_argument_invalid")
})
