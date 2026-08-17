# ss_aipe_r_sensitivity() carries out a Monte Carlo sensitivity study for the
# Pearson correlation AIPE planner: bivariate normal samples at the true
# correlation, the Fisher's Z interval on each replication, and summary rows for
# realized width and coverage. The Monte Carlo blocks are skipped on CRAN to
# stay within the test-time budget; the planner tests in test-ss_aipe_r.R
# still run on CRAN and cover the closed-form and search surface.

test_that("ss_aipe_r_sensitivity() is exported with the documented API", {
  expect_true(is.function(ss_aipe_r_sensitivity))
  expect_true(all(c("true_rho", "estimated_rho", "width", "specified_N",
                    "conf_level", "assurance", "G", "print_iter", "save",
                    "filename") %in% names(formals(ss_aipe_r_sensitivity))))
  expect_equal(eval(formals(ss_aipe_r_sensitivity)$conf_level), 0.95)
  expect_equal(eval(formals(ss_aipe_r_sensitivity)$G), 1000)
  expect_null(eval(formals(ss_aipe_r_sensitivity)$assurance))
})

test_that("ss_aipe_r_sensitivity() validates its inputs", {
  expect_error(ss_aipe_r_sensitivity(true_rho = 0.3, width = 0.2),
               "either 'estimated_rho' or 'specified_N'")
  expect_error(ss_aipe_r_sensitivity(true_rho = 0.3, estimated_rho = 0.3,
                                     specified_N = 100, width = 0.2),
               "not both")
  expect_error(ss_aipe_r_sensitivity(true_rho = 1.2, estimated_rho = 0.3,
                                     width = 0.2),
               "in \\(-1, 1\\)")
  expect_error(ss_aipe_r_sensitivity(true_rho = 0.3, specified_N = 3,
                                     width = 0.2),
               "at least 4")
})

test_that("ss_aipe_r_sensitivity() empirical widths match the planner's promise", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  plan <- ss_aipe_r(rho = 0.30, width = 0.25)
  n    <- plan$value[plan$term == "necessary_N"]

  set.seed(113)
  res <- ss_aipe_r_sensitivity(true_rho = 0.30, estimated_rho = 0.30,
                               width = 0.25, G = 200)
  expect_s3_class(res, "dmar_tbl")
  expect_true(all(c("mean_r", "median_r", "sd_r",
                    "mean_ci_width", "median_ci_width", "sd_ci_width",
                    "pct_ci_less_w",
                    "pct_ci_miss_low", "pct_ci_miss_high",
                    "total_type_I_error",
                    "total_N", "true_rho", "estimated_rho", "width",
                    "conf_level") %in% res$term))
  # The evaluated sample size is the planner's answer.
  expect_equal(res$value[res$term == "total_N"], n)
  # The estimator centers on the population correlation. The tolerance
  # allows for the Monte Carlo error of a G = 200 mean under seed 113.
  expect_equal(res$value[res$term == "mean_r"], 0.30, tolerance = 0.05)
  # Realized mean width sits at the planner's expected width.
  expect_equal(res$value[res$term == "mean_ci_width"],
               plan$value[plan$term == "expected_width"],
               tolerance = 0.02)
  # Coverage is near the nominal 95% (overall miss rate near 0.05).
  expect_lt(res$value[res$term == "total_type_I_error"], 0.10)
})

test_that("ss_aipe_r_sensitivity() assurance raises the fraction of intervals meeting the target", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(113)
  res_50 <- ss_aipe_r_sensitivity(true_rho = 0.30, estimated_rho = 0.30,
                                  width = 0.25, G = 200)
  set.seed(113)
  res_90 <- ss_aipe_r_sensitivity(true_rho = 0.30, estimated_rho = 0.30,
                                  width = 0.25, assurance = 0.90, G = 200)
  expect_gt(res_90$value[res_90$term == "pct_ci_less_w"],
            res_50$value[res_50$term == "pct_ci_less_w"])
  expect_gte(res_90$value[res_90$term == "pct_ci_less_w"], 0.90)
})

test_that("ss_aipe_r_sensitivity() misspecification toward zero widens realized intervals", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Planned at rho = 0.60 but the population correlation is 0; the realized
  # intervals are wider than planned, and the summary shows it.
  set.seed(113)
  res <- ss_aipe_r_sensitivity(true_rho = 0, estimated_rho = 0.60,
                               width = 0.20, G = 200)
  expect_lt(res$value[res$term == "pct_ci_less_w"], 0.25)
  expect_gt(res$value[res$term == "mean_ci_width"], 0.20)
})

test_that("ss_aipe_r_sensitivity() specified_N path echoes the evaluated size", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(113)
  res <- ss_aipe_r_sensitivity(true_rho = 0.30, specified_N = 100,
                               width = 0.25, G = 50)
  expect_equal(res$value[res$term == "total_N"], 100)
  expect_true(is.na(res$value[res$term == "estimated_rho"]))
})
