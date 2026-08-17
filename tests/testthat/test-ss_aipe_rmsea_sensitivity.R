## ss_aipe_rmsea_sensitivity() carries out a Monte Carlo sensitivity study for
## the RMSEA AIPE planner. It fits the proposed model with lavaan (the SEM
## backend used throughout DMAR); the population RMSEA is read from
## lavaan::fitMeasures(), which is computed reliably for the large-N population
## fit. (Earlier versions used the sem package, whose summary()$RMSEA returned
## NA for this population fit, so the function could not run.) The Monte Carlo
## body is skipped on CRAN to stay within the test-time budget.

test_that("ss_aipe_rmsea_sensitivity() is exported with the documented API", {
  expect_true(is.function(ss_aipe_rmsea_sensitivity))
  expect_true(all(c("width", "model", "Sigma", "N", "conf_level", "G") %in%
                    names(formals(ss_aipe_rmsea_sensitivity))))
  expect_equal(eval(formals(ss_aipe_rmsea_sensitivity)$conf_level), 0.95)
  expect_null(eval(formals(ss_aipe_rmsea_sensitivity)$N))
})

test_that("ss_aipe_rmsea_sensitivity() recovers the population RMSEA and runs end to end", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  # True data generating model: two correlated factors (r = 0.5), loadings 0.7.
  # Build the implied population covariance Sigma = Lambda Phi Lambda' + Psi.
  Lambda <- matrix(0, 6, 2); Lambda[1:3, 1] <- 0.7; Lambda[4:6, 2] <- 0.7
  Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
  dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"

  # The population RMSEA from lavaan is finite (not NA, the sem failure mode)
  # and positive, since the one-factor model is misspecified for two-factor data.
  pop_rmsea <- unname(lavaan::fitMeasures(
    lavaan::sem(proposed, sample.cov = Sigma, sample.nobs = 1e6), "rmsea"))
  expect_true(is.finite(pop_rmsea))
  expect_gt(pop_rmsea, 0)

  set.seed(113)
  res <- ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, Sigma = Sigma, G = 10)
  expect_s3_class(res, "dmar_tbl")
  expect_true(all(c("mean_rmsea", "median_rmsea", "sd_rmsea",
                    "mean_ci_width", "median_ci_width", "sd_ci_width",
                    "pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high",
                    "total_type_I_error",
                    "suc_rep", "total_N", "df", "true_rmsea", "width",
                    "conf_level") %in% res$term))
  expect_equal(res$value[res$term == "df"], 9)
  # The reported population RMSEA matches the standalone lavaan population fit.
  expect_equal(res$value[res$term == "true_rmsea"], pop_rmsea, tolerance = 1e-6)
  expect_gt(res$value[res$term == "suc_rep"], 0)
})

test_that("ss_aipe_rmsea_sensitivity() requires named Sigma", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  expect_error(
    ss_aipe_rmsea_sensitivity(width = 0.05, model = "g =~ x1 + x2 + x3",
                              Sigma = diag(3), G = 2),
    "row and column names"
  )
})
