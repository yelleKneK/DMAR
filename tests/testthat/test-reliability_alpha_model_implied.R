# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for the model implied estimator of reliability_alpha(): coefficient
# alpha read off the tau-equivalent single-factor model fit by maximum
# likelihood. The analytic (closed-form equation) estimator is tested in
# test-reliability_alpha.R. The two were separate functions until the
# estimator argument merged them.

# Tests for reliability_alpha(estimator = "model_implied"): the tau-equivalent model
# estimate of coefficient alpha with model-based intervals. Anchored to
# MBESS::ci.reliability(type = "alpha-analytic").

skip_if_not_installed("lavaan")

set.seed(113)
J <- 6
N <- 250
loadings <- seq(0.4, 0.8, length.out = J)
eta <- rnorm(N)
errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) + errors
colnames(items_continuous) <- paste0("y", seq_len(J))
S_continuous <- cov(items_continuous)


test_that("the point estimate matches MBESS type = 'alpha-analytic'", {
  r <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                                  ci_method = "none")
  est <- r$value[r$term == "estimate"]
  expect_identical(attr(r, "coefficient"), "alpha")
  expect_identical(attr(r, "estimator"), "model_implied")

  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(est, 0.8086414182034473, tolerance = 1e-6)
})

test_that("the ML Wald interval is close to the MBESS ml interval", {
  r <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                                  ci_method = "ml")
  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(r$value[r$term == "ci_lower"], 0.7720506085969179,
               tolerance = 2e-3)
  expect_equal(r$value[r$term == "ci_upper"], 0.8452322278099766,
               tolerance = 2e-3)
})

test_that("the likelihood interval matches MBESS interval.type = 'll'", {
  skip_on_cran()  # one likelihood profile; the point estimate anchor above runs on CRAN
  r <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                                  ci_method = "likelihood")
  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(r$value[r$term == "ci_lower"], 0.7691511985556565,
               tolerance = 5e-4)
  expect_equal(r$value[r$term == "ci_upper"], 0.8426612898906286,
               tolerance = 5e-4)
})

test_that("intervals bracket the estimate and covariance input works", {
  for (m in c("mlr", "ml")) {
    r <- reliability_alpha(estimator = "model_implied", data = items_continuous, ci_method = m)
    est <- r$value[r$term == "estimate"]
    expect_lte(r$value[r$term == "ci_lower"], est)
    expect_gte(r$value[r$term == "ci_upper"], est)
  }
  r_S <- reliability_alpha(estimator = "model_implied", S = S_continuous, N = N,
                                    ci_method = "ml")
  r_d <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                                    ci_method = "ml")
  expect_equal(r_S$value[r_S$term == "estimate"],
               r_d$value[r_d$term == "estimate"], tolerance = 1e-6)
})

test_that("the model estimate is near, but not identical to, the sample alpha", {
  r_model  <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                                         ci_method = "none")
  r_sample <- reliability_alpha(data = items_continuous,
                                ci_method = "none")
  a_model  <- r_model$value[r_model$term == "estimate"]
  a_sample <- r_sample$value[r_sample$term == "estimate"]
  expect_lt(abs(a_model - a_sample), 0.02)
  expect_false(isTRUE(all.equal(a_model, a_sample, tolerance = 1e-10)))
})

test_that("raw-data-only methods error with covariance input", {
  for (m in c("mlr", "adf", "percentile", "bca")) {
    expect_error(
      reliability_alpha(estimator = "model_implied", S = S_continuous, N = N, ci_method = m),
      regexp = "requires raw 'data'"
    )
  }
})

test_that("reliability() dispatches the model implied estimator with the S fallback", {
  r_wrap <- reliability(data = items_continuous, type = "alpha",
                       estimator = "model_implied",
                        ci_method = "ml")
  r_direct <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                                         ci_method = "ml")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "alpha")
  expect_equal(attr(r_wrap, "estimator"), "model_implied")

  r_S <- reliability(S = S_continuous, N = N, type = "alpha",
                     estimator = "model_implied")
  expect_equal(attr(r_S, "ci_method"), "ml")
})

test_that("the tau-equivalent bootstrap refits by ML, whose estimate equals MLR's", {
  # ML and MLR share their estimates to the last bit, so the refits by
  # ML (since 1.1.0) reproduce the MLR refits of 1.0.1.
  ml  <- .omega_fit_cfa(data = items_continuous, equal_loadings = TRUE,
                        estimator = "ML", se = "none")
  mlr <- .omega_fit_cfa(data = items_continuous, equal_loadings = TRUE,
                        estimator = "MLR", se = "none")
  expect_identical(ml$omega, mlr$omega)
  expect_identical(ml$errors, mlr$errors)
})

test_that("an improper tau-equivalent solution whose error variances sum below zero is refused", {
  fit <- list(converged = TRUE, improper = TRUE, no_solution = FALSE,
              errors = c(-2, 0.5, 0.5), omega = 1.4,
              items = c("y1", "y2", "y3"))
  expect_error(.omega_point_gate(fit, model = "tau-equivalent single-factor"),
               "The tau-equivalent single-factor solution is improper beyond use: the error variance of y1 \\(-2\\) is negative")
  ok <- fit
  ok$errors <- c(-0.1, 0.5, 0.5)
  expect_warning(.omega_point_gate(ok),
                 "Heywood case\\): the error variance of y1 is -0\\.1")
  # The bootstrap's view: a coefficient from an improper solution
  # carries the reason, and a solution with no coefficient gives none.
  expect_identical(attr(.omega_point_value(ok), "dmar_reason"), "improper")
  expect_identical(as.numeric(.omega_point_value(ok)), 1.4)
  none <- .omega_point_value(fit)
  expect_true(is.na(none))
  expect_identical(attr(none, "dmar_reason"), "improper")
})

test_that("the tau-equivalent bootstrap treats improper replications by the rule", {
  r <- reliability_alpha(estimator = "model_implied", data = items_continuous,
                         ci_method = "percentile", B = 40, seed = 113,
                         improper = "drop")
  bs <- attr(r, "bootstrap")
  expect_identical(bs$improper_rule, "drop")
  expect_identical(bs$improper + bs$improper_kept, 0L)
  ok <- list(converged = TRUE, improper = TRUE, no_solution = FALSE,
             errors = c(-0.1, 0.5, 0.5), omega = 0.8,
             items = c("y1", "y2", "y3"))
  expect_error(.omega_point_gate(ok, improper = "drop", bootstrap = TRUE),
               "improper = \"drop\"")
  expect_warning(.omega_point_gate(ok, improper = "drop", bootstrap = FALSE),
                 "Heywood")
})
