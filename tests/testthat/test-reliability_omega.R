# Tests for reliability_omega(). The point estimate is checked against
# the canonical (sum lambda)^2 / [(sum lambda)^2 + sum psi^2] formula
# applied to parameter estimates from a single-factor CFA fit by lavaan,
# and the various CI methods are checked for invariants.

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


test_that("reliability_omega() matches the closed-form (sum lambda)^2 / total formula", {
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied", ci_method = "none")
  est <- r$value[r$term == "estimate"]

  load_line <- paste(paste0("a", 1:J, "*y", 1:J), collapse = " + ")
  err_line  <- paste(paste0("y", 1:J, " ~~ b", 1:J, "*y", 1:J),
                     collapse = "\n")
  model <- paste0("f1 =~ NA*y1 + ", load_line, "\nf1 ~~ 1*f1\n", err_line)
  fit <- lavaan::cfa(model, data = items_continuous, estimator = "ML",
                     se = "none")
  pe <- lavaan::parameterEstimates(fit)
  var_names <- paste0("y", seq_len(J))
  lams <- pe[pe$op == "=~", "est"]
  psis <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% var_names, "est"]
  manual <- (sum(lams))^2 / ((sum(lams))^2 + sum(psis))

  expect_equal(est, manual, tolerance = 1e-6)
})

test_that("reliability_omega() default is robust omega, no bootstrap, with a message", {
  expect_message(
    r <- reliability_omega(data = items_continuous),
    regexp = "bootstrap based"
  )
  expect_s3_class(r, "data.frame")
  expect_named(r, c("term", "value"))
  expect_true(all(c("estimate", "se", "lower_limit", "upper_limit",
                    "conf_level", "N", "J") %in% r$term))
  expect_equal(attr(r, "coefficient"), "omega")
  expect_equal(attr(r, "denominator"), "observed")
  expect_equal(attr(r, "ci_method"), "none")
  expect_true(is.na(r$value[r$term == "lower_limit"]))
})

test_that("reliability_omega() model implied denominator defaults to the mlr interval", {
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied")
  expect_equal(attr(r, "ci_method"), "mlr")
  expect_equal(attr(r, "denominator"), "model_implied")
  expect_false(is.na(r$value[r$term == "lower_limit"]))
})

test_that("reliability_omega() ML/MLR Wald CIs bracket the estimate", {
  for (m in c("ml", "mlr", "ml_logistic", "mlr_logistic")) {
    r <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
  }
})

test_that("reliability_omega() closed-form alpha CIs apply mechanically to omega", {
  for (m in c("feldt", "fisher", "bonett", "hakstian_whalen")) {
    r <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
  }
})

test_that("reliability_omega() transformation intervals report both standard error scales", {
  # The se row is on the coefficient scale (the delta method
  # back-transform of the cube-root-scale standard error); the
  # transformation-scale value travels in the se_transformed row,
  # named by the se_transform_scale attribute (decided 2026-08-12).
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied",
                         ci_method = "hakstian_whalen")
  est  <- r$value[r$term == "estimate"]
  se   <- r$value[r$term == "se"]
  se_t <- r$value[r$term == "se_transformed"]
  expect_identical(attr(r, "se_transform_scale"), "cube_root")
  df_n <- N - 1
  df_d <- (N - 1) * (J - 1)
  a <- (1 - 2 / (9 * df_n)) / (1 - 2 / (9 * df_d))
  expect_equal(se, se_t * 3 * a^3 * (1 - est)^(2 / 3), tolerance = 1e-12)
  # The delta method (Wald) intervals report the model standard error on
  # the coefficient scale directly, with no transformation row.
  r_ml <- reliability_omega(data = items_continuous,
                            denominator = "model_implied",
                            ci_method = "ml")
  expect_false("se_transformed" %in% r_ml$term)
  expect_null(attr(r_ml, "se_transform_scale"))
})

test_that("reliability_omega() ML CI works from covariance matrix input", {
  r <- reliability_omega(S = S_continuous, N = N,
                         denominator = "model_implied", ci_method = "ml")
  est <- r$value[r$term == "estimate"]
  expect_lte(r$value[r$term == "lower_limit"], est)
  expect_gte(r$value[r$term == "upper_limit"], est)
})

test_that("reliability_omega() bootstrap CIs bracket the estimate", {
  # Skipped on CRAN: each bootstrap iteration refits the single-factor
  # CFA via lavaan, so even at B = 100 the three-method loop runs for
  # several seconds. The closed-form Wald / Feldt / Fisher / Bonett /
  # Hakstian-Whalen CIs above still cover the bracketing invariants.
  skip_on_cran()
  skip_if_not_installed("boot")
  for (m in c("bootstrap_se", "percentile", "bca")) {
    r <- reliability_omega(data = items_continuous,
                           denominator = "model_implied",
                           ci_method = m, B = 100)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
    expect_equal(attr(r, "B"), 100)
  }
})

test_that("reliability_omega() higher conf_level gives a wider CI (MLR)", {
  r80 <- reliability_omega(data = items_continuous,
                           denominator = "model_implied",
                           ci_method = "mlr", conf_level = 0.80)
  r95 <- reliability_omega(data = items_continuous,
                           denominator = "model_implied",
                           ci_method = "mlr", conf_level = 0.95)
  width80 <- r80$value[r80$term == "upper_limit"] -
             r80$value[r80$term == "lower_limit"]
  width95 <- r95$value[r95$term == "upper_limit"] -
             r95$value[r95$term == "lower_limit"]
  expect_gt(width95, width80)
})

test_that("reliability_omega() errors when MLR/ADF/bootstrap requested with S only", {
  for (m in c("mlr", "adf")) {
    expect_error(
      reliability_omega(S = S_continuous, N = N,
                        denominator = "model_implied", ci_method = m),
      regexp = "requires raw 'data'"
    )
  }
  for (m in c("percentile", "bca")) {
    expect_error(
      reliability_omega(S = S_continuous, N = N, ci_method = m),
      regexp = "requires raw 'data'"
    )
  }
})

test_that("reliability_omega() errors when neither data nor S is supplied", {
  expect_error(reliability_omega(), regexp = "Either 'data'.*or.*'S'")
})

test_that("reliability_omega() ci_method = 'none' returns NA bounds", {
  r <- reliability_omega(data = items_continuous, ci_method = "none")
  expect_identical(attr(r, "denominator"), "observed")
  expect_true(is.na(r$value[r$term == "lower_limit"]))
  expect_true(is.na(r$value[r$term == "upper_limit"]))
})



# ---------------------------------------------------------------------------
# denominator = "observed" (the coefficient Kelley & Pornprasertmanit, 2016,
# call hierarchical omega; formerly the separate reliability_omega_h())
# ---------------------------------------------------------------------------

test_that("observed denominator matches the N-divisor (sum lambda)^2 / total-variance formula", {
  r <- reliability_omega(data = items_continuous, denominator = "observed",
                         ci_method = "none")
  est <- r$value[r$term == "estimate"]
  expect_identical(attr(r, "denominator"), "observed")
  expect_identical(attr(r, "coefficient"), "omega")

  var_names <- paste0("y", seq_len(J))
  load_line <- paste(paste0("a", 1:J, "*y", 1:J), collapse = " + ")
  err_line  <- paste(paste0("y", 1:J, " ~~ b", 1:J, "*y", 1:J),
                     collapse = "\n")
  model <- paste0("f1 =~ NA*y1 + ", load_line, "\nf1 ~~ 1*f1\n", err_line)
  fit <- lavaan::cfa(model, data = items_continuous, estimator = "MLR",
                     se = "none", missing = "listwise")
  pe <- lavaan::parameterEstimates(fit)
  lams <- pe[pe$op == "=~", "est"]
  # Denominator on the loadings' maximum-likelihood (N-divisor) metric:
  # sum(cov) uses the N - 1 divisor, so rescale by (n - 1) / n.
  n_cc <- nrow(items_continuous)
  manual <- (sum(lams))^2 / (sum(S_continuous) * (n_cc - 1) / n_cc)
  expect_equal(est, manual, tolerance = 1e-6)

  # Anchor to the reference implementation on the same metric. Pinned from
  # MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  expect_equal(est, 0.7993963239470051, tolerance = 1e-6)
})

test_that("observed and model implied denominators agree closely under good fit", {
  ro <- reliability_omega(data = items_continuous, denominator = "observed",
                          ci_method = "none")
  rm_ <- reliability_omega(data = items_continuous, ci_method = "none")
  expect_equal(ro$value[ro$term == "estimate"],
               rm_$value[rm_$term == "estimate"], tolerance = 0.05)
})

test_that("observed denominator supports S-only input for the point estimate", {
  r <- reliability_omega(S = S_continuous, N = N, denominator = "observed",
                         ci_method = "none")
  est <- r$value[r$term == "estimate"]
  expect_true(!is.na(est) && est > 0 && est < 1)
})

test_that("robust omega never bootstraps by default and messages the request path", {
  expect_message(
    r_S <- reliability_omega(S = S_continuous, N = N,
                             denominator = "observed"),
    regexp = "covariance matrix"
  )
  expect_identical(attr(r_S, "ci_method"), "none")

  expect_message(
    r_d <- reliability_omega(data = items_continuous,
                             denominator = "observed"),
    regexp = "ci_method"
  )
  expect_identical(attr(r_d, "ci_method"), "none")
  expect_true(is.na(r_d$value[r_d$term == "lower_limit"]))
})

test_that("observed denominator rejects closed-form CI methods", {
  for (m in c("mlr", "ml", "adf", "feldt", "fisher", "bonett",
              "hakstian_whalen")) {
    expect_error(
      reliability_omega(data = items_continuous, denominator = "observed",
                        ci_method = m),
      regexp = "model implied total variance"
    )
  }
})

test_that("observed denominator errors when a bootstrap CI is requested with S only", {
  expect_error(
    reliability_omega(S = S_continuous, N = N, denominator = "observed",
                      ci_method = "percentile"),
    regexp = "requires raw 'data'"
  )
})

test_that("observed denominator bootstrap CIs bracket the estimate and are reproducible", {
  skip_on_cran()
  skip_if_not_installed("boot")
  for (m in c("percentile", "bca", "bootstrap_se", "bootstrap_se_logistic")) {
    r <- reliability_omega(data = items_continuous, denominator = "observed",
                           ci_method = m, B = 200, seed = 113)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
  }
  r1 <- reliability_omega(data = items_continuous, denominator = "observed",
                          ci_method = "percentile", B = 100, seed = 113)
  r2 <- reliability_omega(data = items_continuous, denominator = "observed",
                          ci_method = "percentile", B = 100, seed = 113)
  expect_equal(r1$value, r2$value)
  expect_equal(attr(r1, "B"), 100)
})

test_that("reliability_omega() default denominator is observed (robust omega)", {
  r <- reliability_omega(data = items_continuous, ci_method = "none")
  expect_identical(attr(r, "denominator"), "observed")
})


test_that("likelihood CI brackets the estimate and matches MBESS interval.type = 'll'", {
  skip_on_cran()  # two likelihood profiles; the closed-form anchors above run on CRAN
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied",
                         ci_method = "likelihood")
  est <- r$value[r$term == "estimate"]
  lo  <- r$value[r$term == "lower_limit"]
  hi  <- r$value[r$term == "upper_limit"]
  expect_true(lo < est && est < hi)
  expect_gte(lo, 0)
  expect_lte(hi, 1)
  expect_true(is.na(r$value[r$term == "se"]))

  # Covariance input reproduces raw-data input; both are maximum
  # likelihood on the same second moments.
  r_S <- reliability_omega(S = S_continuous, N = N,
                           denominator = "model_implied",
                           ci_method = "likelihood")
  expect_equal(r_S$value[r_S$term == "lower_limit"], lo, tolerance = 1e-5)
  expect_equal(r_S$value[r_S$term == "upper_limit"], hi, tolerance = 1e-5)

  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(lo, 0.7581229465880676, tolerance = 5e-4)
  expect_equal(hi, 0.8351069201742349, tolerance = 5e-4)
})

test_that("likelihood CI is rejected for the observed denominator", {
  expect_error(
    reliability_omega(data = items_continuous, ci_method = "likelihood"),
    regexp = "model implied total variance"
  )
})

test_that("a two-item scale is refused: the congeneric model is not identified", {
  # With two items the model has four free parameters against three
  # observed moments, so any returned omega is solver-start noise
  # (three lavaan starts gave 0.654, 0.652, 0.662 on the same data).
  set.seed(11)
  d2 <- as.data.frame(MASS::mvrnorm(
    n = 120, mu = c(0, 0),
    Sigma = matrix(c(1, 0.49, 0.49, 1), 2, 2)
  ))
  names(d2) <- c("x1", "x2")
  expect_error(reliability_omega(data = d2),
               regexp = "three or more items")
  expect_error(reliability_omega(data = d2, denominator = "model_implied",
                                 ci_method = "ml"),
               regexp = "three or more items")
  expect_error(reliability_omega(S = cov(d2), N = 120,
                                 denominator = "model_implied",
                                 ci_method = "ml"),
               regexp = "three or more items")
  # The two-item tau-equivalent model is just identified and stays legal.
  expect_no_error(reliability_alpha(data = d2, estimator = "model_implied",
                                    ci_method = "ml"))
})

test_that("an improper solution is diagnosed as such, not as nonconvergence", {
  skip_if_not_installed("lavaan")
  # A near-unit loading in a small sample gives a converged fit with a
  # negative error variance. That is a Heywood case, not a failure to
  # converge, and calling it nonconvergence sent the user looking for the
  # wrong remedy. cfa_k() already warned and continued; the single-factor
  # path errored. They now agree.
  set.seed(2)
  J <- 3; n <- 35
  lam <- c(.95, .85, .35)
  eta <- rnorm(n)
  d <- as.data.frame(outer(eta, lam) +
        matrix(rnorm(n * J), n, J) %*% diag(sqrt(1 - lam^2)))
  names(d) <- paste0("y", seq_len(J))

  # Collect every warning rather than matching one, since lavaan emits its
  # own notice about the negative variance alongside DMAR's.
  msgs <- character(0)
  r <- withCallingHandlers(
    reliability_omega(data = d, denominator = "model_implied",
                      ci_method = "ml"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_true(any(grepl("Heywood", msgs)))
  expect_false(any(grepl("did not converge", msgs)))
  # The estimate is returned rather than refused, matching how MBESS and
  # cfa_k() treat an improper but converged solution.
  expect_true(is.finite(r$value[r$term == "estimate"]))
})
