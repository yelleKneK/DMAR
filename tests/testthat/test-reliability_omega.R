# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

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
  expect_true(all(c("estimate", "se", "ci_lower", "ci_upper",
                    "conf_level", "N", "n_items") %in% r$term))
  expect_equal(attr(r, "coefficient"), "omega")
  expect_equal(attr(r, "denominator"), "observed")
  expect_equal(attr(r, "ci_method"), "none")
  expect_true(is.na(r$value[r$term == "ci_lower"]))
})

test_that("reliability_omega() model implied denominator defaults to the mlr_logistic interval", {
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied")
  expect_equal(attr(r, "ci_method"), "mlr_logistic")
  expect_equal(attr(r, "denominator"), "model_implied")
  expect_false(is.na(r$value[r$term == "ci_lower"]))
  expect_identical(r$value,
                   reliability_omega(data = items_continuous,
                                     denominator = "model_implied",
                                     ci_method = "mlr_logistic")$value)
  # The same robust standard error as "mlr", the DMAR 1.0.x default, with
  # the limits formed on the logit scale and transformed back.
  old <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = "mlr")
  v <- function(x, t) x$value[x$term == t]
  expect_equal(v(r, "estimate"), v(old, "estimate"))
  expect_equal(v(r, "se"), v(old, "se"))
  est <- v(r, "estimate")
  half <- stats::qnorm(0.975) * v(r, "se") / (est * (1 - est))
  expect_equal(c(v(r, "ci_lower"), v(r, "ci_upper")),
               stats::plogis(stats::qlogis(est) + c(-1, 1) * half),
               tolerance = 1e-12)
  expect_equal(c(v(old, "ci_lower"), v(old, "ci_upper")),
               est + c(-1, 1) * stats::qnorm(0.975) * v(r, "se"),
               tolerance = 1e-12)
  # Under every treatment of missing values the default reaches.
  d_na <- as.data.frame(items_continuous)
  d_na$y2[seq(1, N, by = 7)] <- NA
  for (m in c("listwise", "fiml", "two_stage")) {
    x <- suppressWarnings(reliability_omega(d_na, missing = m,
                                            denominator = "model_implied"))
    expect_identical(attr(x, "ci_method"), "mlr_logistic", label = m)
  }
})

test_that("omega's logit-scale default gives no warning, and the paths it does not reach keep their default", {
  d <- items_continuous
  mi <- function(...) reliability_omega(d, denominator = "model_implied", ...)
  # An explicit ci_method, robust omega's default (no interval), and a
  # call refused for its data.
  expect_no_warning(mi(ci_method = "mlr"))
  expect_no_warning(mi(ci_method = "mlr_logistic"))
  expect_no_warning(suppressMessages(reliability_omega(d)))
  expect_error(reliability_omega(d[, 1:2], denominator = "model_implied"),
               "three or more items", fixed = TRUE)
  # Under wald_df = "missing_information", calibrated on the scale of omega
  # only, the default stays on that scale, "mlr" ("ml" with a covariance
  # matrix).
  k <- mi(wald_df = "missing_information")
  expect_identical(attr(k, "ci_method"), "mlr")
  expect_identical(k$value, mi(ci_method = "mlr",
                               wald_df = "missing_information")$value)
  ks <- reliability_omega(S = S_continuous, N = N,
                          denominator = "model_implied",
                          wald_df = "missing_information")
  expect_identical(attr(ks, "ci_method"), "ml")
  expect_identical(ks$value,
                   reliability_omega(S = S_continuous, N = N,
                                     denominator = "model_implied",
                                     ci_method = "ml",
                                     wald_df = "missing_information")$value)
  # The default interval gives no warning, and the same limits each time.
  expect_no_warning(r <- mi())
  expect_no_warning(r2 <- mi())
  expect_identical(r2$value, r$value)
  # An explicit ci_method = NULL is the default, as match.arg() made it
  # through DMAR 1.0.x and as reliability() uses it: the choice under the
  # keyword, and robust omega's point estimate with its advice.
  expect_no_warning(n1 <- mi(ci_method = NULL))
  expect_identical(attr(n1, "ci_method"), "mlr_logistic")
  expect_identical(n1$value, r$value)
  n2 <- mi(ci_method = NULL, wald_df = "missing_information")
  expect_identical(attr(n2, "ci_method"), "mlr")
  expect_identical(n2$value, k$value)
  expect_message(n3 <- reliability_omega(d, ci_method = NULL),
                 "bootstrap based", fixed = TRUE)
  expect_identical(attr(n3, "ci_method"), "none")
  # A covariance matrix, which the robust standard error cannot use, takes
  # "ml_logistic" without ci_method, as reliability_H() does (through DMAR
  # 1.0.x such a call stopped).
  expect_no_warning(s1 <- reliability_omega(S = S_continuous, N = N,
                                            denominator = "model_implied"))
  expect_identical(attr(s1, "ci_method"), "ml_logistic")
  expect_identical(s1$value,
                   reliability_omega(S = S_continuous, N = N,
                                     denominator = "model_implied",
                                     ci_method = "ml_logistic")$value)
  s2 <- reliability_omega(S = S_continuous, N = N,
                          denominator = "model_implied", ci_method = NULL)
  expect_identical(s2$value, s1$value)
  # reliability() follows, leaving the choice to reliability_omega(): with
  # raw data "mlr_logistic", and with a covariance matrix "ml_logistic"
  # (the wrapper chose "ml" there through DMAR 1.0.x).
  expect_no_warning(w1 <- reliability(d, type = "omega",
                                      denominator = "model_implied"))
  expect_identical(attr(w1, "ci_method"), "mlr_logistic")
  expect_identical(w1$value, r$value)
  expect_no_warning(w2 <- reliability(S = S_continuous, N = N, type = "omega",
                                      denominator = "model_implied"))
  expect_identical(attr(w2, "ci_method"), "ml_logistic")
  expect_identical(w2$value, s1$value)
  expect_no_warning(w3 <- reliability(S = S_continuous, N = N,
                                      type = "omega",
                                      denominator = "model_implied",
                                      wald_df = "missing_information"))
  expect_identical(attr(w3, "ci_method"), "ml")
  expect_identical(w3$value, ks$value)
  # The model implied alpha keeps its default.
  expect_no_warning(a <- reliability(S = S_continuous, N = N, type = "alpha",
                                     estimator = "model_implied"))
  expect_identical(attr(a, "ci_method"), "ml")
})

test_that("reliability_omega() ML/MLR Wald CIs bracket the estimate", {
  for (m in c("ml", "mlr", "ml_logistic", "mlr_logistic")) {
    r <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
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
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
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
  expect_lte(r$value[r$term == "ci_lower"], est)
  expect_gte(r$value[r$term == "ci_upper"], est)
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
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
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
  width80 <- r80$value[r80$term == "ci_upper"] -
             r80$value[r80$term == "ci_lower"]
  width95 <- r95$value[r95$term == "ci_upper"] -
             r95$value[r95$term == "ci_lower"]
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
  expect_true(is.na(r$value[r$term == "ci_lower"]))
  expect_true(is.na(r$value[r$term == "ci_upper"]))
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
  expect_true(is.na(r_d$value[r_d$term == "ci_lower"]))
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
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
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
  lo  <- r$value[r$term == "ci_lower"]
  hi  <- r$value[r$term == "ci_upper"]
  expect_true(lo < est && est < hi)
  expect_gte(lo, 0)
  expect_lte(hi, 1)
  expect_true(is.na(r$value[r$term == "se"]))

  # Covariance input reproduces raw-data input; both are maximum
  # likelihood on the same second moments.
  r_S <- reliability_omega(S = S_continuous, N = N,
                           denominator = "model_implied",
                           ci_method = "likelihood")
  expect_equal(r_S$value[r_S$term == "ci_lower"], lo, tolerance = 1e-5)
  expect_equal(r_S$value[r_S$term == "ci_upper"], hi, tolerance = 1e-5)

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
  # wrong remedy. fa_confirmatory() already warned and continued; the single-factor
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
  # fa_confirmatory() treat an improper but converged solution.
  expect_true(is.finite(r$value[r$term == "estimate"]))
})


test_that("a ci_method that does not apply to the denominator stops with a message naming the one that does", {
  # The help page states which intervals apply to which denominator;
  # this pins the code to that statement. Every closed form, Wald
  # interval, and the profile likelihood is derived under the model
  # implied ratio, so each stops under the default (observed)
  # denominator with a message that names denominator = "model_implied".
  model_implied_only <- c("ml", "mlr", "ml_logistic", "mlr_logistic",
                          "likelihood", "adf", "adf_logistic",
                          "feldt", "fisher", "bonett", "hakstian_whalen")
  for (m in model_implied_only) {
    expect_error(
      reliability_omega(data = items_continuous, ci_method = m),
      regexp = 'denominator = "model_implied"', fixed = TRUE,
      label = paste("ci_method =", m, "under the observed denominator")
    )
  }
  # "none" and the bootstrap methods apply to both denominators; the
  # bootstrap runs are covered by the bracketing tests above, so only
  # "none" is exercised here.
  for (d in c("observed", "model_implied")) {
    expect_no_error(reliability_omega(data = items_continuous,
                                      ci_method = "none", denominator = d))
  }
  # The menu on the help page is the menu in the signature.
  expect_setequal(eval(formals(reliability_omega)$ci_method),
                  c(model_implied_only, "bootstrap_se",
                    "bootstrap_se_logistic", "percentile", "bca", "none"))
})


test_that("reliability_omega() carries the residual correlations of the single-factor fit", {
  # A shared specific factor on y1 and y2 inflates their loadings under
  # either denominator; the residual correlations are where it shows.
  # (The seed differs from the file's 113 so that the doublet is not a
  # copy of the common factor drawn at the top of the file.)
  set.seed(2026)
  doublet <- rnorm(N)
  items_d <- items_continuous
  items_d[, 1:2] <- items_d[, 1:2] + 0.6 * doublet
  for (den in c("observed", "model_implied")) {
    r <- suppressMessages(reliability_omega(data = items_d, denominator = den))
    rc <- attr(r, "residual_cor")
    expect_type(rc, "list")
    expect_setequal(rc$items, c("y1", "y2"))
    expect_gt(rc$largest, 0)
  }
  # The matrix is lavaan's correlation residuals of the maximum
  # likelihood single-factor fit to the same cases.
  r <- suppressMessages(reliability_omega(data = items_d))
  fit <- lavaan::cfa(paste("f =~", paste(colnames(items_d), collapse = " + ")),
                     data = as.data.frame(items_d), std.lv = TRUE)
  ref <- unclass(lavaan::residuals(fit, type = "cor")$cov)
  expect_equal(attr(r, "residual_cor")$matrix,
               matrix(as.numeric(ref), J, J, dimnames = dimnames(ref)),
               tolerance = 1e-6)
  # Covariance input carries it too.
  r_S <- suppressMessages(reliability_omega(S = cov(items_d), N = N))
  expect_setequal(attr(r_S, "residual_cor")$items, c("y1", "y2"))
})

test_that("the printed table names the denominator and the largest residual correlation", {
  r_obs <- suppressMessages(reliability_omega(data = items_continuous))
  r_mod <- reliability_omega(data = items_continuous,
                             denominator = "model_implied")
  expect_output(print(r_obs), "Denominator: observed composite variance \\(robust omega\\)")
  expect_output(print(r_mod), "Denominator: composite variance implied by the fitted model")
  expect_output(print(r_obs), "Largest residual correlation under the single-factor model")
})

test_that("the residual correlations come from a maximum likelihood fit whatever the interval method", {
  # The 1.1.0 review: with ci_method = "adf" the model implied branch
  # took the residuals from its WLS fit, whose residuals differ from the
  # ML ones and here fall on another pair (y1 with y6) than the planted
  # doublet (y1 with y2).
  set.seed(113)
  n <- 300
  eta <- rnorm(n); s <- rnorm(n)
  Y <- outer(eta, rep(.6, 6)) + (matrix(rexp(n * 6), n, 6) - 1) * 0.7
  Y[, 1:2] <- Y[, 1:2] + 0.5 * s
  colnames(Y) <- paste0("y", 1:6)
  mod <- paste("f =~", paste(colnames(Y), collapse = " + "))
  ref <- function(fit) {
    R <- lavaan::residuals(fit, type = "cor")$cov
    matrix(as.numeric(R), 6, 6, dimnames = dimnames(R))
  }
  ml <- ref(lavaan::cfa(mod, data = as.data.frame(Y), std.lv = TRUE))
  wls <- ref(lavaan::cfa(mod, data = as.data.frame(Y), std.lv = TRUE,
                         estimator = "WLS"))
  expect_gt(max(abs(wls - ml)), 0.1)
  for (m in c("adf", "adf_logistic", "ml", "mlr")) {
    r <- reliability_omega(data = Y, denominator = "model_implied",
                           ci_method = m)
    rc <- attr(r, "residual_cor")
    expect_equal(rc$matrix, ml, tolerance = 1e-6, label = m)
    expect_setequal(rc$items, c("y1", "y2"))
  }
})

test_that("the residual correlations cover the items only, under their own names, with auxiliary variables", {
  # With auxiliary variables under FIML the single-factor fit names the
  # items y1, ..., yJ and the auxiliaries a1, ..., aK, so a lookup by the
  # user's item names matched nothing: the summary kept the auxiliary's
  # row and named its pair by position (here "y1 with y2", beside an a1
  # row). Sentence scores are made missing more often for pupils with
  # high word classification scores, the auxiliary.
  vars <- c("t5_general_information", "t6_paragraph_comprehension",
            "t7_sentence", "t9_word_meaning")
  aux <- "t8_word_classification"
  d <- holzinger_swineford[, c(vars, aux)]
  set.seed(113)
  high <- d[[aux]] > stats::median(d[[aux]])
  d$t7_sentence[stats::runif(nrow(d)) < ifelse(high, 0.40, 0.10)] <- NA
  resid_items <- function(fit) {
    R <- lavaan::residuals(fit, type = "cor")$cov[vars, vars]
    matrix(as.numeric(R), 4, 4, dimnames = list(vars, vars))
  }
  model <- paste("f =~", paste(vars, collapse = " + "))
  # FIML with the auxiliary as a saturated correlate (Graham, 2003): free
  # covariances with every item residual and a free variance. The
  # residual covariance matrix of that parameterization is not positive
  # definite when the auxiliary shares variance with the factor, which
  # lavaan's check after estimation would announce, as the package's own
  # fit knows (.omega_fit_cfa_aux()).
  fiml_ref <- resid_items(lavaan::cfa(
    paste0(model, "\n", aux, " ~~ ", aux, "\n", aux, " ~~ ",
           paste(vars, collapse = " + ")),
    data = d, std.lv = TRUE, missing = "ml", check.post = FALSE))
  ts_ref <- resid_items(lavaan::cfa(model, data = d, std.lv = TRUE,
                                    missing = "two.stage", aux = aux))
  for (m in c("fiml", "two_stage")) {
    ref <- if (m == "fiml") fiml_ref else ts_ref
    for (den in c("observed", "model_implied")) {
      lab <- paste(m, den)
      r <- suppressMessages(reliability_omega(data = d, aux = aux,
                                              missing = m,
                                              denominator = den))
      rc <- attr(r, "residual_cor")
      expect_identical(dimnames(rc$matrix), list(vars, vars), label = lab)
      expect_equal(rc$matrix, ref, tolerance = 1e-5, label = lab)
      k <- which(abs(ref) == max(abs(ref[upper.tri(ref)])) &
                   upper.tri(ref), arr.ind = TRUE)
      expect_identical(rc$items, vars[k[1L, ]], label = lab)
      expect_output(print(r), paste0("(", vars[k[1L, 1L]], " with ",
                                     vars[k[1L, 2L]], ")"),
                    fixed = TRUE)
    }
  }
  # Listwise deletion names the items as it always did.
  r <- reliability_omega(data = d[vars], denominator = "model_implied",
                         missing = "listwise")
  expect_identical(rownames(attr(r, "residual_cor")$matrix), vars)
  # A fit that does not hold every item gives no summary rather than one
  # over other variables.
  fit <- .omega_fit_cfa(data = d, aux = aux, estimator = "ML", se = "none")
  expect_identical(fit$lavaan_items, paste0("y", 1:4))
  expect_null(.omega_residual_cor(fit$fit, vars))
  expect_identical(rownames(.omega_residual_cor(fit$fit, vars,
                                                fit$lavaan_items)$matrix),
                   vars)
})

test_that("with three items the saturated model leaves no residual to name", {
  # A proper three-item solution reproduces every correlation, so the
  # residuals are optimizer noise and the pair carrying the largest of
  # them is arbitrary; the summary says the model is saturated instead.
  for (den in c("observed", "model_implied")) {
    r3 <- suppressMessages(reliability_omega(data = items_continuous[, 1:3],
                                             denominator = den))
    rc <- attr(r3, "residual_cor")
    expect_true(rc$saturated)
    expect_true(is.na(rc$largest))
    expect_null(rc$items)
    expect_output(print(r3), "With three items the single-factor model is saturated (df = 0)",
                  fixed = TRUE)
  }
  r4 <- suppressMessages(reliability_omega(data = items_continuous[, 1:4]))
  expect_false(attr(r4, "residual_cor")$saturated)
  expect_true(is.finite(attr(r4, "residual_cor")$largest))
})

test_that("a three-item fit that cannot reproduce the correlations still reports its residual", {
  # With df = 0 a residual above the tolerance means no real loadings fit
  # the correlations (their product is negative), which is the one case
  # where the residual is the signal, so it is not suppressed.
  R <- matrix(c(0, 0.08, -0.06, 0.08, 0, 1e-9, -0.06, 1e-9, 0), 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  rc <- .residual_cor_summary(R, df = 0)
  expect_false(rc$saturated)
  expect_equal(rc$largest, 0.08)
  expect_identical(rc$items, c("a", "b"))
  # Zero within the tolerance with df = 0 is saturated; without df the
  # check is not made.
  Z <- R * 1e-8
  expect_true(.residual_cor_summary(Z, df = 0)$saturated)
  expect_false(.residual_cor_summary(Z)$saturated)
})


# ---- One maximum likelihood fit, and improper solutions --------------------

hs_items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
              "t4_lozenges")

test_that("the bootstrap refits by ML, whose estimates equal MLR's to the last bit", {
  # MLR changes only the standard errors and the test statistic, so the
  # replications refit by ML (since 1.1.0) with exactly the estimates
  # the MLR refits of 1.0.1 gave, complete data and FIML alike.
  X <- holzinger_swineford[, hs_items]
  set.seed(113)
  for (b in 1:4) {
    d <- X[sample(nrow(X), replace = TRUE), ]
    expect_identical(.omega_observed_from_data(d, "ML"),
                     .omega_observed_from_data(d, "MLR"))
    ml  <- .omega_fit_cfa(data = d, estimator = "ML", se = "none")
    mlr <- .omega_fit_cfa(data = d, estimator = "MLR", se = "none")
    expect_identical(ml$loadings, mlr$loadings)
    expect_identical(ml$errors, mlr$errors)
    expect_identical(ml$omega, mlr$omega)
  }
  Xm <- holzinger_swineford[, c(hs_items, "t25_paper_form_board_r",
                                "t26_flags")]
  for (b in 1:2) {
    d <- Xm[sample(nrow(Xm), replace = TRUE), ]
    ml  <- .omega_fit_cfa(data = d, estimator = "ML", se = "none",
                          missing = "ml")
    mlr <- .omega_fit_cfa(data = d, estimator = "MLR", se = "none",
                          missing = "ml")
    expect_identical(ml$loadings, mlr$loadings)
    expect_identical(ml$omega, mlr$omega)
  }
})

test_that("robust omega fits the single-factor model once per call", {
  # The estimate and the residual correlations come from one ML fit;
  # through 1.0.1 an MLR fit gave the estimate and a second, ML, fit the
  # residual correlations.
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .omega_fit_cfa
  testthat::local_mocked_bindings(
    .omega_fit_cfa = function(...) {
      calls$n <- calls$n + 1L
      orig(...)
    })
  r <- suppressMessages(reliability_omega(data = items_continuous))
  expect_identical(calls$n, 1L)
  expect_false(is.null(attr(r, "residual_cor")))
  # The residual correlations are those of an ML fit to the same data.
  fit <- .omega_fit_cfa(data = items_continuous, estimator = "ML",
                        se = "none")$fit
  expect_equal(attr(r, "residual_cor"), .omega_residual_cor(fit),
               tolerance = 1e-12)
})

# The sample reported on 2026-11-01: three items whose correlations have
# a negative product (population r12 = .4, r13 = .3, r23 = -.05, N = 300,
# seed 3). Through 1.0.1 this returned an omega of 471, from a loading of
# 43.6 and an error variance of -1899, with only a Heywood warning.
set.seed(3)
S_neg <- matrix(c(1, .4, .3, .4, 1, -.05, .3, -.05, 1), 3)
triad_neg <- as.data.frame(matrix(rnorm(300 * 3), 300) %*% chol(S_neg))
names(triad_neg) <- paste0("y", 1:3)

test_that("three items with no single-factor solution are refused under every route", {
  msg <- "no solution for these three items: the product of their correlations is negative"
  expect_error(suppressMessages(reliability_omega(data = triad_neg)), msg)
  expect_error(reliability_omega(data = triad_neg,
                                 denominator = "model_implied",
                                 ci_method = "ml"), msg)
  expect_error(reliability_omega(S = cov(triad_neg), N = 300,
                                 denominator = "model_implied",
                                 ci_method = "ml"), msg)
  expect_error(reliability_omega(data = triad_neg, ci_method = "percentile",
                                 B = 20, seed = 113), msg)
  # The message quotes the correlations.
  expect_error(suppressMessages(reliability_omega(data = triad_neg)),
               "r\\(y2, y3\\) = -0\\.076")
  # Under FIML the fit's saturated covariance matrix shows the same.
  d <- triad_neg
  d$y2[1:30] <- NA
  expect_error(suppressWarnings(
    reliability_omega(data = d, missing = "fiml",
                      denominator = "model_implied", ci_method = "ml")),
    msg)
  # The verdict, and the bootstrap's view of it: no coefficient.
  f <- .omega_fit_cfa(data = triad_neg, se = "none")
  expect_identical(.omega_fit_verdict(f), "no_solution")
  expect_null(f$fit)
  v <- .omega_observed_from_data(triad_neg)
  expect_true(is.na(v))
  # Named as such, so the bootstrap counts it apart from a Heywood case.
  expect_identical(attr(v, "dmar_reason"), "no_solution")
})

test_that("an improper solution whose error variances sum below zero is refused", {
  # Four weakly related items, two of which share a specific factor (a
  # doublet, r = .68 in this sample of 30): the single factor follows the
  # doublet, the error variance of y2 runs to -3.56 while the others stay
  # near one, and the error variances sum to -1.16, so the model implies
  # an omega of 1.24 (McDonald's, 2004, Type III).
  set.seed(1)
  R4 <- matrix(0.05, 4, 4)
  diag(R4) <- 1
  R4[1, 2] <- R4[2, 1] <- 0.6
  Y <- as.data.frame(matrix(rnorm(30 * 4), 30) %*% chol(R4))
  names(Y) <- paste0("y", 1:4)
  # lavaan's own notices of the improper solution do not leak.
  expect_no_warning(f <- .omega_fit_cfa(data = Y, se = "none"))
  expect_identical(.omega_fit_verdict(f), "negative_error")
  expect_lt(sum(f$errors), 0)
  expect_gt(f$omega, 1)
  for (den in c("observed", "model_implied")) {
    expect_error(
      suppressMessages(reliability_omega(data = Y, denominator = den,
                                         ci_method = "none")),
      "improper beyond use: the error variance of y2 \\(-3\\.557\\) is negative, and the error variances sum to -1\\.159")
  }
})

# The Heywood sample of the test above: an improper solution whose
# coefficient is still a proportion.
set.seed(2)
lam_hey <- c(.95, .85, .35)
eta_hey <- rnorm(35)
heywood3 <- as.data.frame(outer(eta_hey, lam_hey) +
  matrix(rnorm(35 * 3), 35, 3) %*% diag(sqrt(1 - lam_hey^2)))
names(heywood3) <- paste0("y", 1:3)

test_that("an improper solution that is still a proportion is returned with a warning naming the item, under either denominator", {
  expect_identical(.omega_fit_verdict(.omega_fit_cfa(data = heywood3,
                                                     se = "none")),
                   "improper")
  expect_warning(
    reliability_omega(data = heywood3, denominator = "model_implied",
                      ci_method = "ml"),
    "Heywood case\\): the error variance of y[1-3] is -[0-9.]+ \\(its observed variance is [0-9.]+\\)")
  # The observed denominator now says so too (through 1.0.1 only
  # lavaan's own notice did).
  expect_warning(suppressMessages(reliability_omega(data = heywood3)),
                 "Heywood")
})

test_that("improper = \"drop\" refuses an improper point estimate only when a bootstrap is requested", {
  expect_error(
    reliability_omega(data = heywood3, denominator = "model_implied",
                      ci_method = "percentile", B = 20, seed = 113,
                      improper = "drop"),
    "improper = \"drop\" sets aside every bootstrap replication")
  expect_warning(
    reliability_omega(data = heywood3, denominator = "model_implied",
                      ci_method = "ml", improper = "drop"),
    "Heywood")
})

test_that("a bootstrap counts improper replications, keeps them by default, and sets them aside under drop", {
  skip_on_cran()
  # A proper sample with one loading of .95 at N = 60, so that some
  # resamples push the loading past one.
  set.seed(4)
  lh <- c(.95, .7, .6, .5, .4)
  f <- rnorm(60)
  Y <- as.data.frame(sapply(lh, function(l) l * f + sqrt(1 - l^2) * rnorm(60)))
  names(Y) <- paste0("y", 1:5)
  expect_identical(.omega_fit_verdict(.omega_fit_cfa(data = Y, se = "none")),
                   "proper")
  w_keep <- testthat::capture_warnings(
    k <- reliability_omega(Y, denominator = "model_implied",
                           ci_method = "percentile", B = 40, seed = 113))
  w_drop <- testthat::capture_warnings(
    d <- reliability_omega(Y, denominator = "model_implied",
                           ci_method = "percentile", B = 40, seed = 113,
                           improper = "drop"))
  bk <- attr(k, "bootstrap")
  bd <- attr(d, "bootstrap")
  expect_gt(bk$improper_kept, 0L)
  expect_identical(bd$improper, bk$improper_kept + bk$improper)
  expect_identical(bd$improper_kept, 0L)
  expect_identical(attr(d, "B_used"), 40L - bd$failed)
  expect_identical(bk$t, bd$t)
  # One improper replication in 40 already reaches the 2.5 percent
  # beyond a limit.
  expect_true(bk$caution)
  expect_match(w_keep, "Caution:", fixed = TRUE)
  expect_output(print(k), "came from an improper \\(Heywood\\) solution and (was|were) kept")
  expect_output(print(k), "caution: [0-9.]+% of the replications came from an improper solution")
  # lavaan's notice of each improper replication does not leak.
  expect_false(any(grepl("ov variances are negative", c(w_keep, w_drop))))
})

test_that("a bootstrap of three items sets aside the resamples with no single-factor solution", {
  skip_on_cran()
  # Three weakly related items: about a third of the resamples have
  # correlations with a negative product. They are set aside and counted
  # as improper under either rule, and no fit is attempted for them.
  set.seed(8)
  l3 <- c(.45, .4, .35)
  f <- rnorm(60)
  Y <- as.data.frame(sapply(l3, function(l) l * f + sqrt(1 - l^2) * rnorm(60)))
  names(Y) <- paste0("y", 1:3)
  w <- testthat::capture_warnings(
    r <- reliability_omega(Y, ci_method = "percentile", B = 40, seed = 113))
  bs <- attr(r, "bootstrap")
  expect_gt(bs$improper, 0L)
  expect_identical(bs$improper, bs$failed)
  expect_true(all(is.na(bs$t[bs$status == "improper"])))
  expect_lte(max(bs$t, na.rm = TRUE), 1)
  # The resamples with no single-factor solution are counted on their own
  # (a subset of the improper ones, reported under that status), named
  # apart from the Heywood cases, and the caution says which way the
  # lower limit is likely wrong.
  expect_gt(bs$no_solution, 0L)
  expect_lte(bs$no_solution, bs$improper)
  expect_identical(levels(bs$status), DMAR:::.BOOT_STATUS)
  expect_match(w, sprintf("%d had no single-factor solution",
                          bs$no_solution), fixed = TRUE, all = FALSE)
  expect_match(w, "belongs on neither side", fixed = TRUE, all = FALSE)
  expect_match(w, "the lower limit is likely too high", fixed = TRUE,
               all = FALSE)
  expect_output(print(r), "had no single-factor solution", fixed = TRUE)
  expect_output(print(r), "the lower limit is likely too high", fixed = TRUE)
  # (boot also warns that so few returned replications put the
  # percentile limits at extreme order statistics.)
  if (bs$improper > bs$no_solution) {
    expect_match(w, "gave an improper (Heywood) solution", fixed = TRUE,
                 all = FALSE)
    expect_match(w, "an error variance fell below zero", fixed = TRUE,
                 all = FALSE)
  }
  expect_false(any(grepl("loading ran past one", w, fixed = TRUE)))
})


# ---- wald_df: a t reference distribution for the Wald intervals ------------

test_that("wald_df = Inf reproduces every Wald interval of reliability_omega() exactly", {
  # The default is the normal reference; qt(p, Inf) is qnorm(p) to the last
  # bit, so the whole result, values and attributes, is unchanged whether
  # or not wald_df is named.
  for (m in c("ml", "mlr", "ml_logistic", "mlr_logistic", "adf",
              "adf_logistic", "bonett")) {
    a <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m)
    b <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m,
                           wald_df = Inf)
    expect_identical(a, b)
    expect_identical(attr(a, "wald_df"), Inf)
  }
  d <- as.data.frame(items_continuous)
  d$y2[seq(1, N, by = 4)] <- NA
  for (m in c("ml", "mlr_logistic")) {
    expect_identical(
      reliability_omega(data = d, denominator = "model_implied",
                        ci_method = m, missing = "fiml"),
      reliability_omega(data = d, denominator = "model_implied",
                        ci_method = m, missing = "fiml", wald_df = Inf))
  }
  expect_identical(
    reliability_omega(S = S_continuous, N = N,
                      denominator = "model_implied", ci_method = "ml"),
    reliability_omega(S = S_continuous, N = N,
                      denominator = "model_implied", ci_method = "ml",
                      wald_df = Inf))
  # The interval itself is the normal one at the default.
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied", ci_method = "ml")
  est <- r$value[r$term == "estimate"]
  se <- r$value[r$term == "se"]
  expect_identical(r$value[r$term == "ci_lower"],
                   est - stats::qnorm(0.975) * se)
})

test_that("a finite wald_df puts the t quantile into the Wald interval", {
  df <- 12
  r0 <- reliability_omega(data = items_continuous,
                          denominator = "model_implied", ci_method = "ml")
  r <- reliability_omega(data = items_continuous,
                         denominator = "model_implied", ci_method = "ml",
                         wald_df = df)
  est <- r$value[r$term == "estimate"]
  se <- r$value[r$term == "se"]
  crit <- stats::qt(0.975, df)
  # The estimate and its standard error do not change; the limits do.
  expect_identical(r$value[r$term %in% c("estimate", "se")],
                   r0$value[r0$term %in% c("estimate", "se")])
  expect_equal(r$value[r$term == "ci_lower"], est - crit * se,
               tolerance = 1e-14)
  expect_equal(r$value[r$term == "ci_upper"], est + crit * se,
               tolerance = 1e-14)
  expect_identical(attr(r, "wald_df"), df)
  expect_true(any(grepl("t reference on 12 degrees of freedom",
                        capture.output(print(r)))))
  expect_false(any(grepl("t reference", capture.output(print(r0)))))

  # The logistic variant applies the t quantile on the logit scale.
  rl <- reliability_omega(data = items_continuous,
                          denominator = "model_implied",
                          ci_method = "mlr_logistic", wald_df = df)
  est <- rl$value[rl$term == "estimate"]
  se <- rl$value[rl$term == "se"]
  logit_se <- se / (est * (1 - est))
  limits <- stats::plogis(stats::qlogis(est) + c(-1, 1) * crit * logit_se)
  expect_equal(rl$value[rl$term %in% c("ci_lower", "ci_upper")],
               limits, tolerance = 1e-12)

  # Under FIML the same: the FIML standard error with the t quantile.
  d <- as.data.frame(items_continuous)
  d$y2[seq(1, N, by = 4)] <- NA
  rf <- reliability_omega(data = d, denominator = "model_implied",
                          ci_method = "ml", missing = "fiml", wald_df = df)
  est <- rf$value[rf$term == "estimate"]
  se <- rf$value[rf$term == "se"]
  expect_equal(rf$value[rf$term == "ci_upper"], est + crit * se,
               tolerance = 1e-14)
})

test_that("wald_df is refused where the interval is not a Wald interval", {
  expect_error(reliability_omega(data = items_continuous, wald_df = 20),
               "robust omega \\(denominator = \"observed\"\\) has a bootstrap")
  expect_error(reliability_omega(data = items_continuous,
                                 denominator = "model_implied",
                                 ci_method = "bonett", wald_df = 20),
               "ci_method = \"bonett\" is not one")
  expect_error(reliability_omega(data = items_continuous,
                                 denominator = "model_implied",
                                 ci_method = "bca", wald_df = 20),
               "ci_method = \"bca\" is not one")
  expect_error(reliability_omega(data = items_continuous,
                                 denominator = "model_implied",
                                 ci_method = "likelihood", wald_df = 20),
               "ci_method = \"likelihood\" is not one")
  expect_error(reliability_omega(data = items_continuous,
                                 denominator = "model_implied",
                                 ci_method = "ml", wald_df = 0),
               "single positive number")
  # Robust omega without an interval records the normal reference.
  r <- suppressMessages(reliability_omega(data = items_continuous))
  expect_identical(attr(r, "wald_df"), Inf)
})


test_that("wald_df = 'missing_information' on complete data is t on N - 1", {
  # With no value missing omega's fraction of missing information is 0
  # exactly, so the calibrated rule max(3, N(1 - FMI) - 1) is N - 1 and the
  # interval is the one wald_df = N - 1 gives, to the last bit: from a
  # covariance matrix ("ml") and from raw data ("ml", and "mlr" with its
  # robust.sem standard error).
  K <- "missing_information"
  a <- reliability_omega(S = S_continuous, N = N,
                         denominator = "model_implied", ci_method = "ml",
                         wald_df = K)
  b <- reliability_omega(S = S_continuous, N = N,
                         denominator = "model_implied", ci_method = "ml",
                         wald_df = N - 1)
  expect_identical(a$value, b$value)
  expect_identical(attr(a, "wald_df"), N - 1)
  tab <- attr(a, "wald_df_table")
  expect_s3_class(tab, "dmar_df_missing")
  expect_identical(tab$term, "omega")
  expect_identical(c(tab$fmi, tab$N, tab$k), c(0, N, 1))
  expect_null(attr(b, "wald_df_table"))
  printed <- capture.output(print(a))
  expect_true(any(grepl(paste0(
    "Wald interval: t reference on 249 degrees of freedom, from von ",
    "Hippel's \\(2016\\) max\\(3, N\\(1 - FMI\\) - 1\\) with N = 250 and ",
    "FMI = 0 \\(wald_df = \"missing_information\"\\), in place of the ",
    "normal"), printed)))
  for (m in c("ml", "mlr")) {
    a <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m,
                           wald_df = K)
    b <- reliability_omega(data = items_continuous,
                           denominator = "model_implied", ci_method = m,
                           wald_df = N - 1)
    expect_identical(a$value, b$value, info = m)
    expect_identical(attr(a, "wald_df"), N - 1, info = m)
  }
  # Listwise deletion analyzes the complete cases, on their N - 1.
  d <- as.data.frame(items_continuous)
  d$y2[seq(1, N, by = 5)] <- NA
  lw <- reliability_omega(data = d, denominator = "model_implied",
                          ci_method = "ml", wald_df = K,
                          missing = "listwise")
  expect_identical(attr(lw, "wald_df"), N - N %/% 5 - 1)
})

test_that("wald_df = 'missing_information' under FIML is max(3, N(1 - FMI) - 1)", {
  skip_on_cran()  # FIML fits with the FMI's two refits
  d <- as.data.frame(items_continuous)
  d$y2[seq(1, N, by = 4)] <- NA
  d$y5[seq(2, N, by = 3)] <- NA
  # The fit behind the coefficient is the single-factor model of cfa_1(),
  # fit by FIML to the items under positional names.
  fit <- cfa_1(d, missing = "fiml", output = "fit")
  fm <- fraction_missing_information(fit)
  fmi <- fm$fmi[fm$term == "omega_f1"]
  df_expected <- max(3, N * (1 - fmi) - 1)
  expect_true(fmi > 0.05 && df_expected < N - 10)
  for (m in c("ml", "mlr")) {
    r <- reliability_omega(data = d, denominator = "model_implied",
                           ci_method = m, missing = "fiml",
                           wald_df = "missing_information")
    expect_equal(attr(r, "wald_df"), df_expected, tolerance = 1e-10,
                 info = m)
    expect_equal(attr(r, "wald_df_table")$fmi, fmi, tolerance = 1e-10)
    expect_output(print(r), "calibration study")
    expect_match(generics::tidy(r)$wald_calibration_note, "affine-collapse screening")
    expect_match(generics::glance(r)$wald_calibration_note, "affine-collapse screening")
    est <- r$value[r$term == "estimate"]
    se <- r$value[r$term == "se"]
    expect_equal(r$value[r$term == "ci_upper"],
                 est + stats::qt(0.975, attr(r, "wald_df")) * se,
                 tolerance = 1e-14)
    expect_true(any(grepl(paste0("Wald interval: t reference on [0-9.]+ ",
                                 "degrees of freedom, from von Hippel's"),
                          capture.output(print(r)))))
  }
  # With an auxiliary variable the FMI is omega's in the fit that carries
  # it, which the auxiliary lowers.
  d$z <- rowMeans(items_continuous) + stats::rnorm(N, sd = 0.3)
  ra <- reliability_omega(data = d, aux = "z", denominator = "model_implied",
                          ci_method = "ml", wald_df = "missing_information")
  expect_true(attr(ra, "wald_df_table")$fmi < fmi)
})

test_that("reliability_omega() refuses wald_df = 'missing_information' outside its calibration", {
  K <- "missing_information"
  expect_error(reliability_omega(data = items_continuous, wald_df = K),
               "robust omega \\(denominator = \"observed\"\\) has a bootstrap")
  for (m in c("ml_logistic", "mlr_logistic")) {
    expect_error(reliability_omega(data = items_continuous,
                                   denominator = "model_implied",
                                   ci_method = m, wald_df = K),
                 "t reference on the logit scale covered up to \\.9645")
  }
  for (m in c("adf", "adf_logistic")) {
    expect_error(reliability_omega(data = items_continuous,
                                   denominator = "model_implied",
                                   ci_method = m, wald_df = K),
                 "not the asymptotically distribution free one")
  }
  for (m in c("bca", "percentile", "bootstrap_se", "likelihood", "feldt",
              "none")) {
    expect_error(reliability_omega(data = items_continuous,
                                   denominator = "model_implied",
                                   ci_method = m, wald_df = K),
                 paste0("ci_method = \"", m, "\" is not one; the keyword ",
                        "applies to ci_method = \"ml\" and \"mlr\""))
  }
  d <- as.data.frame(items_continuous)
  d$y2[seq(1, N, by = 4)] <- NA
  expect_error(reliability_omega(data = d, denominator = "model_implied",
                                 ci_method = "ml", missing = "two_stage",
                                 wald_df = K),
               "DMAR does not compute for a two-stage fit")
})


# ---- Item names lavaan cannot parse ------------------------------------------

test_that("item names that are not lavaan identifiers give the coefficient, not a convergence failure", {
  # Through 1.0.1 the single-factor fit passed the column names to lavaan,
  # which stopped on "item 1" or "1a"; the error was caught and reported
  # as a model that did not converge. The items now enter the fit as y1,
  # ..., yJ, and every name in the result is the user's.
  d0 <- as.data.frame(items_continuous[, 1:5])
  for (nm in list(paste("item", 1:5), paste0(1:5, "a"))) {
    d <- d0
    names(d) <- nm
    a <- suppressMessages(reliability_omega(d))
    b <- suppressMessages(reliability_omega(d0))
    expect_identical(a$value, b$value)
    expect_true(all(attr(a, "residual_cor")$items %in% nm))
    expect_identical(dimnames(attr(a, "residual_cor")$matrix),
                     list(nm, nm))
    a <- reliability_omega(d, denominator = "model_implied")
    b <- reliability_omega(d0, denominator = "model_implied")
    expect_identical(a$value, b$value)
    a <- reliability_omega(S = cov(d), N = N, denominator = "model_implied",
                           ci_method = "ml")
    b <- reliability_omega(S = cov(d0), N = N, denominator = "model_implied",
                           ci_method = "ml")
    expect_identical(a$value, b$value)
  }
})

test_that("item names lavaan cannot parse work under FIML, two-stage estimation, and the bootstrap", {
  skip_on_cran()
  d0 <- as.data.frame(items_continuous[, 1:4])
  d0$y2[1:25] <- NA
  d0$z <- eta + stats::rnorm(N)
  d <- d0
  names(d) <- c(paste("item", 1:4), "aux z")
  for (miss in c("fiml", "two_stage")) {
    a <- reliability_omega(d, missing = miss, aux = "aux z",
                           denominator = "model_implied")
    b <- reliability_omega(d0, missing = miss, aux = "z",
                           denominator = "model_implied")
    expect_identical(a$value, b$value)
    a <- suppressMessages(reliability_omega(d[1:4], missing = miss))
    b <- suppressMessages(reliability_omega(d0[1:4], missing = miss))
    expect_identical(a$value, b$value)
  }
  e <- d[1:4]
  a <- suppressWarnings(reliability_omega(e[stats::complete.cases(e), ],
                                          ci_method = "percentile", B = 30,
                                          seed = 113))
  expect_identical(attr(a, "bootstrap")$failed, 0L)
  # The model implied alpha takes the same route.
  a <- reliability_alpha(e, estimator = "model_implied", missing = "fiml")
  b <- reliability_alpha(d0[1:4], estimator = "model_implied",
                         missing = "fiml")
  expect_identical(a$value, b$value)
  # And coefficient H, which fits the same single-factor model.
  e <- e[stats::complete.cases(e), ]
  a <- reliability_H(data = e)
  b <- reliability_H(data = d0[1:4][stats::complete.cases(d0[1:4]), ])
  expect_identical(a$value, b$value)
})

test_that("duplicated item names are refused by name", {
  d <- as.data.frame(items_continuous[, 1:5])
  names(d) <- c("x", "x", "y", "z", "w")
  expect_error(suppressMessages(reliability_omega(d)),
               "x.{0,3} names more than one column")
  expect_error(reliability_omega(d, denominator = "model_implied"),
               "x.{0,3} names more than one column")
})

test_that("duplicated item names are refused with auxiliary variables too", {
  skip_if_not_installed("lavaan")
  # With aux the items were found by setdiff(), which keeps one column of
  # each name, so the refusal never saw the duplicate: the first of the two
  # x columns was analyzed and the second dropped, and alpha came back as
  # 0.589 with J = 4 for five item columns (0.679 with J = 5 when the two
  # are named apart).
  set.seed(11)
  n <- 200
  f <- stats::rnorm(n)
  X <- sapply(c(.75, .7, .6, .5, .4),
              function(l) l * f + stats::rnorm(n, sd = sqrt(1 - l^2)))
  d <- data.frame(X, a = f + stats::rnorm(n))
  names(d) <- c("x", "x", "y", "z", "w", "a")
  d$y[d$a > 1] <- NA
  msg <- "The item names must be distinct; 'x' names more than one column."
  calls <- list(
    function(d) reliability_alpha(d, aux = "a"),
    function(d) reliability_alpha(d, aux = "a", missing = "two_stage"),
    function(d) reliability_alpha(d, aux = "a", estimator = "model_implied"),
    function(d) reliability_omega(d, aux = "a"),
    function(d) reliability_omega(d, aux = "a", missing = "two_stage"),
    function(d) reliability_omega(d, aux = "a",
                                  denominator = "model_implied"),
    function(d) reliability_H(data = d, aux = "a"))
  for (g in calls) {
    expect_error(g(d), msg, fixed = TRUE)
    # Matrix input keeps the duplicated names through as.data.frame().
    expect_error(g(as.matrix(d)), msg, fixed = TRUE)
  }
  # An auxiliary variable's name given to an item column as well.
  e <- d
  names(e) <- c("a", "x2", "y", "z", "w", "a")
  expect_error(reliability_alpha(e, aux = "a"),
               paste0("The column names of 'data' must be distinct; 'a' ",
                      "names more than one column ('a' named in 'aux')."),
               fixed = TRUE)
  # So does the fit itself.
  expect_error(.omega_fit_cfa(data = e, aux = "a", missing = "ml"),
               "'a' names more than one column", fixed = TRUE)
  # Named apart, the five items are analyzed.
  names(d)[1:2] <- c("x1", "x2")
  r <- suppressWarnings(reliability_alpha(d, aux = "a"))
  expect_identical(r$value[r$term == "n_items"], 5)
})

test_that("lavaan's warnings name the items and auxiliary variables as the data do", {
  skip_if_not_installed("lavaan")
  # The fits enter the items into lavaan as y1, ..., yJ (and the EM of the
  # saturated model as v1, v2, ...), so the warnings lavaan lets through
  # named them that way: "variables involved are: y2 y4" where 1.0.1 said
  # "dep wor" for the listwise and FIML fits, and the FIML saturated
  # estimate said "v2 v4" even then. Each warning is now raised again with
  # the caller's names.
  set.seed(11)
  n <- 120
  f <- stats::rnorm(n)
  d <- as.data.frame(sapply(c(.75, .7, .6, .5), function(l) {
    l * f + stats::rnorm(n, sd = sqrt(1 - l^2))
  }))
  names(d) <- c("anx", "dep", "str", "wor")
  d$wor <- 2 * d$dep + 1
  e <- d
  e$anx[c(3, 9, 20)] <- NA
  e$str[c(5, 40)] <- NA
  e$z <- stats::rnorm(n)
  involved <- function(expr) {
    w <- testthat::capture_warnings(try(suppressMessages(expr),
                                        silent = TRUE))
    hits <- regmatches(w, regexpr("variables involved are:[^\n]*", w))
    unique(trimws(sub("variables involved are:", "", hits)))
  }
  calls <- list(
    function() reliability_omega(d),
    function() reliability_omega(d, denominator = "model_implied"),
    function() reliability_alpha(d, estimator = "model_implied"),
    function() reliability_H(data = d),
    function() reliability_omega(d, denominator = "model_implied",
                                 ci_method = "likelihood"),
    function() .omega_fit_cfa_ts_estimates(e, items = names(e)[1:4]))
  for (g in calls) expect_identical(involved(g()), "dep wor")
  incomplete <- list(
    function() reliability_omega(e[1:4], missing = "fiml"),
    function() reliability_omega(e, aux = "z"),
    function() reliability_omega(e[1:4], missing = "two_stage"),
    function() reliability_omega(e, missing = "two_stage", aux = "z"),
    function() reliability_alpha(e[1:4], missing = "fiml"),
    function() reliability_alpha(e, aux = "z"))
  for (g in incomplete) {
    err <- tryCatch(g(), error = conditionMessage)
    expect_match(err, "dep", fixed = TRUE)
    expect_match(err, "wor", fixed = TRUE)
    expect_no_match(err, "\\b[vy][1-9]\\b")
    expect_match(err, "affine", fixed = TRUE)
  }
  # Each name once, whole words, all at once: items the caller named y2
  # and y1, in that order, are not exchanged twice.
  sw <- d
  names(sw) <- c("y2", "y1", "y3", "y4")
  expect_identical(involved(reliability_omega(sw)), "y1 y4")
  sw$y2[1:5] <- NA
  expect_error(reliability_alpha(sw, missing = "fiml"), "y1.*y4")
  # The categorical fit maps the item names only.
  set.seed(113)
  f <- stats::rnorm(300)
  o <- as.data.frame(sapply(c(.8, .7, .6, .6), function(l) {
    cut(l * f + stats::rnorm(300, sd = sqrt(1 - l^2)),
        c(-Inf, -1, 0, 1, Inf), labels = FALSE)
  }))
  names(o) <- c("anx", "dep", "str", "wor")
  o$wor <- o$dep
  o$str[1:20] <- NA
  w <- testthat::capture_warnings(.omega_c_parts(o, missing = "pairwise"))
  expect_match(w, "correlation between variables wor and dep is",
               all = FALSE, fixed = TRUE)
  expect_false(any(grepl("\\by[0-9]\\b", w)))
})

test_that("a fit that stops with an error is not reported as one that did not converge", {
  # lavaan's own message is quoted; a fit that ran but did not converge
  # keeps the message it had.
  failed <- DMAR:::.omega_fit_none(5L, paste0("y", 1:5),
                                   error = "lavaan->lav_x():\n  something")
  expect_identical(
    DMAR:::.relia_fit_failure(failed, estimator = "MLR"),
    paste0("The single-factor CFA model could not be fit with estimator = ",
           "'MLR'; lavaan stopped with this error: lavaan->lav_x(): ",
           "something"))
  expect_identical(
    DMAR:::.relia_fit_failure(DMAR:::.omega_fit_none(5L, paste0("y", 1:5)),
                              missing = "fiml"),
    paste0("The single-factor CFA model did not converge.",
           DMAR:::.relia_convergence_note("fiml")))
  # An error that is not lavaan's is not attributed to it.
  expect_identical(
    DMAR:::.relia_fit_failure(DMAR:::.omega_fit_none(
      5L, paste0("y", 1:5), error = "a check", from_lavaan = FALSE)),
    "The single-factor CFA model could not be fit: a check")
})

test_that("a fit that ran and did not converge is reported as one, with the advice for its method", {
  skip_if_not_installed("lavaan")
  # fa_confirmatory() stops on a fit that did not converge, and the reliability fit
  # recorded every error as lavaan's: through the first versions of 1.1.0
  # such a fit was reported as "could not be fit; lavaan stopped with this
  # error: The model did not converge." and lost the advice for an
  # incomplete data method. fa_confirmatory() now classes its two stops, and the
  # fit tells them apart.
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  # lavaan's own error is quoted without fa_confirmatory()'s "Model fitting with
  # lavaan failed: ".
  verbal <- c("t5_general_information", "t6_paragraph_comprehension",
              "t7_sentence", "t8_word_classification")
  ch <- holzinger_swineford[, verbal]
  ch$t7_sentence <- as.character(ch$t7_sentence)
  fc <- .omega_fit_cfa(data = ch, se = "none", missing = "two.stage")
  expect_match(fc$fit_error, "^lavaan->lav_data_full\\(\\):\\s+Some model")
  expect_true(fc$fit_error_lavaan)
  # A fit lavaan reports as not converged, on every path that fits
  # through fa_confirmatory(), with lavaan's convergence flag set to FALSE (lavaan's
  # optimizer reaches that only after thousands of iterations, which a
  # real case below shows once).
  d <- holzinger_swineford[, verbal]
  d$t6_paragraph_comprehension[1:30] <- NA
  real <- lavaan::lavInspect
  testthat::local_mocked_bindings(
    lavInspect = function(object, what = "free", ...) {
      if (identical(what, "converged")) FALSE else real(object, what, ...)
    }, .package = "lavaan")
  fit <- quiet(.omega_fit_cfa(data = d, missing = "ml", se = "none"))
  expect_false(fit$converged)
  expect_null(fit$fit_error)
  fiml <- DMAR:::.relia_convergence_note("fiml")
  two <- DMAR:::.relia_convergence_note("two_stage")
  expect_error(quiet(reliability_omega(d, missing = "fiml")),
               paste0("The single-factor CFA model did not converge.", fiml),
               fixed = TRUE)
  expect_error(quiet(reliability_omega(d, missing = "fiml",
                                       denominator = "model_implied")),
               paste0("The single-factor CFA model did not converge with ",
                      "estimator = 'MLR'.", fiml), fixed = TRUE)
  expect_error(quiet(reliability_omega(d, missing = "two_stage")),
               paste0("The single-factor CFA model did not converge.", two),
               fixed = TRUE)
  expect_error(quiet(reliability_alpha(d, estimator = "model_implied",
                                       missing = "two_stage")),
               paste0("The tau-equivalent single-factor model did not ",
                      "converge with estimator = 'ML'.", two), fixed = TRUE)
  expect_error(quiet(reliability_H(data = d, missing = "fiml")),
               paste0("The single-factor model did not converge.", fiml),
               fixed = TRUE)
  # Listwise deletion gives no advice for an incomplete data method.
  expect_error(quiet(reliability_omega(d, missing = "listwise")),
               "^The single-factor CFA model did not converge\\.$")
})

test_that("a fit lavaan's optimizer leaves unconverged is reported as one", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # Four items on ten cases: lavaan's single-factor fit stops at its
  # iteration limit without converging (a few seconds).
  set.seed(10)
  d <- as.data.frame(matrix(stats::rnorm(40), 10, 4))
  expect_error(suppressWarnings(suppressMessages(reliability_omega(d))),
               "^The single-factor CFA model did not converge\\.$")
})

test_that("robust omega's advice about its interval waits for the estimate", {
  # A call that stops on its fit no longer first advises how to request
  # the interval.
  expect_error(reliability_omega(data = triad_neg),
               "no solution for these three items")
  msgs <- testthat::capture_messages(
    try(reliability_omega(data = triad_neg), silent = TRUE))
  expect_length(msgs, 0L)
  expect_message(reliability_omega(data = items_continuous),
                 "Robust omega is reported without a confidence interval")
})

test_that("reliability_omega() validates B when a bootstrap is asked for", {
  expect_error(reliability_omega(data = items_continuous,
                                 ci_method = "percentile", B = 1),
               "'B' must be a single whole number of at least 2",
               fixed = TRUE)
  expect_error(reliability_omega(data = items_continuous,
                                 ci_method = "bca", B = NA),
               "'B' must be a single whole number of at least 2",
               fixed = TRUE)
})

test_that("listwise deletion that leaves fewer than two rows says so before any fit, for every coefficient fitted from data", {
  # No complete row: lavaan printed its table of the variables, and the
  # family reported its error (through 1.0.1, as a model that did not
  # converge); the analytic alpha returned NA with a warning that NaNs
  # were produced. Each now stops before any fit, names the cause, and
  # offers full information maximum likelihood.
  # Listwise deletion is asked for by name: since 1.1.0 the default
  # analyzes such data by FIML.
  e <- as.data.frame(items_continuous[1:60, ])
  for (i in seq_len(nrow(e))) e[i, (i %% ncol(e)) + 1L] <- NA
  lw <- "listwise"
  calls <- list(
    function(d) reliability_omega(d, missing = lw),
    function(d) reliability_omega(d, denominator = "model_implied",
                                  missing = lw),
    function(d) reliability_alpha(d, missing = lw),
    function(d) reliability_alpha(d, estimator = "model_implied",
                                  missing = lw),
    function(d) reliability_H(data = d, missing = lw))
  msg0 <- paste0("No row of 'data' is observed on every item, so listwise ",
                 "deletion leaves nothing to analyze. Use ",
                 "missing = \"fiml\", which uses the incomplete rows as well ",
                 "and assumes the values missing at random")
  for (f in calls) {
    out <- utils::capture.output(err <- tryCatch(
      withCallingHandlers(f(e), warning = function(w) stop(w)),
      error = conditionMessage))
    expect_length(out, 0L)
    expect_match(err, msg0, fixed = TRUE)
  }
  # One complete row is no better.
  e1 <- e
  e1[1, ] <- items_continuous[1, ]
  for (f in calls) {
    expect_error(f(e1), "Only one row of 'data' is observed on every item",
                 fixed = TRUE)
  }
  # The data path FIML uses is the one the message offers.
  expect_true(is.finite(suppressMessages(
    reliability_omega(e, missing = "fiml"))$value[1]))
  # FIML is not offered when a pair of items is never observed together.
  g <- as.data.frame(items_continuous[1:60, ])
  g$y1[1:30] <- NA
  g$y2[31:60] <- NA
  expect_error(reliability_omega(g, missing = lw),
               paste0("Full information maximum likelihood (missing = ",
                      "\"fiml\") cannot help either: no row observes both ",
                      "y1 and y2"), fixed = TRUE)
  # An item no row observes is named, whatever the treatment.
  h <- as.data.frame(items_continuous[1:60, ])
  h$y4 <- NA_real_
  expect_error(reliability_alpha(h),
               "No row observes y4, so neither its variance", fixed = TRUE)
  # Two complete rows are no more than the items, which the single-factor
  # fit needs (lavaan refused the singular matrix, and the refusal quoted
  # it); the next test says more.
  e2 <- e1
  e2[2, ] <- items_continuous[2, ]
  expect_error(reliability_omega(e2, missing = lw),
               "Only 2 rows of 'data' are observed on every item, no more ",
               fixed = TRUE)
})

test_that("no more complete rows than items stops the coefficients fitted by the single-factor model, and alpha and KR-20 take two", {
  skip_if_not_installed("lavaan")
  # N complete rows give a covariance matrix of rank at most N - 1, so with
  # N <= J it is singular and the maximum likelihood fit of the
  # single-factor model is not defined: lavaan refused the matrix (quoted
  # through the first versions of 1.1.0), stopped short of a maximum (a
  # report that the model did not converge), or now and then returned a
  # coefficient. Alpha and KR-20 need only the trace and the sum of the
  # matrix and are computed from two or more complete rows.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence")
  for (nc in c(2L, 5L)) {
    d <- holzinger_swineford[1:40, items]
    for (i in (nc + 1L):40L) d[i, (i %% 5L) + 1L] <- NA
    lead <- paste0("Only ", nc, " rows of 'data' are observed on every ",
                   "item, no more than the 5 items, so the covariance ",
                   "matrix of the complete rows is singular (its rank is ",
                   "at most ", nc - 1L, ")")
    # Listwise deletion is asked for by name: since 1.1.0 the default
    # analyzes such data by FIML.
    lw <- "listwise"
    calls <- list(
      omega = function() reliability_omega(d, missing = lw),
      omega = function() reliability_omega(d, missing = lw,
                                           denominator = "model_implied"),
      `the model implied alpha` = function() {
        reliability_alpha(d, estimator = "model_implied", missing = lw)
      },
      H = function() reliability_H(data = d, missing = lw))
    for (k in seq_along(calls)) {
      out <- utils::capture.output(
        err <- tryCatch(withCallingHandlers(calls[[k]](), warning = stop),
                        error = conditionMessage))
      expect_length(out, 0L)
      expect_match(err, lead, fixed = TRUE)
      expect_match(err, paste0("single-factor model behind ",
                               names(calls)[k], " needs it positive"),
                   fixed = TRUE)
      if (k == 1L) {
        expect_match(err, "cannot help either: The unrestricted normal likelihood", fixed = TRUE)
      } else {
        expect_match(err, 'Use missing = "fiml"', fixed = TRUE)
      }
    }
    # Alpha is computed, with Feldt's interval at two rows.
    a <- suppressWarnings(reliability_alpha(d, ci_method = "feldt",
                                            missing = lw))
    expect_identical(a$value[a$term == "N"], as.numeric(nc))
  }
  # Complete data with no more cases than items.
  expect_error(reliability_omega(holzinger_swineford[1:4, items]),
               paste0("'data' has only 4 rows, no more than the 5 items, so ",
                      "the covariance matrix of the items is singular"),
               fixed = TRUE)
  expect_error(reliability_omega(holzinger_swineford[1:4, items]),
               "The model needs more cases than items.", fixed = TRUE)
  # Six complete rows for five items reach the fit.
  d6 <- holzinger_swineford[1:40, items]
  for (i in 7:40) d6[i, (i %% 5L) + 1L] <- NA
  r6 <- suppressWarnings(reliability_omega(d6, missing = "listwise",
                                           denominator = "model_implied"))
  expect_identical(r6$value[r6$term == "N"], 6)
  # At two rows the intervals of alpha that are not defined there are
  # refused, each for its reason, and every refusal names Feldt's
  # interval, exact under multivariate normality and parallel items (for
  # KR-20, an approximation). Hakstian and Whalen's remains, with an upper
  # limit of one for every data set at a 95% level, and KR-20, alpha of
  # 0/1 items, has the same Feldt interval.
  d2 <- holzinger_swineford[c(1, 50), items]
  why <- c(bonett = "infinite at N = 2",
           ml = "zero for any data", ml_logistic = "zero for any data",
           adf = paste0("from two rows each case's centered cross products ",
                        "are proportional to the covariance matrix, and ",
                        "alpha does not change when that matrix is ",
                        "rescaled, so the distribution free variance of ",
                        "alpha is zero for any data"),
           adf_logistic = "is zero for any data (a computed negative value",
           percentile = "A bootstrap interval requires at least 3 rows",
           bootstrap_se = "A bootstrap interval requires at least 3 rows",
           fisher = "requires 'N' of at least 4")
  feldt <- paste0("With two rows use ci_method = \"feldt\", whose F ",
                  "reference is exact for any 'N' of at least 2 under ",
                  "multivariate normality and parallel items.")
  for (m in names(why)) {
    err <- tryCatch(reliability_alpha(d2, ci_method = m, B = 20),
                    error = conditionMessage)
    expect_match(err, why[[m]], fixed = TRUE)
    if (m != "fisher") {
      expect_true(endsWith(err, feldt))
      expect_false(grepl("hakstian_whalen", err, fixed = TRUE))
    }
  }
  expect_error(reliability_alpha(d2), feldt, fixed = TRUE)
  expect_error(reliability_alpha(S = stats::cov(d2), N = 2,
                                 ci_method = "ml"),
               "zero for any data", fixed = TRUE)
  for (m in c("feldt", "hakstian_whalen")) {
    r <- suppressWarnings(reliability_alpha(d2, ci_method = m))
    expect_true(all(is.finite(r$value[r$term %in% c("ci_lower",
                                                   "ci_upper")])))
  }
  hw <- suppressWarnings(reliability_alpha(d2, ci_method = "hakstian_whalen"))
  expect_identical(hw$value[hw$term == "ci_upper"], 1)
  # The cube root upper limit is negative whatever the data: its factor
  # 1 - z sqrt(2J/(9(J - 1)))/(7/9) is below zero for every J at .95.
  J <- 2:200
  expect_lt(max(1 - stats::qnorm(0.975) *
                  sqrt(2 * J / (9 * (J - 1))) / (7 / 9)), 0)
  # The same through the interval itself, unclamped, for every J from 2 to
  # 50 and estimates below one; at an estimate of one, which two rows can
  # give, the interval is [1, 1].
  for (J in 2:50) {
    for (est in c(-3, -0.5, 0, 0.5, 0.95, 0.999)) {
      expect_gt(.ci_hakstian_whalen(est, 2, J, 0.95)$upper, 1)
    }
    one <- .ci_hakstian_whalen(1, 2, J, 0.95)
    expect_equal(c(one$lower, one$upper), c(1, 1))
  }
  b <- data.frame(a = c(1, 0), b = c(0, 1), c = c(1, 0))
  k <- suppressWarnings(reliability_kr20(b))
  a <- suppressWarnings(reliability_alpha(b, ci_method = "feldt"))
  expect_equal(k$value, a$value, tolerance = 1e-12)
  kr20_feldt <- paste0("With two rows use ci_method = \"feldt\", whose F ",
                       "reference approximates the sampling distribution of ",
                       "KR-20 (Feldt, 1965) for any 'N' of at least 2.")
  for (m in c("bonett", "ml", "adf", "percentile")) {
    err <- tryCatch(reliability_kr20(b, ci_method = m, B = 20),
                    error = conditionMessage)
    expect_true(endsWith(err, kr20_feldt))
  }
  # A bootstrap whose caller has no interval for two rows names none.
  err <- tryCatch(DMAR:::.bootstrap_ci(d2, function(x) 1, B = 10,
                                       conf_level = 0.95,
                                       kind = "percentile", seed = NULL),
                  error = conditionMessage)
  expect_true(endsWith(err, "so the interval would have no width."))
})

test_that("lavaan's error is quoted with the caller's names for the items and the auxiliary variables", {
  # The fits enter the items into lavaan as y1, ..., yJ and the auxiliary
  # variables as a1, ..., aK, so a column name lavaan cannot parse does not
  # break the syntax, and lavaan names the variables that way when it
  # stops. Through the first versions of 1.1.0 its error was quoted as it
  # came (an item of character type was "y2"); each name is now mapped
  # back to the caller's, as a whole word, so the tenth item's y10 is not
  # read as y1 followed by 0.
  nm <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
          "t4_lozenges", "t5_general_information",
          "t6_paragraph_comprehension", "t7_sentence",
          "t8_word_classification", "t9_word_meaning", "t10_addition",
          "t11_code")
  hs <- holzinger_swineford[, nm]
  quoted <- function(err, item) {
    expect_match(err, paste0("lavaan stopped with this error: .*Some model ",
                             "variables are of type character: \"", item,
                             "\""))
    expect_false(grepl("\\b[ya][0-9]+\\b", err))
  }
  ch <- hs[, 1:4]
  ch$t2_cubes <- as.character(ch$t2_cubes)
  # The functions refuse an item that is not numeric before any fit, by
  # name (.relia_check_numeric_items()), so the fits are called directly,
  # as the functions call them.
  fits <- list(
    function(d) .omega_fit_cfa(data = d, se = "none",
                               missing = "two.stage"),
    function(d) .omega_fit_cfa(data = d, se = "none", equal_loadings = TRUE,
                               missing = "two.stage"),
    function(d) .omega_fit_cfa(data = cbind(d, t9_word_meaning =
                                              hs$t9_word_meaning),
                               aux = "t9_word_meaning", se = "none",
                               missing = "ml"))
  for (f in fits) {
    quoted(.relia_fit_failure(suppressWarnings(f(ch))), "t2_cubes")
  }
  # The tenth of eleven items: lavaan's y10 is t10_addition, not
  # t1_visual_perception followed by 0.
  ten <- hs
  ten$t10_addition <- as.character(ten$t10_addition)
  quoted(.relia_fit_failure(.omega_fit_cfa(data = ten, se = "none",
                                           missing = "two.stage")),
         "t10_addition")
  # Names are replaced all at once: items the caller named y2 and y1, in
  # that order, are not exchanged twice.
  sw <- hs[, 1:4]
  names(sw) <- c("y2", "y1", "y3", "y4")
  sw$y1 <- as.character(sw$y1)
  expect_match(.omega_fit_cfa(data = sw, se = "none",
                              missing = "two.stage")$fit_error,
               "of type character: \"y1\"", fixed = TRUE)
  # An auxiliary variable is named too. The functions refuse an auxiliary
  # that is not numeric before any fit, so the fit is called directly.
  ca <- cbind(ch[, -2L], t2_cubes = hs$t2_cubes,
              t9_word_meaning = as.character(hs$t9_word_meaning))
  fa <- suppressWarnings(.omega_fit_cfa(data = ca, aux = "t9_word_meaning",
                                        missing = "ml"))
  expect_false(fa$converged)
  expect_match(fa$fit_error, "of type character: \"t9_word_meaning\"",
               fixed = TRUE)
  # The mapping itself: whole words only, simultaneous, and anything else
  # left alone.
  expect_identical(
    .relia_lavaan_names_back(
      "lav_data_full(): y1 y10 \"y2\". y1_x a1,y10 f1",
      c("y1", "y2", "y10", "a1"), c("first", "second", "tenth", "aux")),
    "lav_data_full(): first tenth \"second\". y1_x aux,tenth f1")
  expect_identical(.relia_lavaan_names_back("y1 and y2", c("y1", "y2"),
                                            c("y2", "y1")),
                   "y2 and y1")
  expect_null(.relia_lavaan_names_back(NULL, "y1", "first"))
})


test_that("a bootstrap replication in which an item takes a single value is set aside with the refusal, and lavaan prints nothing", {
  # Regression test. Through version 1.0.1 a resample that left out the
  # one case on which an item differs from the rest reached lavaan, which
  # printed its table of the variables (y1 to y5, with the item's zero
  # variance marked) to the console for every such resample, about 130
  # lines for these data, and the replication was counted as non-finite
  # without a reason. It is now refused before the fit, by name, as the
  # data would be, and counted as an error that quotes the refusal.
  skip_if_not_installed("lavaan")
  set.seed(1)
  n <- 60
  f <- stats::rnorm(n)
  d <- as.data.frame(sapply(c(.8, .7, .6, .5), function(l) {
    l * f + stats::rnorm(n, sd = sqrt(1 - l^2))
  }))
  names(d) <- paste0("x", 1:4)
  d$x5 <- 0
  d$x5[which.max(f)] <- 1
  set.seed(3)
  bb <- boot::boot(d, function(x, i) sum(x$x5[i]), R = 20)
  lost <- sum(bb$t[, 1] == 0)
  expect_gt(lost, 0L)
  for (den in c("observed", "model_implied")) {
    out <- utils::capture.output(
      r <- suppressWarnings(reliability_omega(d, denominator = den,
                                              ci_method = "percentile",
                                              B = 20, seed = 3)))
    expect_identical(out, character(0))
    bt <- attr(r, "bootstrap")
    expect_identical(bt$reasons[["error"]], lost)
    expect_identical(
      names(bt$error_messages),
      paste0("x5 takes a single value in this resample, so its variance ",
             "is zero, and an item that does not vary is refused."))
  }
})


test_that("the warnings behind wald_df = \"missing_information\" name the caller's variables and a remedy the caller can follow", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # Isolate the public name/remedy translation from the independently
  # tested diagnostic. The former numerical fixture has no finite H1
  # maximum and now correctly stops before reaching this warning.
  testthat::local_mocked_bindings(.wald_df_resolve = function(...) {
    warning(structure(
      list(message = "maximum diagnostic", call = NULL, lhs = "y4",
           op = "~~", rhs = "a1", from = .2, to = .21, gain = .001,
           units = TRUE),
      class = c("dmar_fmi_not_at_maximum", "warning", "condition")))
    list(df = 42, table = NULL)
  })
  w <- testthat::capture_warnings(
    r <- .relia_wald_df("missing_information", "reliability_omega",
      fit = list(fit = TRUE, coef_name = "omega"), N = 100, J = 5,
      items = paste0("q", 1:5, "_resp"), aux = "aux_income"))
  expect_length(w, 1L)
  expect_match(w, "the covariance of q4_resp and aux_income", fixed = TRUE)
  expect_match(w, "Divide the items by one common constant", fixed = TRUE)
  expect_no_match(w, "rescale_data", fixed = TRUE)
  expect_false(any(grepl("\\b[ya][0-9]+\\b", w)))
  expect_identical(r$df, 42)
})



test_that("a flat-likelihood FMI warning is not recast as a units remedy", {
  # The FMI check distinguishes a material change on a flat likelihood
  # from one caused by units. The reliability caller must retain that
  # distinction while translating the fit's internal item names.
  testthat::local_mocked_bindings(.wald_df_resolve = function(...) {
    warning(structure(
      list(message = "maximum diagnostic", call = NULL, lhs = "y1",
           op = "~~", rhs = "a1", from = .2, to = .21, gain = .001,
           units = FALSE,
           change = list(term = "y2 ~~ y2", from = .1, to = .11)),
      class = c("dmar_fmi_not_at_maximum", "warning", "condition")))
    list(df = 42, table = NULL)
  })
  w <- testthat::capture_warnings(
    r <- .relia_wald_df("missing_information", "reliability_omega",
      fit = list(fit = TRUE, coef_name = "omega"), N = 100, J = 3,
      items = c("response1", "response2", "response3"), aux = "income"))
  expect_length(w, 1L)
  expect_match(w, "the covariance of response1 and income", fixed = TRUE)
  expect_match(w, "response2 ~~ response2", fixed = TRUE)
  expect_match(w, "flat in some direction", fixed = TRUE)
  expect_false(any(grepl("Divide the items|rescale_data|y1|a1", w)))
  expect_identical(r$df, 42)
})

test_that("the default 'missing' is FIML for incomplete items and leaves complete data as they stand", {
  # DMAR 1.0.x deleted incomplete rows listwise by default. Since 1.1.0
  # the default fits incomplete items by FIML and analyzes complete data
  # and a covariance matrix as they stand, so that no complete data result
  # changes.
  d <- as.data.frame(items_continuous)
  d_na <- d
  d_na$y2[1:40] <- NA
  mi <- function(x, ci_method = "ml", ...) {
    reliability_omega(x, denominator = "model_implied",
                      ci_method = ci_method, ...)
  }
  a <- mi(d)
  expect_identical(attr(a, "missing"), "listwise")
  expect_identical(a$value, mi(d, missing = "listwise")$value)
  # FIML would give the same estimate; its standard error comes from the
  # observed information rather than the expected.
  f <- mi(d, missing = "fiml")
  expect_equal(f$value[f$term == "estimate"], a$value[a$term == "estimate"],
               tolerance = 1e-6)
  expect_identical(attr(reliability_omega(S = S_continuous, N = N,
                                          denominator = "model_implied",
                                          ci_method = "ml"), "missing"),
                   "listwise")
  # Incomplete items: FIML, every case, without a warning.
  expect_no_warning(x <- mi(d_na))
  expect_identical(attr(x, "missing"), "fiml")
  expect_equal(x$value[x$term == "N"], nrow(d_na))
  expect_equal(x$value, mi(d_na, missing = "fiml")$value)
  # Listwise deletion asked for by name is the DMAR 1.0.x analysis.
  expect_no_warning(l <- mi(d_na, missing = "listwise"))
  expect_equal(l$value[l$term == "N"], sum(stats::complete.cases(d_na)))
  # Robust omega takes no closed-form interval.
  expect_error(reliability_omega(d_na, ci_method = "ml"),
               "is available only for denominator = \"model_implied\"",
               fixed = TRUE)
  # Robust omega follows the same rule.
  ro <- suppressMessages(reliability_omega(d_na, ci_method = "none"))
  expect_identical(attr(ro, "missing"), "fiml")
  # An interval FIML cannot give is computed as before with complete data
  # and, with values missing, stops and names the choices rather than
  # deleting cases unasked; an explicit FIML request keeps its own
  # refusal.
  for (m in c("feldt", "adf", "likelihood")) {
    expect_identical(attr(mi(d, ci_method = m), "missing"), "listwise")
    expect_error(mi(d_na, ci_method = m),
                 paste0("not available with ci_method = \"", m, "\", and ",
                        "the items have missing values"), fixed = TRUE)
  }
  expect_error(mi(d_na, ci_method = "feldt"),
               paste0("Use missing = \"listwise\" to analyze the complete ",
                      "cases, or keep FIML with ci_method = \"ml\", ",
                      "\"mlr\", their _logistic variants, or a bootstrap ",
                      "method."), fixed = TRUE)
  expect_error(mi(d_na, ci_method = "feldt", missing = "fiml"),
               "is a complete-data closed form", fixed = TRUE)
  # Robust omega takes no closed-form interval under any treatment, and
  # that is said first, so that the remedy named is never refused next:
  # the FIML refusal would have offered missing = "listwise", which then
  # stops on the denominator.
  for (miss in list(NULL, "fiml", "listwise")) {
    expect_error(do.call(reliability_omega,
                         c(list(d_na, ci_method = "feldt"),
                           if (!is.null(miss)) list(missing = miss))),
                 "is available only for denominator = \"model_implied\"",
                 fixed = TRUE)
  }
  # A row with every item missing does not make the data incomplete:
  # listwise deletion and FIML both drop it, so complete data with a blank
  # row keep the complete data analysis, with no warning.
  blank <- rbind(d, NA)
  expect_no_warning(b <- mi(blank))
  expect_identical(attr(b, "missing"), "listwise")
  expect_identical(b$value, a$value)
  expect_identical(attr(mi(blank, ci_method = "feldt"), "missing"),
                   "listwise")
  expect_identical(attr(mi(d_na, ci_method = "feldt", missing = "listwise"),
                        "missing"), "listwise")
  # Auxiliary variables imply FIML whether or not a value is missing.
  z <- data.frame(d, z = d$y1 + stats::rnorm(nrow(d)))
  expect_no_warning(za <- mi(z, aux = "z"))
  expect_identical(attr(za, "missing"), "fiml")
  # reliability() forwards an unsupplied 'missing' unsupplied.
  expect_no_warning(r <- reliability(d_na, type = "omega",
                                     denominator = "model_implied",
                                     ci_method = "ml"))
  expect_equal(r$value, x$value)
  expect_identical(attr(reliability(d, type = "omega",
                                    denominator = "model_implied",
                                    ci_method = "ml"), "missing"),
                   "listwise")
})

test_that("by default the model implied omega and fa_confirmatory()'s agree on incomplete data", {
  spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
               "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
  d <- holzinger_swineford[spatial]
  quiet <- function(x) suppressWarnings(suppressMessages(x))
  r <- quiet(reliability_omega(d, denominator = "model_implied",
                               ci_method = "ml_logistic"))
  f <- quiet(fa_confirmatory(d, list(f1 = spatial), output = "measurement"))
  w <- f[f$term == "omega_f1", ]
  expect_equal(r$value[match(c("estimate", "se", "ci_lower",
                               "ci_upper", "N"), r$term)],
               c(w$estimate, w$se, w$ci_lower, w$ci_upper, nrow(d)),
               tolerance = 1e-5)
  # And both are lavaan's own FIML fit of the single-factor model, with
  # omega defined in the model and its limits formed on the logit scale
  # from lavaan's delta method standard error.
  skip_if_not_installed("lavaan")
  l <- paste0("l", seq_along(spatial))
  e <- paste0("e", seq_along(spatial))
  syntax <- paste0(
    "f =~ ", paste0(l, "*", spatial, collapse = " + "), "\n",
    paste0(spatial, " ~~ ", e, "*", spatial, collapse = "\n"), "\n",
    "omega := (", paste(l, collapse = " + "), ")^2 / ((",
    paste(l, collapse = " + "), ")^2 + ", paste(e, collapse = " + "), ")")
  fit <- lavaan::cfa(syntax, data = d, missing = "ml", std.lv = TRUE)
  pe <- lavaan::parameterEstimates(fit)
  om <- pe[pe$label == "omega", ]
  limits <- stats::plogis(stats::qlogis(om$est) + c(-1, 1) *
                            stats::qnorm(0.975) * om$se /
                            (om$est * (1 - om$est)))
  expect_equal(r$value[match(c("estimate", "se", "ci_lower",
                               "ci_upper", "N"), r$term)],
               c(om$est, om$se, limits, lavaan::lavInspect(fit, "ntotal")),
               tolerance = 1e-7)
  # With complete items the default here analyzes them as they stand, with
  # the expected information, and fa_confirmatory() fits them by FIML, with
  # the observed information: the estimates agree, and the standard errors
  # differ by 0.4 percent on these four tests (as the help page says).
  four <- spatial[1:4]
  r4 <- reliability_omega(d[four], denominator = "model_implied",
                          ci_method = "ml_logistic")
  f4 <- quiet(fa_confirmatory(d[four], list(f1 = four),
                              output = "measurement"))
  w4 <- f4[f4$term == "omega_f1", ]
  expect_equal(r4$value[r4$term == "estimate"], w4$estimate,
               tolerance = 1e-6)
  expect_equal(r4$value[r4$term == "se"] / w4$se, 0.9956, tolerance = 1e-3)
})
