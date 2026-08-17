# Tests for reliability_alpha(). The point estimate is benchmarked against
# the closed-form Guttman (1945) / Cronbach (1951) formula computed
# directly from a covariance matrix. Confidence interval methods are
# checked for invariants (estimate inside the CI, CI bounded in [0, 1],
# lower < upper, sensitivity to conf_level), not for exact reproduction
# of any single reference number.

set.seed(113)
J <- 6
N <- 200
loadings <- seq(0.4, 0.8, length.out = J)
eta <- rnorm(N)
errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) + errors
colnames(items_continuous) <- paste0("y", seq_len(J))
S_continuous <- cov(items_continuous)

manual_alpha <- function(S) {
  J <- nrow(S)
  (J / (J - 1)) * (1 - sum(diag(S)) / sum(S))
}
alpha_truth <- manual_alpha(S_continuous)


test_that("reliability_alpha() matches the closed-form Guttman/Cronbach formula", {
  r <- reliability_alpha(data = items_continuous, ci_method = "none")
  expect_equal(r$value[r$term == "estimate"], alpha_truth, tolerance = 1e-10)
})

test_that("reliability_alpha() returns the same point estimate from data and S", {
  from_data <- reliability_alpha(data = items_continuous, ci_method = "none")
  from_S    <- reliability_alpha(S = S_continuous, N = N, ci_method = "none")
  expect_equal(
    from_data$value[from_data$term == "estimate"],
    from_S$value[from_S$term == "estimate"],
    tolerance = 1e-10
  )
})

test_that("reliability_alpha() returns a tidy data.frame with the expected rows", {
  r <- reliability_alpha(data = items_continuous)
  expect_s3_class(r, "data.frame")
  expect_named(r, c("term", "value"))
  expect_true(all(c("estimate", "se", "lower_limit", "upper_limit",
                    "conf_level", "N", "J") %in% r$term))
  expect_equal(attr(r, "coefficient"), "alpha")
  expect_equal(attr(r, "ci_method"), "bonett")
})

test_that("reliability_alpha() CI brackets the point estimate (closed-form methods)", {
  for (m in c("feldt", "fisher", "bonett", "hakstian_whalen", "ml", "ml_logistic")) {
    r <- reliability_alpha(data = items_continuous, ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
  }
})

test_that("reliability_alpha() ADF CI brackets the estimate and is in [0, 1]", {
  r <- reliability_alpha(data = items_continuous, ci_method = "adf")
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "lower_limit"]
  ul  <- r$value[r$term == "upper_limit"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
})

test_that("reliability_alpha() bootstrap CIs bracket the estimate", {
  skip_if_not_installed("boot")
  for (m in c("bootstrap_se", "bootstrap_se_logistic", "percentile", "bca")) {
    r <- reliability_alpha(data = items_continuous, ci_method = m, B = 200)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
    expect_equal(attr(r, "B"), 200)
  }
})

test_that("reliability_alpha() higher conf_level gives a wider CI (Bonett)", {
  r80 <- reliability_alpha(data = items_continuous, ci_method = "bonett",
                        conf_level = 0.80)
  r95 <- reliability_alpha(data = items_continuous, ci_method = "bonett",
                        conf_level = 0.95)
  width80 <- r80$value[r80$term == "upper_limit"] -
             r80$value[r80$term == "lower_limit"]
  width95 <- r95$value[r95$term == "upper_limit"] -
             r95$value[r95$term == "lower_limit"]
  expect_gt(width95, width80)
})

test_that("reliability_alpha() bootstrap is reproducible given the same seed", {
  skip_if_not_installed("boot")
  r1 <- reliability_alpha(data = items_continuous, ci_method = "percentile",
                       B = 100, seed = 113)
  r2 <- reliability_alpha(data = items_continuous, ci_method = "percentile",
                       B = 100, seed = 113)
  expect_equal(r1$value, r2$value)
})

test_that("reliability_alpha() ci_method = 'none' returns NA bounds", {
  r <- reliability_alpha(data = items_continuous, ci_method = "none")
  expect_true(is.na(r$value[r$term == "se"]))
  expect_true(is.na(r$value[r$term == "lower_limit"]))
  expect_true(is.na(r$value[r$term == "upper_limit"]))
})

test_that("transformation-based intervals report the se on the coefficient scale", {
  # Decided 2026-08-12: the se row is the delta method back-transform of
  # the transformation-scale standard error, evaluated at the estimate,
  # so it reads on the same scale as the estimate; the transformation
  # scale value travels in its own se_transformed row, with the scale
  # named in the se_transform_scale attribute.
  scales <- c(fisher          = "fisher_z",
              bonett          = "log(1-alpha)",
              hakstian_whalen = "cube_root")
  df_n <- N - 1
  df_d <- (N - 1) * (J - 1)
  for (m in names(scales)) {
    r <- reliability_alpha(data = items_continuous, ci_method = m)
    est  <- r$value[r$term == "estimate"]
    se   <- r$value[r$term == "se"]
    se_t <- r$value[r$term == "se_transformed"]
    expect_identical(attr(r, "se_transform_scale"), scales[[m]])
    ref <- switch(
      m,
      fisher = list(se_t  = sqrt(1 / (N - 3)),
                    deriv = 1 - est^2),
      bonett = list(se_t  = sqrt(2 * J / ((J - 1) * (N - 2))),
                    deriv = 1 - est),
      hakstian_whalen = {
        a <- (1 - 2 / (9 * df_n)) / (1 - 2 / (9 * df_d))
        list(se_t  = sqrt((2 * J / (9 * df_d)) * (1 - est)^(2 / 3) /
                            (1 - 2 / (9 * df_n))^2),
             deriv = 3 * a^3 * (1 - est)^(2 / 3))
      })
    expect_equal(se_t, ref$se_t, tolerance = 1e-12,
                 label = paste("transformation-scale se for", m))
    expect_equal(se, se_t * ref$deriv, tolerance = 1e-12,
                 label = paste("coefficient-scale se for", m))
  }
})

test_that("the coefficient-scale se is pinned per transformation method", {
  # Anchors computed from this file's seed-113 dataset (N = 200, J = 6)
  # at the closed forms above; a drift here means the back-transform
  # (or the transformation-scale formula behind it) changed.
  anchor_se <- c(fisher          = 0.0293876760674221,
                 bonett          = 0.0257073559344799,
                 hakstian_whalen = 0.0256025954609992)
  anchor_se_t <- c(fisher          = 0.0712470499879096,
                   bonett          = 0.1100963765126361,
                   hakstian_whalen = 0.0225668634756272)
  for (m in names(anchor_se)) {
    r <- reliability_alpha(data = items_continuous, ci_method = m)
    expect_equal(r$value[r$term == "se"], anchor_se[[m]],
                 tolerance = 1e-8, label = paste("se anchor for", m))
    expect_equal(r$value[r$term == "se_transformed"], anchor_se_t[[m]],
                 tolerance = 1e-8,
                 label = paste("se_transformed anchor for", m))
  }
})

test_that("non-transformation methods carry no se_transformed row", {
  for (m in c("none", "feldt", "ml", "ml_logistic", "adf")) {
    r <- reliability_alpha(data = items_continuous, ci_method = m)
    expect_false("se_transformed" %in% r$term,
                 label = paste("no se_transformed row for", m))
    expect_null(attr(r, "se_transform_scale"))
  }
  # Feldt is an F pivot with no standard error on either scale.
  r <- reliability_alpha(data = items_continuous, ci_method = "feldt")
  expect_true(is.na(r$value[r$term == "se"]))
})

test_that("the back-transformed se matches the Monte Carlo SD of the coefficient", {
  skip_on_cran()  # a 2 x 3000-replication simulation; the algebraic scale checks above run on CRAN
  # The se row claims the coefficient scale, so across repeated samples
  # the standard deviation of the sample coefficient must match it
  # within Monte Carlo error. Under a parallel-items population the
  # asymptotic standard deviation of the sample coefficient is
  # (1 - alpha) * sqrt(2J / ((J - 1) N)); the Bonett and Hakstian-Whalen
  # transformation-scale standard errors are correct there, so their
  # back-transforms must reproduce the Monte Carlo SD. Fisher's Z scale
  # standard error 1 / sqrt(N - 3) is the correlation formula and equals
  # the asymptotic SD of the transformed coefficient only where
  # (1 + alpha)^2 = 2J / (J - 1), so the Fisher check runs at that
  # population value: there a mismatch would implicate the
  # back-transform rather than the method's documented overcoverage.
  G <- 3000
  N_mc <- 400
  J_mc <- 6
  alpha_hats <- function(l, G) {
    vapply(seq_len(G), function(g) {
      eta <- rnorm(N_mc)
      X <- outer(eta, rep(l, J_mc)) +
        matrix(rnorm(N_mc * J_mc, sd = sqrt(1 - l^2)), N_mc, J_mc)
      DMAR:::.alpha_from_S(cov(X))
    }, numeric(1))
  }
  pop_alpha <- function(l) J_mc * l^2 / (1 + (J_mc - 1) * l^2)

  set.seed(113)
  l1 <- 0.6                     # population alpha 0.771
  a1 <- pop_alpha(l1)
  sd1 <- sd(alpha_hats(l1, G))
  expect_lt(abs(DMAR:::.ci_bonett(a1, N_mc, J_mc, 0.95)$se / sd1 - 1),
            0.05)
  expect_lt(abs(DMAR:::.ci_hakstian_whalen(a1, N_mc, J_mc, 0.95)$se /
                  sd1 - 1),
            0.05)

  a2 <- sqrt(2 * J_mc / (J_mc - 1)) - 1   # 0.549: Fisher's se is correct here
  l2 <- sqrt(a2 / (J_mc - (J_mc - 1) * a2))
  sd2 <- sd(alpha_hats(l2, G))
  expect_lt(abs(DMAR:::.ci_fisher(a2, N_mc, 0.95)$se / sd2 - 1), 0.05)
})

test_that("reliability_alpha() errors when S is supplied without N", {
  expect_error(reliability_alpha(S = S_continuous),
               regexp = "'N' must be supplied")
})

test_that("reliability_alpha() errors when neither data nor S is supplied", {
  expect_error(reliability_alpha(),
               regexp = "Either 'data'.*or.*'S'")
})

test_that("reliability_alpha() errors when data and S are both supplied", {
  expect_error(reliability_alpha(data = items_continuous, S = S_continuous, N = N),
               regexp = "not both")
})

test_that("reliability_alpha() errors when fewer than two items are supplied", {
  expect_error(reliability_alpha(data = items_continuous[, 1, drop = FALSE]),
               regexp = "At least two items")
})

test_that("reliability_alpha() errors when an ADF CI is requested with S only", {
  expect_error(
    reliability_alpha(S = S_continuous, N = N, ci_method = "adf"),
    regexp = "requires raw 'data'"
  )
})

test_that("reliability_alpha() Fisher CI stops for N < 4", {
  # The Fisher's Z variance is 1/(N - 3): infinite at N = 3 (the interval
  # would be vacuous) and undefined below it.
  expect_error(
    reliability_alpha(S = S_continuous, N = 3, ci_method = "fisher"),
    regexp = "at least 4"
  )
  expect_error(
    reliability_alpha(S = S_continuous, N = 3, ci_method = "fisher"),
    regexp = "vacuous"
  )
  expect_error(
    reliability_alpha(S = S_continuous, N = 2, ci_method = "fisher"),
    regexp = "1/\\(N - 3\\)"
  )
})

test_that("reliability_alpha() rejects out-of-range conf_level", {
  expect_error(reliability_alpha(data = items_continuous, conf_level = 1.2),
               regexp = "conf_level")
  expect_error(reliability_alpha(data = items_continuous, conf_level = 0),
               regexp = "conf_level")
})


test_that("the profile likelihood interval brackets the model implied estimate", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two likelihood profiles, each refitting the model many times; the closed-form anchors above run on CRAN
  # The profile likelihood belongs to the model implied estimator, whose
  # point estimate it profiles. Pairing it with the analytic estimate was
  # the defect that merging the two functions removed, so the interval is
  # now requested where the estimate it refers to actually lives, and it
  # must contain that estimate.
  r <- reliability_alpha(data = items_continuous, estimator = "model_implied",
                         ci_method = "likelihood")
  est <- r$value[r$term == "estimate"]
  lo  <- r$value[r$term == "lower_limit"]
  hi  <- r$value[r$term == "upper_limit"]
  expect_true(lo < est && est < hi)
  expect_gte(lo, 0)
  expect_lte(hi, 1)

  # Raw data and covariance input profile the same likelihood.
  r_S <- reliability_alpha(S = S_continuous, N = N,
                           estimator = "model_implied",
                           ci_method = "likelihood")
  expect_equal(r_S$value[r_S$term == "lower_limit"], lo, tolerance = 1e-5)
  expect_equal(r_S$value[r_S$term == "upper_limit"], hi, tolerance = 1e-5)

  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(lo, 0.7409301493224807, tolerance = 5e-4)
  expect_equal(hi, 0.831808628222755, tolerance = 5e-4)
})

test_that("the estimator argument selects between the two routes to alpha", {
  skip_if_not_installed("lavaan")
  # Deliberately unequal loadings, so tau-equivalence is false and the two
  # estimators must disagree; under tau-equivalence they would converge.
  set.seed(113)
  J <- 5; n <- 300
  lam <- c(.8, .7, .6, .5, .4)
  eta <- rnorm(n)
  d <- as.data.frame(outer(eta, lam) +
         matrix(rnorm(n * J), n, J) %*% diag(sqrt(1 - lam^2)))
  names(d) <- paste0("y", seq_len(J))

  a <- reliability_alpha(data = d, estimator = "analytic", ci_method = "none")
  m <- reliability_alpha(data = d, estimator = "model_implied",
                         ci_method = "none")
  est <- function(x) x$value[x$term == "estimate"]

  expect_identical(attr(a, "estimator"), "analytic")
  expect_identical(attr(m, "estimator"), "model_implied")
  expect_identical(attr(a, "coefficient"), "alpha")
  expect_false(isTRUE(all.equal(est(a), est(m))))

  # Each route matches the reference implementation for that estimator.
  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(est(a), 0.7448879415292547, tolerance = 1e-6)
  expect_equal(est(m), 0.7624058296878771, tolerance = 1e-6)
})

test_that("an interval method is refused when its estimator cannot supply it", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  d <- as.data.frame(matrix(rnorm(200 * 4), 200, 4))
  names(d) <- paste0("y", 1:4)

  # The profile likelihood profiles the model implied coefficient, so pairing
  # it with the analytic point estimate would report an interval and an
  # estimate for different quantities: under a misspecified model the
  # interval can exclude the estimate it accompanies. Refused, with the
  # remedy named.
  expect_error(
    reliability_alpha(data = d, estimator = "analytic",
                      ci_method = "likelihood"),
    'estimator = "model_implied"', fixed = TRUE
  )
  # The closed-form intervals for the sample coefficient are equally
  # unavailable to the model implied estimator.
  for (m in c("bonett", "feldt", "fisher", "hakstian_whalen")) {
    expect_error(
      reliability_alpha(data = d, estimator = "model_implied", ci_method = m),
      'estimator = "analytic"', fixed = TRUE
    )
  }
  # As is the robust standard error to the analytic estimator.
  expect_error(
    reliability_alpha(data = d, estimator = "analytic", ci_method = "mlr"),
    'estimator = "model_implied"', fixed = TRUE
  )
  # A value that belongs to neither is reported as such.
  expect_error(reliability_alpha(data = d, ci_method = "not_a_method"),
               "must be one of", fixed = TRUE)

  # Defaults are per estimator and unchanged from the two former functions.
  expect_identical(attr(reliability_alpha(data = d), "ci_method"), "bonett")
  expect_identical(
    attr(reliability_alpha(data = d, estimator = "model_implied"),
         "ci_method"), "mlr")
  # "mlr" needs raw data, so covariance input falls back to "ml".
  expect_identical(
    attr(reliability_alpha(S = cov(d), N = 200, estimator = "model_implied"),
         "ci_method"), "ml")
})
