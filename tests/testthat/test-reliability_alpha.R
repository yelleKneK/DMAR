# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

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
  skip_on_cran()  # The extended validation grid runs in the full local suite.
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
  expect_true(all(c("estimate", "se", "ci_lower", "ci_upper",
                    "conf_level", "N", "n_items") %in% r$term))
  expect_equal(attr(r, "coefficient"), "alpha")
  expect_equal(attr(r, "ci_method"), "bonett")
})

test_that("reliability_alpha() CI brackets the point estimate (closed-form methods)", {
  for (m in c("feldt", "fisher", "bonett", "hakstian_whalen", "ml", "ml_logistic")) {
    r <- reliability_alpha(data = items_continuous, ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
  }
})

test_that("reliability_alpha() ADF CI brackets the estimate and is in [0, 1]", {
  r <- reliability_alpha(data = items_continuous, ci_method = "adf")
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "ci_lower"]
  ul  <- r$value[r$term == "ci_upper"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
})

test_that("reliability_alpha() bootstrap CIs bracket the estimate", {
  skip_on_cran()
  skip_if_not_installed("boot")
  for (m in c("bootstrap_se", "bootstrap_se_logistic", "percentile", "bca")) {
    r <- reliability_alpha(data = items_continuous, ci_method = m, B = 200)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
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
  width80 <- r80$value[r80$term == "ci_upper"] -
             r80$value[r80$term == "ci_lower"]
  width95 <- r95$value[r95$term == "ci_upper"] -
             r95$value[r95$term == "ci_lower"]
  expect_gt(width95, width80)
})

test_that("reliability_alpha() bootstrap is reproducible given the same seed", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
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
  expect_true(is.na(r$value[r$term == "ci_lower"]))
  expect_true(is.na(r$value[r$term == "ci_upper"]))
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
  lo  <- r$value[r$term == "ci_lower"]
  hi  <- r$value[r$term == "ci_upper"]
  expect_true(lo < est && est < hi)
  expect_gte(lo, 0)
  expect_lte(hi, 1)

  # Raw data and covariance input profile the same likelihood.
  r_S <- reliability_alpha(S = S_continuous, N = N,
                           estimator = "model_implied",
                           ci_method = "likelihood")
  expect_equal(r_S$value[r_S$term == "ci_lower"], lo, tolerance = 1e-5)
  expect_equal(r_S$value[r_S$term == "ci_upper"], hi, tolerance = 1e-5)

  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(lo, 0.7409301493224807, tolerance = 5e-4)
  expect_equal(hi, 0.831808628222755, tolerance = 5e-4)
})

test_that("the estimator argument selects between the two routes to alpha", {
  skip_on_cran()
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
  skip_on_cran()  # The extended validation grid runs in the full local suite.
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
  # (The items are independent noise, so the sample alpha here is
  # slightly negative, which is warned about.)
  expect_warning(r <- reliability_alpha(data = d),
                 class = "dmar_negative_reliability")
  expect_identical(attr(r, "ci_method"), "bonett")
  expect_identical(
    attr(reliability_alpha(data = d, estimator = "model_implied"),
         "ci_method"), "mlr")
  # "mlr" needs raw data, so covariance input falls back to "ml".
  expect_identical(
    attr(reliability_alpha(S = cov(d), N = 200, estimator = "model_implied"),
         "ci_method"), "ml")
})


# ---- wald_df: a t reference distribution for the Wald intervals ------------

test_that("wald_df = Inf reproduces the closed-form Wald intervals of alpha exactly", {
  for (m in c("ml", "ml_logistic", "adf", "adf_logistic", "bonett",
              "feldt")) {
    a <- reliability_alpha(data = items_continuous, ci_method = m)
    b <- reliability_alpha(data = items_continuous, ci_method = m,
                           wald_df = Inf)
    expect_identical(a, b)
    expect_identical(attr(a, "wald_df"), Inf)
  }
  r <- reliability_alpha(S = S_continuous, N = N, ci_method = "ml")
  est <- r$value[r$term == "estimate"]
  se <- r$value[r$term == "se"]
  expect_identical(r$value[r$term == "ci_upper"],
                   est + stats::qnorm(0.975) * se)
})

test_that("a finite wald_df puts the t quantile into alpha's Wald intervals", {
  crit <- stats::qt(0.975, 15)
  for (m in c("ml", "adf")) {
    r <- reliability_alpha(data = items_continuous, ci_method = m,
                           wald_df = 15)
    est <- r$value[r$term == "estimate"]
    se <- r$value[r$term == "se"]
    expect_equal(r$value[r$term %in% c("ci_lower", "ci_upper")],
                 est + c(-1, 1) * crit * se, tolerance = 1e-14)
    expect_identical(attr(r, "wald_df"), 15)
  }
  rl <- reliability_alpha(S = S_continuous, N = N,
                          ci_method = "ml_logistic", wald_df = 15)
  est <- rl$value[rl$term == "estimate"]
  se <- rl$value[rl$term == "se"]
  limits <- stats::plogis(stats::qlogis(est) +
                            c(-1, 1) * crit * se / (est * (1 - est)))
  expect_equal(rl$value[rl$term %in% c("ci_lower", "ci_upper")],
               limits, tolerance = 1e-12)
  expect_true(any(grepl("Wald interval: t reference on 15 degrees",
                        capture.output(print(rl)))))
})

test_that("wald_df is refused for alpha's intervals that are not Wald intervals", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  for (m in c("bonett", "feldt", "fisher", "hakstian_whalen",
              "percentile", "none")) {
    expect_error(reliability_alpha(data = items_continuous, ci_method = m,
                                   wald_df = 30),
                 paste0("ci_method = \"", m, "\" is not one"))
  }
  # The default interval is Bonett's, so a finite wald_df needs a Wald
  # interval named.
  expect_error(reliability_alpha(data = items_continuous, wald_df = 30),
               "ci_method = \"bonett\" is not one")
  # The keyword is refused for every interval, with the calibration's
  # verdict and its basis: no rule passed, the t reference came barely
  # closer to .95, and its coverage reached .966 in the one cell where the
  # study's saturated fit converged in 53 percent of the replications.
  for (m in c("ml", "ml_logistic", "adf", "bonett")) {
    expect_error(reliability_alpha(data = items_continuous, ci_method = m,
                                   wald_df = "missing_information"),
                 paste0("not available for reliability_alpha\\(\\): in the ",
                        "calibration study behind the keyword no rule ",
                        "passed for alpha.*\\.0042 against \\.0043.*",
                        "converged in 53 percent.*reached \\.966"))
  }
  expect_error(reliability_alpha(data = items_continuous,
                                 estimator = "model_implied",
                                 wald_df = "missing_information"),
               "not available for reliability_alpha")
})

test_that("wald_df reaches alpha's FIML and model implied Wald intervals", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan fits per case; the closed forms run on CRAN
  d <- as.data.frame(items_continuous)
  d$y3[seq(2, N, by = 5)] <- NA
  crit <- stats::qt(0.975, 15)
  cases <- list(list(data = d, ci_method = "ml", missing = "fiml"),
                list(data = items_continuous, ci_method = "mlr",
                     estimator = "model_implied"))
  for (args in cases) {
    a <- do.call(reliability_alpha, args)
    expect_identical(a, do.call(reliability_alpha,
                                c(args, list(wald_df = Inf))))
    r <- do.call(reliability_alpha, c(args, list(wald_df = 15)))
    est <- r$value[r$term == "estimate"]
    se <- r$value[r$term == "se"]
    expect_identical(c(est, se), a$value[a$term %in% c("estimate", "se")])
    expect_equal(r$value[r$term %in% c("ci_lower", "ci_upper")],
                 est + c(-1, 1) * crit * se, tolerance = 1e-14)
  }
})


# ---- Edge cases of the estimate and of the bootstrap ------------------------

test_that("a negative alpha is warned about and its limits are clamped in order", {
  skip_on_cran()
  # One item keyed in the opposite direction to the other two. The sample
  # alpha is -1.25, and Bonett's unclamped interval, about [-1.86, -0.77],
  # lies wholly below zero. Through 1.0.1 only the lower limit was
  # clamped, so the interval came back inverted, [0, -0.77], with the
  # estimate outside it; each limit is now clamped into [0, 1].
  X <- items_continuous[, c(1, 5, 6)]
  X[, 3] <- -X[, 3]
  est <- manual_alpha(cov(X))
  expect_lt(est, 0)
  expect_lt(DMAR:::.ci_bonett(est, N, 3, 0.95)$upper, 0)
  for (m in c("bonett", "feldt", "ml", "adf")) {
    expect_warning(r <- reliability_alpha(data = X, ci_method = m),
                   class = "dmar_negative_reliability")
    v <- stats::setNames(r$value, r$term)
    expect_equal(v[["estimate"]], est, tolerance = 1e-12)
    expect_identical(unname(v[c("ci_lower", "ci_upper")]), c(0, 0))
  }
  # The warning says what a negative estimate usually means and what
  # became of the limits.
  w <- testthat::capture_warnings(reliability_alpha(data = X))
  expect_length(w, 1L)
  expect_match(w, "is negative (-1.25)", fixed = TRUE)
  expect_match(w, "needs reverse scoring", fixed = TRUE)
  expect_match(w, "reported as [0, 0]", fixed = TRUE)
  # The logit scale has no negative values, so a logistic interval has NA
  # limits, and the warning says why.
  expect_warning(r <- reliability_alpha(data = X, ci_method = "ml_logistic"),
                 "not defined for a negative estimate")
  expect_true(all(is.na(r$value[r$term %in% c("ci_lower",
                                               "ci_upper")])))
  # A bootstrap interval is clamped the same way.
  w <- testthat::capture_warnings(
    r <- reliability_alpha(data = X, ci_method = "percentile", B = 200,
                           seed = 113))
  expect_identical(r$value[r$term %in% c("ci_lower", "ci_upper")],
                   c(0, 0))
  expect_true(any(grepl("is negative", w, fixed = TRUE)))
})

test_that("the result formatter clamps each limit into [0, 1] and keeps the limits in order", {
  lim <- function(est, lo, hi) {
    r <- suppressWarnings(DMAR:::.relia_result(
      estimate = est, lower = lo, upper = hi, conf_level = 0.95, N = 50,
      J = 4, coefficient = "alpha", ci_method = "bonett"))
    r$value[r$term %in% c("ci_lower", "ci_upper")]
  }
  expect_identical(lim(-0.2, -0.5, -0.1), c(0, 0))
  expect_identical(lim(-0.1, -0.4, 0.3), c(0, 0.3))
  expect_identical(lim(0.5, 0.2, 0.7), c(0.2, 0.7))
  expect_identical(lim(0.9, 0.85, 1.2), c(0.85, 1))
  expect_identical(lim(NA_real_, NA_real_, NA_real_), c(NA_real_, NA_real_))
  # Only a negative estimate is warned about.
  expect_no_warning(DMAR:::.relia_result(
    estimate = 0.5, lower = -0.1, upper = 0.8, conf_level = 0.95, N = 50,
    J = 4, coefficient = "alpha", ci_method = "bonett"))
})

test_that("replicates that are all equal give a zero-width interval, not an error", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Three items that are exact linear functions of one variable have the
  # same alpha in every resample. Through 1.0.1 boot.ci() printed "All
  # values of t are equal" and the call stopped with an error about
  # vapply(); the bootstrap distribution is a single point, whose
  # quantiles are that point.
  set.seed(1)
  x <- rnorm(40)
  X <- cbind(a = x, b = x, c = 2 * x)
  est <- manual_alpha(cov(X))
  for (m in c("percentile", "bca", "bootstrap_se")) {
    w <- testthat::capture_warnings(
      out <- utils::capture.output(
        r <- reliability_alpha(data = X, ci_method = m, B = 30, seed = 1)))
    expect_false(any(grepl("All values of t are equal", out, fixed = TRUE)))
    expect_equal(r$value[r$term %in% c("ci_lower", "ci_upper")],
                 c(est, est), tolerance = 1e-12)
    expect_true(any(grepl("gave the same estimate", w, fixed = TRUE)))
  }
  # The BCa request says that it fell back to the percentile interval.
  r <- suppressWarnings(reliability_alpha(data = X, ci_method = "bca",
                                          B = 30, seed = 1))
  expect_match(attr(r, "ci_method"), "^percentile \\(every returned replicate has the same value")
})

test_that("the logistic bootstrap interval of replicates all equal to one is that value, and the warning names the cause", {
  skip_on_cran()
  # Identical items have alpha (and KR-20) of one in every resample. The
  # logit of one is infinite, so through the first versions of 1.1.0 the
  # bootstrap_se_logistic limits were NA, under a warning that the
  # interval had no width. In 1.0.1 they were NA, or one where rounding
  # left the replicates just below one (as for KR-20 of the 0/1 item
  # below; three copies of an item split 25 and 25 gave NA). A standard
  # deviation of zero now gives the estimate at both limits, on either
  # scale, as the help says every method does.
  set.seed(1)
  x <- rnorm(50)
  for (m in c("bootstrap_se_logistic", "bootstrap_se", "percentile")) {
    w <- testthat::capture_warnings(
      r <- reliability_alpha(data = data.frame(a = x, b = x, c = x),
                             ci_method = m, B = 30, seed = 1))
    expect_equal(r$value[r$term %in% c("estimate", "ci_lower",
                                       "ci_upper")], c(1, 1, 1),
                 tolerance = 1e-12)
    expect_match(w, paste0("which happens when every item is an exact ",
                           "linear function of the same variable"),
                 all = FALSE, fixed = TRUE)
  }
  b <- stats::rbinom(50, 1, 0.5)
  k <- suppressWarnings(reliability_kr20(data.frame(a = b, b = b, c = b),
                                         ci_method = "bootstrap_se_logistic",
                                         B = 30, seed = 1))
  expect_equal(k$value[k$term %in% c("ci_lower", "ci_upper")],
               c(1, 1), tolerance = 1e-12)
  # Below one, the logistic limits of linear copies are the common value.
  r <- suppressWarnings(reliability_alpha(
    data = data.frame(a = x, b = 2 * x, c = 3 * x + 1),
    ci_method = "bootstrap_se_logistic", B = 30, seed = 1))
  est <- r$value[r$term == "estimate"]
  expect_equal(r$value[r$term %in% c("ci_lower", "ci_upper")],
               c(est, est), tolerance = 1e-12)
  # An item entered twice beside other items does not make the replicates
  # equal, and no such warning is given.
  w <- testthat::capture_warnings(
    r <- reliability_alpha(data = data.frame(a = x, b = x,
                                             c = stats::rnorm(50),
                                             e = x + stats::rnorm(50)),
                           ci_method = "percentile", B = 30, seed = 1))
  expect_false(any(grepl("gave the same estimate", w, fixed = TRUE)))
  expect_lt(r$value[r$term == "ci_lower"],
            r$value[r$term == "ci_upper"])
})

test_that("B is validated whenever a bootstrap interval is asked for", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  for (bad in list(0, 1, 2.5, -3, NA, "100", c(10, 20), Inf)) {
    expect_error(reliability_alpha(data = items_continuous,
                                   ci_method = "percentile", B = bad),
                 "'B' must be a single whole number of at least 2",
                 fixed = TRUE)
  }
  # B is not read without a bootstrap.
  expect_no_error(reliability_alpha(data = items_continuous, B = 0))
})

test_that("the analytic bootstrap records the improper rule it was given", {
  r <- reliability_alpha(data = items_continuous, ci_method = "percentile",
                         B = 50, seed = 113, improper = "drop")
  bs <- attr(r, "bootstrap")
  expect_identical(bs$improper_rule, "drop")
  expect_identical(bs$improper, 0L)
  expect_identical(bs$no_solution, 0L)
})

test_that("duplicated item names are refused on the incomplete data paths", {
  # Those paths select the items by name, and a name given to two columns
  # would select the first of them twice.
  d <- as.data.frame(items_continuous[, 1:4])
  names(d) <- c("y1", "y1", "y3", "y4")
  d$y3[1:10] <- NA
  expect_error(reliability_alpha(data = d, missing = "fiml"),
               "y1.{0,3} names more than one column")
})


# ---- The ml standard error is computed on the covariance matrix (1.1.0) ----

# The van Zyl, Neudecker, and Nel (2000) variance of alpha times N - 1,
# written out from Terry and Kelley (2012, Eqs. 7 and 8) and Kelley and
# Pornprasertmanit (2016, Eq. 20) on whatever matrix it is given.
.vz_alpha <- function(S) {
  J <- ncol(S)
  g <- sum(S)
  J^2 / (J - 1)^2 * 2 / g^3 *
    (g * (sum(diag(S %*% S)) + sum(diag(S))^2) -
       2 * sum(diag(S)) * sum(S %*% S))
}

test_that("the ml standard error of alpha is van Zyl's variance on the covariance matrix", {
  skip_on_cran()
  # 1.0.x evaluated the variance on cov2cor(S), the variance of alpha of
  # standardized items, and paired it with alpha of the raw items, as
  # MBESS's ci.reliability() does. Equal loadings .7 and error variances .1,
  # .2, 1, 2 (item standard deviations .77 to 1.58), N = 206: the standard
  # error was .0255 where alpha-hat varies with a standard deviation of
  # about .034.
  Sig <- tcrossprod(rep(.7, 4)) + diag(c(.1, .2, 1, 2))
  n <- 206
  for (m in c("ml", "ml_logistic")) {
    r <- reliability_alpha(S = Sig, N = n, ci_method = m)
    expect_equal(r$value[r$term == "se"], sqrt(.vz_alpha(Sig) / (n - 1)),
                 tolerance = 1e-12)
  }
  expect_equal(sqrt(.vz_alpha(Sig) / (n - 1)), 0.03405282, tolerance = 1e-6)
  expect_equal(sqrt(.vz_alpha(stats::cov2cor(Sig)) / (n - 1)), 0.02548508,
               tolerance = 1e-6)
  # A common rescaling of the items changes nothing.
  r9 <- reliability_alpha(S = 9 * Sig, N = n, ci_method = "ml")
  r1 <- reliability_alpha(S = Sig, N = n, ci_method = "ml")
  expect_equal(r9$value, r1$value, tolerance = 1e-12)
  # The listwise closed form now agrees with the FIML route, a lavaan delta
  # method on the maximum likelihood covariance matrix, on complete data:
  # the two differ only by the divisor, N - 1 against N. Before 1.1.0 they
  # disagreed by the ratio of the two variances above.
  skip_if_not_installed("lavaan")
  lw <- reliability_alpha(data = items_continuous, ci_method = "ml")
  fi <- reliability_alpha(data = items_continuous, ci_method = "ml",
                          missing = "fiml")
  expect_equal(fi$value[fi$term == "se"],
               lw$value[lw$term == "se"] * sqrt((N - 1) / N),
               tolerance = 1e-4)
})

test_that("the ml interval of alpha covers at the nominal rate when the item variances differ", {
  skip_on_cran()
  # Population above, N = 206, 1,000 samples. The 1.0.x interval covered
  # .86 of the time (4,000 samples); a population whose high loading items
  # carry the large variances was covered .9998 of the time, its standard
  # error about twice the standard deviation of alpha-hat.
  Sig <- tcrossprod(rep(.7, 4)) + diag(c(.1, .2, 1, 2))
  a0 <- manual_alpha(Sig)
  U <- chol(Sig)
  n <- 206
  set.seed(2026)
  covered <- replicate(1000, {
    S <- stats::cov(matrix(stats::rnorm(n * 4), n, 4) %*% U)
    r <- reliability_alpha(S = S, N = n, ci_method = "ml")
    v <- stats::setNames(r$value, r$term)
    v[["ci_lower"]] <= a0 && a0 <= v[["ci_upper"]]
  })
  expect_gt(mean(covered), .93)
  expect_lt(mean(covered), .97)
})

test_that("an item that is not numeric, or that does not vary in the rows analyzed, is refused by name before anything is computed", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford[1:100, c("t1_visual_perception", "t2_cubes",
                                     "t4_lozenges",
                                     "t6_paragraph_comprehension")]
  # Through 1.0.1 a constant item was counted among the items of the
  # analytic alpha (J(J - 2)/(J - 1)^2 times the alpha of the other items,
  # the alpha of the same composite), and the fitted coefficients stopped
  # with lavaan's "some variables have no values (only missings) or no
  # variance", which names none, after lavaan printed its table of the
  # variables.
  k <- hs
  k$t2_cubes <- 5
  msg <- paste0("t2_cubes takes a single value in the 100 complete rows, ",
                "so its variance is zero.")
  calls <- list(
    function(d, ...) reliability_alpha(d, ...),
    function(d, ...) reliability_alpha(d, estimator = "model_implied", ...),
    function(d, ...) reliability_omega(d, ...),
    function(d, ...) reliability_omega(d, denominator = "model_implied",
                                       ...),
    function(d, ...) reliability_H(data = d, ...))
  for (g in calls) {
    out <- utils::capture.output(
      err <- tryCatch(g(k), error = conditionMessage))
    expect_length(out, 0L)
    expect_match(err, msg, fixed = TRUE)
    expect_match(err, "Leave it out of 'data'.", fixed = TRUE)
  }
  # Under FIML and two-stage estimation, the rows that observe it.
  kf <- k
  kf$t1_visual_perception[1:5] <- NA
  for (m in c("fiml", "two_stage")) {
    expect_error(reliability_alpha(kf, missing = m),
                 paste0("t2_cubes takes a single value in the 100 rows that ",
                        "observe it"), fixed = TRUE)
    expect_error(reliability_omega(kf, missing = m),
                 "t2_cubes takes a single value", fixed = TRUE)
  }
  # Constant in the complete rows, but not in the rows deletion removes.
  # FIML is offered only when it could use those rows: here t2 varies
  # only in rows that miss t1, so the rows that observe both say nothing
  # about how the two vary together, and FIML would refuse the pair.
  # Listwise deletion is asked for by name, since the default analyzes
  # incomplete items by FIML.
  kv <- hs
  kv$t2_cubes[1:95] <- 5
  kv$t1_visual_perception[96:100] <- NA
  expect_error(reliability_alpha(kv, missing = "listwise"),
               paste0("t2_cubes takes a single value in the 95 complete ",
                      "rows"), fixed = TRUE)
  expect_error(reliability_alpha(kv, missing = "listwise"),
               paste0("the item varies in the rows it removes. Full ",
                      "information maximum likelihood (missing = \"fiml\") ",
                      "cannot help either: it would be refused for the pair ",
                      "t1_visual_perception with t2_cubes"), fixed = TRUE)
  kw <- hs
  kw$t2_cubes[1:95] <- 5
  kw$t1_visual_perception[96:98] <- NA
  kw$t4_lozenges[99:100] <- NA
  expect_error(reliability_alpha(kw, missing = "listwise"),
               "missing = \"fiml\" uses those rows as well", fixed = TRUE)
  # A zero variance in a covariance matrix.
  S <- stats::cov(hs)
  S[2, ] <- S[, 2] <- 0
  expect_error(reliability_alpha(S = S, N = 100),
               "The variance of t2_cubes in 'S' is zero.", fixed = TRUE)
  expect_error(reliability_omega(S = S, N = 100,
                                 denominator = "model_implied"),
               "The variance of t2_cubes in 'S' is zero.", fixed = TRUE)
  S[2, 2] <- -1
  expect_error(reliability_alpha(S = S, N = 100),
               "the entry for t2_cubes is not (-1)", fixed = TRUE)
  S[3, 3] <- 0
  expect_error(reliability_alpha(S = S, N = 100),
               "the entries for t2_cubes and t4_lozenges are not (-1, 0)",
               fixed = TRUE)
  # KR-20, which is alpha on 0/1 items, refuses it too.
  b <- as.data.frame(lapply(hs, function(x) as.numeric(x > stats::median(x))))
  b$t2_cubes <- 1
  expect_error(reliability_kr20(b),
               "t2_cubes takes a single value in the 100 complete rows",
               fixed = TRUE)
  # The first 20 pupils on the nine tests t1 to t9 scored 1 above their
  # 90th percentiles: no pupil passes four of them, and KR-20 was 0.196,
  # 0.9 of the 0.217 of the other five (the four were counted among J).
  nine <- grep("^t[1-9]_", names(holzinger_swineford), value = TRUE)[1:9]
  h9 <- as.data.frame(lapply(holzinger_swineford[nine], function(x) {
    as.numeric(x > stats::quantile(x, 0.9))
  }))[1:20, ]
  expect_error(reliability_kr20(h9),
               paste0("t1_visual_perception, t3_paper_form_board, ",
                      "t7_sentence, and t9_word_meaning each take a single ",
                      "value in the 20 complete rows, so their variances ",
                      "are zero."), fixed = TRUE)
  b$t2_cubes[1] <- 0
  b$t1_visual_perception[1] <- NA
  expect_error(reliability_kr20(b),
               paste0("Full information maximum likelihood ",
                      "(reliability_alpha() with missing = \"fiml\") cannot ",
                      "help either: it would be refused for the pair ",
                      "t1_visual_perception with t2_cubes"), fixed = TRUE)
  b$t2_cubes[2] <- 0
  b$t4_lozenges[2] <- NA
  expect_error(reliability_kr20(b),
               "reliability_alpha() with missing = \"fiml\" uses those rows",
               fixed = TRUE)
  # An item that is not numeric: stats::cov() stopped listwise deletion
  # with "is.numeric(x) || is.logical(x) is not TRUE", and FIML said that
  # the saturated model did not converge.
  ch <- hs
  ch$t2_cubes <- as.character(ch$t2_cubes)
  fa <- hs
  fa$t2_cubes <- factor(fa$t2_cubes)
  for (d in list(ch, fa)) {
    cls <- class(d$t2_cubes)[1L]
    named <- paste0("this column of 'data' is not: 't2_cubes' (", cls, ")")
    for (g in calls) expect_error(g(d), named, fixed = TRUE)
    expect_error(reliability_alpha(d, missing = "fiml"), named, fixed = TRUE)
    expect_error(reliability_omega(d, missing = "two_stage"), named,
                 fixed = TRUE)
  }
  bf <- as.data.frame(lapply(hs, function(x) {
    factor(as.numeric(x > stats::median(x)))
  }))
  expect_error(reliability_kr20(bf), "these columns of 'data' are not",
               fixed = TRUE)
  # Listwise deletion without auxiliary variables reads the columns by
  # position, so two columns may share a name there; each is checked.
  dn <- hs
  names(dn) <- c("x", "x", "y", "z")
  dn[[2]] <- as.character(dn[[2]])
  expect_error(reliability_alpha(dn), "this column of 'data' is not: 'x'",
               fixed = TRUE)
  dn[[2]] <- 5
  expect_error(reliability_alpha(dn), "x takes a single value", fixed = TRUE)
  # A logical item is read as 0 and 1, as stats::cov() reads it.
  lg <- hs
  lg$t2_cubes <- lg$t2_cubes > stats::median(lg$t2_cubes)
  expect_equal(reliability_alpha(lg)$value[1],
               DMAR:::.alpha_from_S(stats::cov(lg)), tolerance = 1e-12)
  # Categorical omega checks its ordered items itself.
  of <- data.frame(a = factor(c(1, 2, 3, 1)), b = c(1, 1, 1, 1),
                   c = c(2, 1, 2, 3))
  expect_no_error(.relia_resolve_inputs(of, S = NULL, N = NULL,
                                        allow_S = FALSE, missing = "fiml",
                                        continuous = FALSE))
})


test_that("a bootstrap replication in which an item takes a single value is set aside with the refusal", {
  # About 1.7 s of lavaan and EM refits; the KR-20 test of the same check,
  # which runs on CRAN, takes a tenth of a second.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # Regression test. An item that takes a single value is refused in the
  # data (alpha would count it among the J items), and a resample that
  # leaves out the one case on which an item differs from the rest is
  # such data. Through version 1.0.1 the analytic alpha of that resample
  # was computed with the item and kept, and the model implied refit
  # stopped in lavaan, which printed its table of the variables for every
  # such resample, and was counted as non-finite. Both now set it aside as
  # an error that quotes the refusal, and nothing is printed.
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
  bb <- boot::boot(d, function(x, i) sum(x$x5[i]), R = 30)
  lost <- sum(bb$t[, 1] == 0)
  expect_gt(lost, 0L)
  msg <- paste0("x5 takes a single value in this resample, so its variance ",
                "is zero, and an item that does not vary is refused.")
  for (est in c("analytic", "model_implied")) {
    out <- utils::capture.output(
      r <- suppressWarnings(reliability_alpha(d, estimator = est,
                                              ci_method = "percentile",
                                              B = 30, seed = 3)))
    expect_identical(out, character(0))
    bt <- attr(r, "bootstrap")
    expect_identical(bt$reasons[["error"]], lost)
    expect_identical(names(bt$error_messages), msg)
  }
  # The jackknife deletion behind the BCa acceleration of that case is
  # refused the same way.
  r <- suppressWarnings(reliability_alpha(d, ci_method = "bca", B = 30,
                                          seed = 3))
  jk <- attr(r, "bootstrap")$jackknife
  expect_identical(jk$reasons[["error"]], 1L)
  expect_identical(names(jk$error_messages), msg)
  # Under full information maximum likelihood the item is checked on the
  # rows of the resample that observe it.
  d$x1[1:5] <- NA
  r <- suppressWarnings(reliability_alpha(d, missing = "fiml",
                                          ci_method = "percentile", B = 30,
                                          seed = 3))
  bt <- attr(r, "bootstrap")
  expect_identical(bt$reasons[["error"]], lost)
  expect_identical(
    names(bt$error_messages),
    paste0("x5 takes a single value in the rows of this resample that ",
           "observe it, so its variance is zero, and an item that does not ",
           "vary is refused."))
})


test_that("a resample in which no row observes an item is refused as the data are, not as a constant item", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  # Regression test. .relia_single_valued() flags a column with no
  # observed value, and the check of a replication reported it as taking a
  # single value, so its variance zero, in the rows that observe it, and
  # tallied it with the items that do; the data are refused with "No row
  # observes".
  set.seed(1)
  n <- 60
  f <- rnorm(n)
  d <- as.data.frame(sapply(c(.8, .7, .6, .5), function(l) {
    l * f + rnorm(n, sd = sqrt(1 - l^2))
  }))
  names(d) <- paste0("x", 1:4)
  d$x4[-(1:10)] <- NA
  res <- d[rep(11:60, length.out = 60), ]
  chk <- DMAR:::.relia_checked_point(function(x) 1, missing = "fiml",
                                     items = names(d))
  expect_error(chk(res), paste0("No row of this resample observes x4, so ",
                                "neither its variance nor its covariances ",
                                "with the other items can be estimated"),
               fixed = TRUE)
  chk_aux <- DMAR:::.relia_checked_point(function(x) 1, missing = "fiml",
                                         items = names(d)[1:3], aux = "x4")
  expect_error(chk_aux(res), "covariances with the other variables",
               fixed = TRUE)
  # An item observed once still takes a single value.
  one <- res
  one$x4[1] <- 2
  expect_error(chk(one), "x4 takes a single value in the rows of this resample that observe it",
               fixed = TRUE)
  # The bootstrap of a lavaan fit in reliability_H() reads the same check.
  fit <- suppressWarnings(lavaan::cfa("f =~ x1 + x2 + x3 + x4", data = d,
                                      std.lv = TRUE, missing = "ml"))
  expect_error(DMAR:::.H_refit_fun(fit, NULL)(res),
               "No row of this resample observes x4", fixed = TRUE)
})


test_that("a composite that does not vary is refused, from data, from S, and in a bootstrap replication", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Regression test. Ipsative items, here three holzinger_swineford tests
  # each centered on the pupil's mean of the three, sum to the same total
  # in every row, so the composite has no variance and no reliability.
  # Through version 1.0.1 alpha divided the item variances by the sum of
  # the covariance matrix, which rounding error left at -7.1e-15, and
  # returned 1.42e16 with the interval [0.100, NaN]; KR-20 of two rows on
  # which as many items are passed in one as in the other returned -Inf
  # with the limits 0 and -Inf.
  x <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t4_lozenges")]
  ip <- as.data.frame(as.matrix(x) - rowMeans(x))
  expect_lt(abs(sum(stats::cov(ip))), 1e-12)
  msg <- paste0("The composite, the sum of the items, does not vary: its ",
                "variance, the sum of the entries of the covariance ",
                "matrix, is zero up to rounding error, as when the items ",
                "are ipsative (their scores sum to the same total in every ",
                "row), so the reliability of the composite is not defined ",
                "(alpha, J/(J - 1) times one minus the ratio of the item ",
                "variances to that variance, has no finite value, and the ",
                "single-factor model cannot be fit to the singular ",
                "covariance matrix).")
  ipm <- ip
  ipm[1:20, 1] <- NA
  calls <- list(
    function() reliability_alpha(ip),
    function() reliability_alpha(ip, ci_method = "feldt"),
    function() reliability_alpha(S = stats::cov(ip), N = 301),
    function() reliability_alpha(ip, estimator = "model_implied"),
    function() reliability_omega(ip),
    function() reliability_omega(ip, denominator = "model_implied"))
  for (f in calls) {
    expect_identical(
      tryCatch(suppressWarnings(suppressMessages(f())),
               error = conditionMessage), msg)
  }
  k <- data.frame(a = c(1, 0), b = c(0, 1), c = c(1, 0), e = c(0, 1))
  expect_error(reliability_kr20(k), msg, fixed = TRUE)
  # A replication whose rows all have the same total is refused the same
  # way, and counted.
  k4 <- data.frame(a = c(1, 0, 1, 1), b = c(0, 1, 1, 1), c = c(1, 0, 0, 1),
                   e = c(0, 1, 0, 0))
  r <- suppressWarnings(reliability_kr20(k4, ci_method = "percentile",
                                         B = 50, seed = 1))
  em <- attr(r, "bootstrap")$error_messages
  expect_true(msg %in% names(em))
  # Under FIML the rows that observe every item are checked.
  expect_error(reliability_alpha(ipm, missing = "fiml"),
               "does not vary in the 281 rows that observe every item",
               fixed = TRUE)
})


test_that("a supplied S whose composite variance is not above zero is refused as not a covariance matrix when it is indefinite", {
  # Regression test. Every S whose sum was not above the tolerance was
  # refused as a composite that does not vary, its variance zero up to
  # rounding error and its matrix singular, which is false of an
  # indefinite matrix: three items that correlate -.9 have the sum -2.4
  # and the eigenvalues 1.9, 1.9, and -0.8 (1.0.1 returned alpha 3.38 with
  # the limits [4.31, 1]).
  S <- matrix(-.9, 3, 3)
  diag(S) <- 1
  dimnames(S) <- list(letters[1:3], letters[1:3])
  msg <- paste0("'S' is not a covariance matrix: its smallest eigenvalue after standardizing item units is ",
                "-0.8, so some weighted sum of the items would have a ",
                "negative variance, which no population has, and the ",
                "composite, the sum of the items, would have the variance ",
                "-2.4 (the sum of the entries of 'S').")
  expect_identical(tryCatch(reliability_alpha(S = S, N = 100,
                                              ci_method = "feldt"),
                            error = conditionMessage), msg)
  expect_identical(tryCatch(reliability_omega(S = S, N = 100),
                            error = conditionMessage), msg)
  # An indefinite matrix whose sum is zero to rounding is named the same
  # way, without a composite variance.
  S0 <- matrix(c(1, -.9, -.9, -.9, 1, .3, -.9, .3, 1), 3, 3)
  expect_lt(abs(sum(S0)), 1e-12)
  err <- tryCatch(reliability_alpha(S = S0, N = 100, ci_method = "feldt"),
                  error = conditionMessage)
  expect_match(err, "'S' is not a covariance matrix: its smallest eigenvalue after standardizing item units is -0.132,",
               fixed = TRUE)
  expect_no_match(err, "would have the variance", fixed = TRUE)
  # The covariance matrix of ipsative items keeps the composite's message.
  ip <- matrix(-.5, 3, 3)
  diag(ip) <- 1
  expect_error(reliability_alpha(S = ip, N = 100, ci_method = "feldt"),
               "The composite, the sum of the items, does not vary",
               fixed = TRUE)
})


test_that("a composite that does not vary is refused under FIML and two-stage estimation for every coefficient of the composite", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  # Regression test. The refusal was applied only under listwise deletion
  # and to S, and under FIML and two-stage estimation only where a path
  # happened to compute alpha from the EM estimate, whose composite
  # variance the EM leaves at its convergence error rather than at
  # rounding error. With ten values of the first of four ipsative tests
  # missing, the model implied alpha under FIML was 1.5e-13 with the
  # interval [0, 9.2e-8] (as in 1.0.1), under two-stage estimation
  # 1.9e-12, and two-stage omega reported a model that did not converge;
  # with a tenth of the values missing across the items the analytic
  # alpha under FIML was of the order of -5e7.
  x <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t4_lozenges", "t6_paragraph_comprehension")]
  ip <- as.data.frame(as.matrix(x) - rowMeans(x))
  ipm <- ip
  ipm[1:10, 1] <- NA
  rows <- function(n, label) {
    paste0("The composite, the sum of the items, does not vary in the ", n,
           " rows that observe every item: its variance there is zero up ",
           "to rounding error, as when the items are ipsative (their scores ",
           "sum to the same total in every row), so the reliability of the ",
           "composite is not defined (alpha, J/(J - 1) times one minus the ",
           "ratio of the item variances to that variance, has no finite ",
           "value), and the likelihood of the saturated model that ", label,
           " estimates first grows without bound as its estimate of that ",
           "variance goes to zero, so the method has no estimate to give ",
           "either.")
  }
  fiml <- "full information maximum likelihood"
  ts <- "two-stage maximum likelihood"
  calls <- list(
    list(function() reliability_alpha(ipm, missing = "fiml"), 291, fiml),
    list(function() reliability_alpha(ipm, missing = "two_stage"), 291, ts),
    list(function() reliability_alpha(ipm, missing = "fiml",
                                      estimator = "model_implied",
                                      ci_method = "ml"), 291, fiml),
    list(function() reliability_alpha(ipm, missing = "two_stage",
                                      estimator = "model_implied",
                                      ci_method = "ml"), 291, ts),
    list(function() reliability_omega(ipm, missing = "fiml"), 291, fiml),
    list(function() reliability_omega(ipm, missing = "fiml",
                                      denominator = "model_implied"),
         291, fiml),
    list(function() reliability_omega(ipm, missing = "two_stage"), 291, ts),
    list(function() reliability_omega(ipm, missing = "two_stage",
                                      denominator = "model_implied"),
         291, ts),
    list(function() reliability_alpha(ip, missing = "fiml",
                                      estimator = "model_implied",
                                      ci_method = "ml"), 301, fiml))
  for (cl in calls) {
    expect_identical(tryCatch(suppressWarnings(suppressMessages(cl[[1L]]())),
                              error = conditionMessage),
                     rows(cl[[2L]], cl[[3L]]))
  }
  # A tenth of the values missing across the items.
  set.seed(1)
  ipx <- ip
  ipx[matrix(stats::runif(301 * 4) < .10, 301, 4)] <- NA
  for (m in c("fiml", "two_stage")) {
    expect_error(suppressWarnings(reliability_alpha(ipx, missing = m)),
                 "does not vary in the 198 rows that observe every item",
                 fixed = TRUE)
  }
  # A bootstrap replication is checked the same way.
  chk <- DMAR:::.relia_checked_point(function(d) 1, missing = "fiml",
                                     items = names(ip), composite = TRUE)
  expect_error(chk(ip), rows(301, fiml), fixed = TRUE)
  # An EM estimate is held to a tolerance of 1e-6 times the sum of the item
  # variances, not to rounding error: a composite variance of 2e-8 times
  # that sum passes the second and fails the first.
  S0 <- stats::cov(ip)
  S1 <- S0 + 2e-8 * sum(diag(S0)) / 16
  expect_true(DMAR:::.relia_composite_varies(S1))
  expect_false(DMAR:::.relia_composite_varies(S1, em = TRUE))
  expect_error(DMAR:::.alpha_from_S(S1, em = TRUE),
               "is below a millionth of the sum of the item variances",
               fixed = TRUE)
})


test_that("the default 'missing' is FIML for incomplete items and keeps Bonett's interval for complete data", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # DMAR 1.0.x deleted incomplete rows listwise by default, with Bonett's
  # interval. Since 1.1.0 the default estimates alpha from incomplete
  # items by FIML, with the "ml" interval (Bonett's needs complete cases);
  # complete data and a covariance matrix are analyzed as they stand, with
  # the interval they had.
  d <- as.data.frame(items_continuous)
  d_na <- d
  d_na$y3[seq(1, 120, by = 3)] <- NA
  a <- reliability_alpha(d)
  expect_identical(attr(a, "missing"), "listwise")
  expect_identical(attr(a, "ci_method"), "bonett")
  expect_identical(a$value, reliability_alpha(d, missing = "listwise")$value)
  s <- reliability_alpha(S = S_continuous, N = N)
  expect_identical(attr(s, "missing"), "listwise")
  expect_identical(attr(s, "ci_method"), "bonett")
  expect_no_warning(x <- reliability_alpha(d_na))
  expect_identical(attr(x, "missing"), "fiml")
  expect_identical(attr(x, "ci_method"), "ml")
  expect_equal(x$value[x$term == "N"], nrow(d_na))
  expect_equal(x$value, reliability_alpha(d_na, missing = "fiml")$value)
  # The model implied estimator follows the same rule.
  m <- reliability_alpha(d_na, estimator = "model_implied")
  expect_identical(attr(m, "missing"), "fiml")
  expect_identical(attr(reliability_alpha(d, estimator = "model_implied"),
                        "missing"), "listwise")
  # Listwise deletion asked for by name is the DMAR 1.0.x analysis.
  expect_no_warning(l <- reliability_alpha(d_na, missing = "listwise"))
  expect_identical(attr(l, "ci_method"), "bonett")
  expect_equal(l$value[l$term == "N"], sum(stats::complete.cases(d_na)))
  # The robust interval is refused with the analytic estimator.
  expect_error(reliability_alpha(d_na, ci_method = "mlr"),
               "is not available with estimator = \"analytic\"",
               fixed = TRUE)
  # An interval FIML cannot give is computed as before with complete data
  # and, with values missing, stops and names the choices for the
  # estimator in use.
  for (m in c("bonett", "feldt", "fisher", "hakstian_whalen", "adf")) {
    expect_identical(attr(reliability_alpha(d, ci_method = m), "missing"),
                     "listwise")
    expect_error(reliability_alpha(d_na, ci_method = m),
                 paste0("not available with ci_method = \"", m, "\", and ",
                        "the items have missing values"), fixed = TRUE)
  }
  expect_error(reliability_alpha(d_na, ci_method = "bonett"),
               paste0("Use missing = \"listwise\" to analyze the complete ",
                      "cases, or keep FIML with ci_method = \"ml\", ",
                      "\"ml_logistic\", or a bootstrap method."),
               fixed = TRUE)
  expect_error(reliability_alpha(d_na, estimator = "model_implied",
                                 ci_method = "likelihood"),
               paste0("keep FIML with ci_method = \"ml\", \"ml_logistic\", ",
                      "\"mlr\", \"mlr_logistic\", or a bootstrap method."),
               fixed = TRUE)
  expect_error(reliability_alpha(d_na, ci_method = "bonett",
                                 missing = "fiml"),
               "is a complete-data closed form", fixed = TRUE)
  # An interval the estimator does not offer is said first, so that the
  # remedy named is never refused next: the FIML refusal would have
  # offered missing = "listwise", which then stops on the estimator.
  for (miss in list(NULL, "fiml", "listwise")) {
    expect_error(do.call(reliability_alpha,
                         c(list(d_na, ci_method = "likelihood"),
                           if (!is.null(miss)) list(missing = miss))),
                 "is not available with estimator = \"analytic\"",
                 fixed = TRUE)
    expect_error(do.call(reliability_alpha,
                         c(list(d_na, ci_method = "feldt",
                                estimator = "model_implied"),
                           if (!is.null(miss)) list(missing = miss))),
                 "is not available with estimator = \"model_implied\"",
                 fixed = TRUE)
  }
  # A row with every item missing does not make the data incomplete:
  # listwise deletion and FIML both drop it, so complete data with a blank
  # row keep the complete data analysis and Bonett's interval, with no
  # warning.
  blank <- rbind(d, NA)
  expect_no_warning(b <- reliability_alpha(blank))
  expect_identical(attr(b, "missing"), "listwise")
  expect_identical(attr(b, "ci_method"), "bonett")
  expect_identical(b$value, a$value)
  expect_identical(attr(reliability_H(data = blank), "missing"), "listwise")
  # Auxiliary variables imply FIML whether or not a value is missing.
  z <- data.frame(d, z = d$y1 + stats::rnorm(nrow(d)))
  expect_no_warning(za <- reliability_alpha(z, aux = "z"))
  expect_identical(attr(za, "missing"), "fiml")
  # reliability() forwards an unsupplied 'missing' unsupplied.
  expect_no_warning(r <- reliability(d_na, type = "alpha"))
  expect_equal(r$value, x$value)
  # The default on incomplete items is the classical coefficient of the
  # saturated FIML covariance matrix, with its delta method interval:
  # lavaan's own FIML fit of the saturated model, with alpha defined in
  # it, gives the same estimate, standard error, and limits.
  k <- ncol(d_na)
  ij <- which(upper.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  s <- paste0("s", ij[, 1L], "_", ij[, 2L])
  v <- s[ij[, 1L] == ij[, 2L]]
  cv <- s[ij[, 1L] != ij[, 2L]]
  syntax <- paste0(
    paste0(names(d_na)[ij[, 1L]], " ~~ ", s, "*", names(d_na)[ij[, 2L]],
           collapse = "\n"),
    "\nalpha := ", k, " / ", k - 1L, " * (1 - (", paste(v, collapse = " + "),
    ") / (", paste(v, collapse = " + "), " + 2 * (",
    paste(cv, collapse = " + "), ")))")
  fit <- lavaan::sem(syntax, data = d_na, missing = "ml")
  pe <- lavaan::parameterEstimates(fit)
  al <- pe[pe$label == "alpha", ]
  expect_equal(x$value[match(c("estimate", "se", "ci_lower",
                               "ci_upper", "N"), x$term)],
               c(al$est, al$se, al$ci.lower, al$ci.upper,
                 lavaan::lavInspect(fit, "ntotal")),
               tolerance = 1e-6)
})
