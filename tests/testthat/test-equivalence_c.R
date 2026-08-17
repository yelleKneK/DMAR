test_that("equivalence_c() returns documented rows", {
  res <- equivalence_c(psi_hat = -5.28, se = 2.49, df_error = 399, delta_upper = 5)
  expect_setequal(res$term,
                  c("psi_hat", "se", "df",
                    "t_lower", "t_upper",
                    "p_lower", "p_upper", "p_tost", "p_noninferiority",
                    "lower_limit", "upper_limit",
                    "delta_lower", "delta_upper",
                    "equivalent", "noninferior", "superior", "inferior"))
  expect_true(attr(res, "verdict") %in%
                c("Equivalent", "Superior", "Inferior",
                  "Noninferior only", "Inconclusive"))
})

test_that("equivalence_c() reproduces a worked contrast to published precision", {
  # A two-group contrast sharing a five-group pooled error term:
  # psi_hat = -5.28, SE = 2.49, nu = 399, bounds of 5. The 90% CI is
  # [-9.38, -1.17] and the interval straddles the lower bound.
  res <- equivalence_c(psi_hat = -5.28, se = 2.49, df_error = 399, delta_upper = 5)
  expect_equal(res$value[res$term == "lower_limit"], -9.386, tolerance = 1e-3)
  expect_equal(res$value[res$term == "upper_limit"], -1.174, tolerance = 1e-3)
  expect_equal(res$value[res$term == "p_noninferiority"], 0.5447,
               tolerance = 1e-3)
  expect_identical(attr(res, "verdict"), "Inconclusive")
  expect_equal(res$value[res$term == "equivalent"], 0)
  expect_equal(res$value[res$term == "noninferior"], 0)
})

test_that("equivalence_c() summary-statistic interface matches the direct interface", {
  means <- c(70.40, 55.61, 51.91, 65.66, 65.12)
  n     <- c(113, 74, 76, 80, 61)
  w     <- c(-1, 0, 0, 0, 1)
  r1 <- equivalence_c(means = means, s_anova = 15.673, c_weights = w, n = n,
               delta_upper = 5)
  se <- 15.673 * sqrt(sum(w^2 / n))
  r2 <- equivalence_c(psi_hat = sum(w * means), se = se, df_error = sum(n) - 5,
               delta_upper = 5)
  expect_equal(r1$value[r1$term == "p_tost"],
               r2$value[r2$term == "p_tost"], tolerance = 1e-12)
  expect_equal(r1$value[r1$term == "df"], 399)
})

test_that("equivalence_c() agrees with pinned emmeans p-values to machine precision", {
  set.seed(113)
  d <- data.frame(
    g = factor(rep(c("a", "b", "c"), times = c(20, 25, 30))),
    y = rnorm(75, mean = rep(c(10, 10.5, 12), times = c(20, 25, 30)), sd = 4)
  )
  fit <- lm(y ~ g, data = d)
  # The means below are the emmeans reference-grid means for this fit, and the
  # two anchor p-values are the equivalence and noninferiority tests of the
  # b-versus-a contrast with delta = 2 and adjust = "none". Pinned from
  # emmeans::emmeans and emmeans::test (emmeans 2.0.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  res <- equivalence_c(means = c(10.07931054038358, 11.67660032745882,
                          11.80114880502362),
                s_anova = summary(fit)$sigma,
                c_weights = c(-1, 1, 0), n = c(20, 25, 30), delta_upper = 2)
  expect_equal(res$value[res$term == "p_tost"], 0.3788590420812228,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "p_noninferiority"], 0.003602147664032475,
               tolerance = 1e-12)
})

test_that("equivalence_c() produces each of the five verdicts", {
  v <- function(est, lo_hi) {
    # Solve for se so the 90% CI is exactly [lo, hi] around est.
    se <- (lo_hi[2] - lo_hi[1]) / (2 * qt(0.95, 100))
    attr(equivalence_c(psi_hat = est, se = se, df_error = 100, delta_upper = 5),
         "verdict")
  }
  expect_identical(v(-1.0, c(-3.2, 1.2)),   "Equivalent")
  expect_identical(v( 8.0, c( 5.5, 10.5)),  "Superior")
  expect_identical(v(-8.0, c(-10.5, -5.5)), "Inferior")
  expect_identical(v( 3.5, c(-1.4, 8.4)),   "Noninferior only")
  expect_identical(v(-1.5, c(-6.6, 3.6)),   "Inconclusive")
})

test_that("equivalence_c() benchmark interface removes the constant's variability", {
  res <- equivalence_c(means = c(70.40, 60), s_anova = 15.673,
                c_weights = c(1, 0), n = c(113, 50),
                benchmark = 68, delta_upper = 5)
  expect_equal(res$value[res$term == "psi_hat"], 2.40, tolerance = 1e-10)
  expect_equal(res$value[res$term == "se"], 15.673 * sqrt(1 / 113),
               tolerance = 1e-10)
})

test_that("equivalence_c() enforces the raw-scale weight condition", {
  expect_error(equivalence_c(means = c(1, 2), s_anova = 1, c_weights = c(2, -2),
                      n = 10, delta_upper = 1),
               "positive weights must sum to 1")
  expect_error(equivalence_c(means = c(1, 2), s_anova = 1, c_weights = c(1, -0.5),
                      n = 10, delta_upper = 1),
               "sum of the contrast weights")
})

test_that("equivalence_c() rejects mixed and incomplete interfaces", {
  expect_error(equivalence_c(means = c(1, 2), s_anova = 1, c_weights = c(1, -1),
                      n = 10, psi_hat = 0.5, se = 1, df_error = 10,
                      delta_upper = 1),
               "not both")
  expect_error(equivalence_c(psi_hat = 0.5, se = 1, delta_upper = 1),
               "df_error")
  expect_error(equivalence_c(psi_hat = 0.5, se = 1, df_error = 10),
               "delta_upper")
})

test_that("equivalence_c() rejects non-finite inputs (MEDIUM-07)", {
  # An infinite estimate is not an admissible contrast; it previously produced a
  # "Superior" verdict with infinite interval endpoints.
  expect_error(equivalence_c(psi_hat = Inf, se = 1, df_error = 20, delta_upper = 2),
               "finite")
  expect_error(equivalence_c(psi_hat = 0.5, se = Inf, df_error = 20, delta_upper = 2),
               "finite")
  expect_error(equivalence_c(psi_hat = 0.5, se = 1, df_error = Inf, delta_upper = 2),
               "finite")
  expect_error(equivalence_c(psi_hat = 0.5, se = 1, df_error = 20, delta_upper = Inf),
               "finite")
  # Ordinary finite inputs are unaffected.
  r <- equivalence_c(psi_hat = 0.5, se = 1, df_error = 20, delta_upper = 2)
  expect_true(is.finite(r$value[r$term == "lower_limit"]))
})
