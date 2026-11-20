test_that("equivalence_c() returns documented rows", {
  res <- equivalence_c(contrast = -5.28, se = 2.49, df_error = 399, contrast_upper = 5)
  # The 17 rows of DMAR 1.0.x keep their positions; nonsuperior, added in
  # 1.1.0, is the 18th, so code that reads the rows by position still finds
  # every earlier row where it was.
  expect_identical(res$term,
                   c("contrast", "se", "df",
                     "t_lower", "t_upper",
                     "p_lower", "p_upper", "p_tost", "p_noninferiority",
                     "ci_lower", "ci_upper",
                     "contrast_lower", "contrast_upper",
                     "equivalent", "noninferior", "superior", "inferior",
                     "nonsuperior"))
  expect_type(res$value, "double")
  expect_true(attr(res, "verdict") %in%
                c("Equivalent", "Superior", "Inferior",
                  "Noninferior only", "Not relevantly superior only",
                  "Inconclusive"))
  expect_identical(.EQUIVALENCE_VERDICTS,
                   c("Equivalent", "Superior", "Inferior",
                     "Noninferior only", "Not relevantly superior only",
                     "Inconclusive"))
})

test_that("equivalence_c() reproduces a worked contrast to published precision", {
  # A two-group contrast sharing a five-group pooled error term:
  # psi = -5.28, SE = 2.49, nu = 399, bounds of 5. The 90% CI is
  # [-9.38, -1.17]: it straddles the lower bound and lies below the upper
  # one, so the contrast is not relevantly superior and nothing more. DMAR
  # 1.0.x, which had no nonsuperior flag, labeled it "Inconclusive".
  res <- equivalence_c(contrast = -5.28, se = 2.49, df_error = 399, contrast_upper = 5)
  expect_equal(res$value[res$term == "ci_lower"], -9.386, tolerance = 1e-3)
  expect_equal(res$value[res$term == "ci_upper"], -1.174, tolerance = 1e-3)
  expect_equal(res$value[res$term == "p_noninferiority"], 0.5447,
               tolerance = 1e-3)
  expect_identical(attr(res, "verdict"), "Not relevantly superior only")
  expect_equal(res$value[res$term == "equivalent"], 0)
  expect_equal(res$value[res$term == "noninferior"], 0)
  expect_equal(res$value[res$term == "nonsuperior"], 1)
})

test_that("equivalence_c() flags a switch that is not warranted", {
  # A retrieval evaluation contrast of -.004 (SE .002, 200 df) against a
  # noninferiority margin of .005 and a switching threshold of .0172. The
  # 90% limits are -.00731 and -.000695; the upper limit is below the
  # threshold, so the improvement does not justify a switch, while the
  # lower limit leaves noninferiority open.
  e <- equivalence_c(contrast = -.004, se = .002, df_error = 200,
                     contrast_lower = .005, contrast_upper = .0172)
  v <- setNames(e$value, e$term)
  t_crit <- qt(.95, 200)
  expect_equal(v[["ci_lower"]], -.004 - t_crit * .002, tolerance = 1e-12)
  expect_equal(v[["ci_upper"]], -.004 + t_crit * .002, tolerance = 1e-12)
  expect_equal(v[["ci_upper"]], -0.000695, tolerance = 1e-3)
  expect_equal(v[["nonsuperior"]], 1)
  expect_equal(v[c("equivalent", "noninferior", "superior", "inferior")],
               c(equivalent = 0, noninferior = 0, superior = 0, inferior = 0))
  expect_identical(attr(e, "verdict"), "Not relevantly superior only")
})

test_that("equivalence_c()'s one-sided flags are the decisions of the two one-sided tests", {
  # noninferior is 1 exactly when p_lower < alpha_level, and nonsuperior
  # exactly when p_upper < alpha_level; superior and inferior are the
  # one-sided tests against the opposite bounds. Checked over a grid that
  # crosses every bound, at two alpha levels and asymmetric bounds.
  for (alpha in c(.05, .025)) {
    for (psi in seq(-9, 9, by = 0.75)) {
      for (se in c(0.5, 1.5, 4)) {
        r <- equivalence_c(contrast = psi, se = se, df_error = 40,
                           contrast_lower = 3, contrast_upper = 5,
                           alpha_level = alpha)
        v <- setNames(r$value, r$term)
        expect_identical(v[["noninferior"]] == 1, v[["p_lower"]] < alpha)
        expect_identical(v[["nonsuperior"]] == 1, v[["p_upper"]] < alpha)
        expect_identical(v[["superior"]] == 1,
                         pt((psi - 5) / se, 40, lower.tail = FALSE) < alpha)
        expect_identical(v[["inferior"]] == 1,
                         pt((psi + 3) / se, 40) < alpha)
        expect_identical(v[["equivalent"]],
                         v[["noninferior"]] * v[["nonsuperior"]])
        # All five flags are 0 exactly when the verdict is inconclusive.
        flags <- v[c("equivalent", "noninferior", "superior", "inferior",
                     "nonsuperior")]
        expect_identical(all(flags == 0),
                         attr(r, "verdict") == "Inconclusive")
      }
    }
  }
})

test_that("equivalence_c() summary-statistic interface matches the direct interface", {
  means <- c(70.40, 55.61, 51.91, 65.66, 65.12)
  n     <- c(113, 74, 76, 80, 61)
  w     <- c(-1, 0, 0, 0, 1)
  r1 <- equivalence_c(means = means, s_anova = 15.673, c_weights = w, n_per_group = n,
               contrast_upper = 5)
  se <- 15.673 * sqrt(sum(w^2 / n))
  r2 <- equivalence_c(contrast = sum(w * means), se = se, df_error = sum(n) - 5,
               contrast_upper = 5)
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
                c_weights = c(-1, 1, 0), n_per_group = c(20, 25, 30), contrast_upper = 2)
  expect_equal(res$value[res$term == "p_tost"], 0.3788590420812228,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "p_noninferiority"], 0.003602147664032475,
               tolerance = 1e-12)
  # emmeans::test(..., delta = 2, side = "nonsuperiority", adjust = "none")
  # gives the upper one-sided test, the one behind the nonsuperior flag
  # (emmeans 2.0.4, 2026-11-05; live comparison in
  # tools/oracle_snippets/emmeans.R, run by tools/oracle_checks.R).
  expect_equal(res$value[res$term == "p_upper"], 0.3788590420812226,
               tolerance = 1e-12)
})

test_that("equivalence_c() produces each of the six verdicts", {
  v <- function(est, lo_hi) {
    # Solve for se so the 90% CI is exactly [lo, hi] around est.
    se <- (lo_hi[2] - lo_hi[1]) / (2 * qt(0.95, 100))
    attr(equivalence_c(contrast = est, se = se, df_error = 100, contrast_upper = 5),
         "verdict")
  }
  expect_identical(v(-1.0, c(-3.2, 1.2)),   "Equivalent")
  expect_identical(v( 8.0, c( 5.5, 10.5)),  "Superior")
  expect_identical(v(-8.0, c(-10.5, -5.5)), "Inferior")
  expect_identical(v( 3.5, c(-1.4, 8.4)),   "Noninferior only")
  # An interval below the upper bound that straddles the lower one was
  # "Inconclusive" through DMAR 1.0.x.
  expect_identical(v(-1.5, c(-6.6, 3.6)),   "Not relevantly superior only")
  expect_identical(v( 0.5, c(-6.0, 7.0)),   "Inconclusive")
})

test_that("the verdict helper is strict at the bounds and vectorized", {
  # A limit exactly on a bound establishes nothing about that bound, as
  # the one-sided tests' strict inequalities say.
  expect_identical(.equivalence_verdict(-5, 4, 5, 5), "Not relevantly superior only")
  expect_identical(.equivalence_verdict(-4, 5, 5, 5), "Noninferior only")
  expect_identical(.equivalence_verdict(-5, 5, 5, 5), "Inconclusive")
  expect_identical(.equivalence_verdict(5, 7, 5, 5), "Noninferior only")
  expect_identical(.equivalence_verdict(-7, -5, 5, 5), "Not relevantly superior only")
  # Asymmetric bounds: the region is (-2, 8).
  lo <- c(-1, 9, -12, -1, -3, -3)
  hi <- c( 7, 12, -3, 10,  7, 9)
  expect_identical(.equivalence_verdict(lo, hi, 2, 8),
                   c("Equivalent", "Superior", "Inferior",
                     "Noninferior only", "Not relevantly superior only",
                     "Inconclusive"))
  # The bounds are named for what they bound (lower_bound and upper_bound,
  # delta_lower and delta_upper through the first 1.1.0 builds); the callers
  # pass them by position, and by name the verdict is the same.
  expect_identical(names(formals(.equivalence_verdict)),
                   c("lower", "upper", "lower_bound", "upper_bound"))
  expect_identical(.equivalence_verdict(lower = lo, upper = hi,
                                        lower_bound = 2, upper_bound = 8),
                   .equivalence_verdict(lo, hi, 2, 8))
})

test_that("equivalence_c() benchmark interface removes the constant's variability", {
  res <- equivalence_c(means = c(70.40, 60), s_anova = 15.673,
                c_weights = c(1, 0), n_per_group = c(113, 50),
                benchmark = 68, contrast_upper = 5)
  expect_equal(res$value[res$term == "contrast"], 2.40, tolerance = 1e-10)
  expect_equal(res$value[res$term == "se"], 15.673 * sqrt(1 / 113),
               tolerance = 1e-10)
})

test_that("equivalence_c() enforces the raw-scale weight condition", {
  expect_error(equivalence_c(means = c(1, 2), s_anova = 1, c_weights = c(2, -2),
                      n_per_group = 10, contrast_upper = 1),
               "positive weights must sum to 1")
  expect_error(equivalence_c(means = c(1, 2), s_anova = 1, c_weights = c(1, -0.5),
                      n_per_group = 10, contrast_upper = 1),
               "sum of the contrast weights")
})

test_that("equivalence_c() rejects mixed and incomplete interfaces", {
  expect_error(equivalence_c(means = c(1, 2), s_anova = 1, c_weights = c(1, -1),
                      n_per_group = 10, contrast = 0.5, se = 1, df_error = 10,
                      contrast_upper = 1),
               "not both")
  expect_error(equivalence_c(contrast = 0.5, se = 1, contrast_upper = 1),
               "df_error")
  expect_error(equivalence_c(contrast = 0.5, se = 1, df_error = 10),
               "contrast_upper")
})

test_that("equivalence_c() rejects non-finite inputs (MEDIUM-07)", {
  # An infinite estimate is not an admissible contrast; it previously produced a
  # "Superior" verdict with infinite interval endpoints.
  expect_error(equivalence_c(contrast = Inf, se = 1, df_error = 20, contrast_upper = 2),
               "finite")
  expect_error(equivalence_c(contrast = 0.5, se = Inf, df_error = 20, contrast_upper = 2),
               "finite")
  expect_error(equivalence_c(contrast = 0.5, se = 1, df_error = Inf, contrast_upper = 2),
               "finite")
  expect_error(equivalence_c(contrast = 0.5, se = 1, df_error = 20, contrast_upper = Inf),
               "finite")
  # Ordinary finite inputs are unaffected.
  r <- equivalence_c(contrast = 0.5, se = 1, df_error = 20, contrast_upper = 2)
  expect_true(is.finite(r$value[r$term == "ci_lower"]))
})

test_that("results_sentence() reports the TOST p-value, on the table and on its tidy()", {
  # The 1.1.0 review: the sentence reported p_lower (3.9e-07), one of the
  # two one-sided tests, where the TOST p-value is .039. Example 3 of the
  # help page, the first group against a benchmark of 68.
  y <- equivalence_c(means = c(70.40, 55.61, 51.91, 65.66, 65.12),
                     s_anova = 15.67, c_weights = c(1, 0, 0, 0, 0),
                     n_per_group = c(113, 74, 76, 80, 61), benchmark = 68,
                     contrast_upper = 5)
  p <- setNames(y$value, y$term)
  expect_lt(p[["p_lower"]], 1e-4)
  want <- "contrast = 2.40, 90% CI [-0.03, 4.83], TOST p = 0.0393"
  expect_identical(want, paste0("contrast = 2.40, 90% CI [-0.03, 4.83], ",
                                "TOST p = ", format_p(p[["p_tost"]])))
  expect_identical(results_sentence(y), want)
  expect_identical(results_sentence(generics::tidy(y)), want)
})
