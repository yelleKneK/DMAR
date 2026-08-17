test_that("ancova() returns expected term rows", {
  set.seed(113)
  n <- 30; grp <- rep(c("A","B"), each = 15)
  x <- rnorm(n, 50, 10); y <- 0.6 * x + 5 * (grp == "B") + rnorm(n, 0, 4)
  d <- data.frame(y, group = grp, x)
  res <- ancova(d, "y", "group", "x")
  expect_true("F_value" %in% res$term)
  expect_true("eta_squared_partial" %in% res$term)
  expect_true("omega_squared_partial" %in% res$term)
  expect_true("adjusted_mean[A]" %in% res$term)
  expect_true("adjusted_mean[B]" %in% res$term)
  expect_true("F_homogeneity_of_regression" %in% res$term)
})

test_that("ancova() omnibus F is the covariate-adjusted treatment effect", {
  set.seed(113)
  n <- 60
  grp <- factor(rep(c("ctrl", "drug", "placebo"), each = 20))
  x   <- rnorm(n, 50, 10)
  y   <- 0.5 * x + 3 * (grp == "drug") + rnorm(n, 0, 3)
  d <- data.frame(y, grp, x)
  res <- ancova(d, "y", "grp", "x")
  # The ANCOVA test is the covariate-adjusted treatment effect: enter the
  # covariate first and the treatment last so its sequential SS is adjusted.
  ref_adj <- stats::anova(stats::lm(y ~ x + grp, data = d))
  expect_equal(res$value[res$term == "F_value"],
               ref_adj["grp", "F value"], tolerance = 1e-10)
  # And it is NOT the treatment-first (unadjusted) F.
  ref_unadj <- stats::anova(stats::lm(y ~ grp + x, data = d))
  expect_false(isTRUE(all.equal(ref_adj["grp", "F value"],
                                ref_unadj["grp", "F value"])))
})

test_that("ancova() rejects missing columns", {
  d <- data.frame(y = 1:5, grp = factor(1:5), x = rnorm(5))
  expect_error(ancova(d, "missing", "grp", "x"), "column name")
})

test_that("ancova() reports Type III SS and matches car::Anova for the adjusted F", {
  skip_if_not_installed("car")
  data(pygmalion, package = "DMAR")
  res <- ancova(outcome = "iq_8", treatment = "treatment",
                covariates = "iq_pre", data = pygmalion)
  expect_equal(res$value[res$term == "sum_of_squares_type"], 3)
  # every reported effect size carries a confidence interval
  expect_true(all(c("eta_squared_partial_lower", "eta_squared_partial_upper",
                    "omega_squared_partial_lower", "omega_squared_partial_upper")
                  %in% res$term))

  d <- pygmalion[stats::complete.cases(
    pygmalion[, c("iq_8", "treatment", "iq_pre")]), ]
  d$treatment <- factor(d$treatment)
  # Type II needs no special contrasts; for the additive one-way ANCOVA it
  # equals Type III, and both equal DMAR's covariate-adjusted F.
  fit2 <- stats::lm(iq_8 ~ treatment + iq_pre, data = d)
  # suppressWarnings: car's internals partially match summary(corr =) and
  # warn under warnPartialMatchArgs; not a property under test.
  f_car2 <- suppressWarnings(car::Anova(fit2, type = 2))["treatment", "F value"]
  fit3 <- stats::lm(iq_8 ~ treatment + iq_pre, data = d,
                    contrasts = list(treatment = stats::contr.sum))
  f_car3 <- suppressWarnings(car::Anova(fit3, type = 3))["treatment", "F value"]
  expect_equal(res$value[res$term == "F_value"], f_car2, tolerance = 1e-8)
  expect_equal(res$value[res$term == "F_value"], f_car3, tolerance = 1e-8)
})
