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
  expect_true(all(c("eta_squared_partial_ci_lower", "eta_squared_partial_ci_upper",
                    "omega_squared_partial_ci_lower", "omega_squared_partial_ci_upper")
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

test_that("ancova() returns the DMAR table of 1.0.x through the shared helper", {
  res <- ancova(pygmalion, outcome = "iq_8", treatment = "treatment",
                covariates = "iq_pre")

  # DMAR 1.0.x set the class and the conf_level and p_terms attributes by
  # hand; .as_dmar_tbl() sets the same ones, so the object is as it was.
  expect_identical(class(res), c("dmar_tbl", "data.frame"))
  expect_identical(attr(res, "conf_level", exact = TRUE), 0.95)
  expect_identical(attr(res, "p_terms", exact = TRUE),
                   c("p_value", "p_homogeneity_of_regression"))
  expect_setequal(names(attributes(res)),
                  c("names", "row.names", "class", "conf_level", "p_terms"))
  plain <- as.data.frame(res)
  expect_identical(class(plain), "data.frame")
  expect_identical(unclass(res), unclass(plain))

  # The p-value rows print to four decimals, never in scientific notation.
  out <- capture.output(print(res))
  expect_false(any(grepl("[0-9]e[-+][0-9]", out)))
  expect_true(any(grepl("^ *p_value +0\\.0210 *$", out)))
  expect_true(any(grepl("^ *p_homogeneity_of_regression +0\\.0498 *$", out)))

  # tidy() is the term and estimate of each row, with the p-value rows
  # still named; glance() is one row with a column for each term.
  td <- generics::tidy(res)
  expect_identical(td$term, res$term)
  expect_identical(td$estimate, res$value)
  expect_identical(attr(td, "p_terms", exact = TRUE),
                   attr(res, "p_terms", exact = TRUE))
  gl <- generics::glance(res)
  expect_identical(nrow(gl), 1L)
  expect_identical(unname(unlist(gl)), res$value)
})
