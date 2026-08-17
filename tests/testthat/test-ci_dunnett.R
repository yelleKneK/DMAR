test_that("ci_dunnett() returns documented columns", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- ci_dunnett(fit, control = "ctrl")
  expect_setequal(colnames(res),
                  c("contrast", "mean_difference", "se", "t_statistic",
                    "lower_limit", "upper_limit", "p_adjusted"))
  expect_equal(nrow(res), 2)  # a - 1 comparisons
  # Convention: the wide table is a dmar_tbl with the confidence footer, and
  # the adjusted p column renders to fixed decimals like every p-value.
  expect_s3_class(res, "dmar_tbl")
  expect_identical(attr(res, "conf_level"), 0.95)
  disp <- format(res)
  expect_true(all(grepl("^(< )?[01]\\.[0-9]{4}$", disp$p_adjusted) |
                    grepl("^< 0\\.0001$", disp$p_adjusted)))
})

test_that("ci_dunnett() control flag rejects bad inputs", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  expect_error(ci_dunnett(fit, control = "missing"), "not one of")
})

test_that("ci_dunnett() matches multcomp::glht() on PlantGrowth", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- ci_dunnett(fit, control = "ctrl")

  # Pinned from confint(multcomp::glht()) under set.seed(113), since
  # multcomp's simulated quantile is stochastic (multcomp 1.4.30,
  # 2026-08-09); live comparison in tools/oracle_checks.R. Rows are
  # trt1 - ctrl, trt2 - ctrl.
  ci_estimate <- c(-0.3709999999999998, 0.494)
  ci_lwr      <- c(-1.021547602034485, -0.1565476020344853)
  ci_upr      <- c(0.2795476020344857, 1.144547602034485)

  # The estimated mean differences are the ordinary contrasts of group
  # means and must reproduce multcomp's estimates exactly.
  expect_equal(res$mean_difference, ci_estimate,
               tolerance = 1e-8)
  # The simultaneous Dunnett limits use the same multivariate-t critical
  # value, so the interval endpoints agree with glht()'s confint().
  expect_equal(res$lower_limit, ci_lwr, tolerance = 1e-4)
  expect_equal(res$upper_limit, ci_upr, tolerance = 1e-4)
})

test_that("ci_dunnett() adjusted p-values are deterministic", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  expect_identical(ci_dunnett(fit, control = "ctrl")$p_adjusted,
                   ci_dunnett(fit, control = "ctrl")$p_adjusted)
})

test_that("ci_dunnett() with one treatment reduces to the ordinary t p-value", {
  # With a single treatment-vs-control comparison there is no multiplicity, so
  # the adjusted p-value equals the unadjusted two-sided t p-value.
  d2 <- PlantGrowth[PlantGrowth$group %in% c("ctrl", "trt1"), ]
  d2$group <- factor(d2$group)
  res <- ci_dunnett(lm(weight ~ group, data = d2), control = "ctrl")
  expect_equal(res$p_adjusted, 2 * pt(-abs(res$t_statistic),
               df = nrow(d2) - 2), tolerance = 1e-7)
})

test_that("ci_dunnett() adjusted p-values match multcomp::glht()", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- ci_dunnett(fit, control = "ctrl")
  # Pinned from summary(multcomp::glht())$test$pvalues under set.seed(113),
  # since multcomp's adjusted p-values are stochastic (multcomp 1.4.30,
  # 2026-08-09); live comparison in tools/oracle_checks.R.
  p_glht <- c(0.322695685768782, 0.1534858615151239)
  expect_equal(res$p_adjusted, p_glht, tolerance = 1e-3)
})

test_that("ci_dunnett() one-sided alternatives return half-open CIs", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  r_g <- ci_dunnett(fit, control = "ctrl", alternative = "greater")
  r_l <- ci_dunnett(fit, control = "ctrl", alternative = "less")
  expect_true(all(is.infinite(r_g$upper_limit) & r_g$upper_limit > 0))
  expect_true(all(is.infinite(r_l$lower_limit) & r_l$lower_limit < 0))
})

test_that("ci_dunnett() one-sided intervals contain their point estimates", {
  # A one-sided confidence bound lies on the far side of the point estimate:
  # the "greater" lower limit at or below it, the "less" upper limit at or
  # above it. Before the sign fix the "less" bound was computed with the
  # negative critical value cv_dunnett() reports for that alternative, placing
  # the bound at diff - |d| * se and excluding the point estimate itself.
  fit <- lm(weight ~ group, data = PlantGrowth)
  r_g <- ci_dunnett(fit, control = "ctrl", alternative = "greater")
  r_l <- ci_dunnett(fit, control = "ctrl", alternative = "less")
  expect_true(all(r_g$lower_limit <= r_g$mean_difference))
  expect_true(all(r_l$upper_limit >= r_l$mean_difference))
  # The two alternatives use the same critical value magnitude, so the "less"
  # margin above the estimate mirrors the "greater" margin below it.
  expect_equal(r_l$upper_limit - r_l$mean_difference,
               r_g$mean_difference - r_g$lower_limit, tolerance = 1e-10)
})

test_that("ci_dunnett() one-sided p-values agree with the bound's decision", {
  # Test-interval duality: at family-wise alpha = 1 - conf_level, the adjusted
  # p-value falls below alpha exactly when the one-sided bound excludes zero,
  # because both come from the same equicorrelated multivariate t reference.
  # Before the sign fix both PlantGrowth "less" rows contradicted this (upper
  # limits below zero beside adjusted p-values of 0.16 and 0.99).
  fit <- lm(weight ~ group, data = PlantGrowth)
  r_l <- ci_dunnett(fit, control = "ctrl", alternative = "less")
  r_g <- ci_dunnett(fit, control = "ctrl", alternative = "greater")
  expect_identical(r_l$p_adjusted < 0.05, r_l$upper_limit < 0)
  expect_identical(r_g$p_adjusted < 0.05, r_g$lower_limit > 0)
  # At conf_level = 0.90 the trt2 "greater" comparison rejects (adjusted p =
  # 0.0768 < 0.10), exercising the duality on the rejection side as well.
  r_g90 <- ci_dunnett(fit, control = "ctrl", alternative = "greater",
                      conf_level = 0.90)
  expect_true(any(r_g90$p_adjusted < 0.10))
  expect_identical(r_g90$p_adjusted < 0.10, r_g90$lower_limit > 0)
})

test_that("ci_dunnett() one-sided limits match multcomp::glht()", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  r_l <- ci_dunnett(fit, control = "ctrl", alternative = "less")
  r_g <- ci_dunnett(fit, control = "ctrl", alternative = "greater")
  # Pinned from confint(multcomp::glht()) under set.seed(113), "less"
  # computed before "greater" so the RNG sequence is reproducible (multcomp
  # 1.4.30, 2026-08-09); live comparison in tools/oracle_checks.R. Rows are
  # trt1 - ctrl, trt2 - ctrl.
  ci_l_upr <- c(0.1858958819319817, 1.050895881931981)
  ci_g_lwr <- c(-0.9278958819319814, -0.06289588193198137)
  # Compare on the absolute scale: multcomp's simulated critical value carries
  # Monte Carlo error of a few units in the fourth decimal, which a relative
  # tolerance would magnify for the near-zero trt2 "greater" limit.
  expect_lt(max(abs(r_l$upper_limit - ci_l_upr)), 5e-4)
  expect_lt(max(abs(r_g$lower_limit - ci_g_lwr)), 5e-4)
})

test_that("ci_dunnett() one-sided simultaneous coverage is near nominal", {
  # Monte Carlo confirmation of the audit that found the "less" sign error:
  # with the bound at diff - |d| * se the simultaneous coverage of the true
  # differences measured 0.00115 against the nominal 0.95; sign-fixed coverage
  # measured 0.94685. The fast anchors that stay on CRAN are the
  # point-estimate containment, p-value duality, and multcomp agreement tests
  # above.
  skip_on_cran()
  set.seed(113)
  # G = 300 with bounds of 0.90 and 0.99 keeps the guard against the
  # catastrophic sign error (measured coverage 0.00115) at a quarter of the
  # earlier G = 1200 runtime; each iteration prices two multivariate t
  # quantiles, which is what made this the slowest test in the file. The
  # fine-calibration evidence remains the G = 1200 run recorded above.
  G <- 300
  grp <- factor(rep(paste0("g", 1:4), each = 10))  # 4 groups, n = 10
  covered_l <- covered_g <- logical(G)
  for (i in seq_len(G)) {
    y <- rnorm(40)  # all population means equal, so every true difference is 0
    r_l <- ci_dunnett(y, group = grp, alternative = "less")
    r_g <- ci_dunnett(y, group = grp, alternative = "greater")
    covered_l[i] <- all(r_l$upper_limit >= 0)
    covered_g[i] <- all(r_g$lower_limit <= 0)
  }
  expect_gt(mean(covered_l), 0.90)
  expect_lt(mean(covered_l), 0.99)
  expect_gt(mean(covered_g), 0.90)
  expect_lt(mean(covered_g), 0.99)
})
