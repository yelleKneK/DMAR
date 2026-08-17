test_that("contrast_test() returns a tidy data.frame with the documented columns", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("contrast", "estimate", "se", "t", "df",
                      "p_value", "p_adjusted", "ci_lower", "ci_upper"))
})

test_that("pairwise label convention matches TukeyHSD()", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise")
  # TukeyHSD reports trt1-ctrl, trt2-ctrl, trt2-trt1 (later level minus earlier).
  expect_equal(res$contrast,
               c("trt1 - ctrl", "trt2 - ctrl", "trt2 - trt1"))
})

test_that("adjust='tukey' reproduces TukeyHSD() p-values and CI bounds", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise", adjust = "tukey")

  hsd <- TukeyHSD(fit)$group
  expect_equal(res$estimate,   unname(hsd[, "diff"]),  tolerance = 1e-10)
  expect_equal(res$ci_lower, unname(hsd[, "lwr"]),   tolerance = 1e-6)
  expect_equal(res$ci_upper, unname(hsd[, "upr"]),   tolerance = 1e-6)
  expect_equal(res$p_adjusted,      unname(hsd[, "p adj"]), tolerance = 1e-6)
})

test_that("unadjusted pairwise p-values match pairwise.t.test(p.adjust='none')", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise")
  ref <- pairwise.t.test(PlantGrowth$weight, PlantGrowth$group,
                         p.adjust.method = "none")$p.value
  # ref is a triangular matrix: ref[i, j] is p for group_i vs group_j (j < i)
  # Map our 'trt1 - ctrl' etc. to entries in ref.
  expect_equal(res$p_value[res$contrast == "trt1 - ctrl"], ref["trt1", "ctrl"])
  expect_equal(res$p_value[res$contrast == "trt2 - ctrl"], ref["trt2", "ctrl"])
  expect_equal(res$p_value[res$contrast == "trt2 - trt1"], ref["trt2", "trt1"])
})

test_that("Bonferroni p_adj equals min(1, m * p) for m contrasts", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise", adjust = "bonferroni")
  m   <- nrow(res)
  expect_equal(res$p_adjusted, pmin(1, res$p_value * m))
})

test_that("Bonferroni CI is wider than unadjusted CI", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  none <- contrast_test(fit, contrasts = "pairwise", adjust = "none")
  bonf <- contrast_test(fit, contrasts = "pairwise", adjust = "bonferroni")
  expect_true(all((bonf$ci_upper - bonf$ci_lower) >
                  (none$ci_upper - none$ci_lower)))
})

test_that("Scheffé CI is wider than Tukey for the same pairwise contrasts", {
  # Scheffé is conservative; for k = 3 pairwise comparisons it gives wider CIs
  # than Tukey on the same data.
  fit <- aov(weight ~ group, data = PlantGrowth)
  tk <- contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
  sf <- contrast_test(fit, contrasts = "pairwise", adjust = "scheffe")
  expect_true(all((sf$ci_upper - sf$ci_lower) >
                  (tk$ci_upper - tk$ci_lower)))
})

test_that("a t-statistic squared at one-df equals the F for a single contrast", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(
    fit,
    contrasts = list("ctrl vs trts" = c(1, -0.5, -0.5))
  )
  # SS_psi / MS_error has F ~ F(1, df_residual); equals t^2.
  expect_equal(res$t^2,
               res$estimate^2 / res$se^2,
               tolerance = 1e-12)
})

test_that("Welch (var_equal = FALSE) gives a different df per contrast", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise", var_equal = FALSE)
  # When variances differ across groups, Welch dfs vary across pairs.
  expect_true(length(unique(res$df)) > 1L)
  # And they should not all equal the residual df 27.
  expect_false(all(res$df == 27))
})

test_that("Welch SE matches a hand calculation for one pairwise contrast", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise", var_equal = FALSE)
  # ctrl vs trt1 by hand
  v <- tapply(PlantGrowth$weight, PlantGrowth$group, var)
  n <- as.integer(table(PlantGrowth$group))
  expected_se <- sqrt(v[["ctrl"]] / n[1] + v[["trt1"]] / n[2])
  expect_equal(res$se[res$contrast == "trt1 - ctrl"],
               unname(expected_se), tolerance = 1e-10)
})

test_that("user-supplied list of contrasts is preserved with names", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(
    fit,
    contrasts = list("ctrl vs trts" = c(1, -0.5, -0.5),
                     "trt1 vs trt2" = c(0, 1, -1))
  )
  expect_equal(res$contrast, c("ctrl vs trts", "trt1 vs trt2"))
})

test_that("matrix input is accepted and unnamed rows get C1, C2, ...", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  cm  <- rbind(c(1, -0.5, -0.5), c(0, 1, -1))
  res <- contrast_test(fit, contrasts = cm)
  expect_equal(res$contrast, c("C1", "C2"))
})

test_that("a single numeric vector is treated as one contrast", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = c(1, -0.5, -0.5))
  expect_equal(nrow(res), 1L)
  expect_equal(res$contrast, "C1")
})

test_that("non-pairwise contrasts with adjust='tukey' raise an informative error", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  expect_error(
    contrast_test(fit,
                  contrasts = list("ctrl vs trts" = c(1, -0.5, -0.5)),
                  adjust = "tukey"),
    "pairwise"
  )
})

test_that("contrast weights not summing to zero produce a warning but still compute", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  expect_warning(
    res <- contrast_test(fit, contrasts = c(1, 1, 0)),
    "do not sum to zero"
  )
  expect_equal(nrow(res), 1L)
})

test_that("multi-way models are rejected with an informative error", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  expect_error(contrast_test(fit), "one-way")
})

test_that("attributes record the user's choices", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise",
                      adjust = "bonferroni", conf_level = 0.99,
                      var_equal = FALSE)
  expect_equal(attr(res, "adjust"),     "bonferroni")
  expect_equal(attr(res, "conf_level"), 0.99)
  expect_equal(attr(res, "var_equal"),  FALSE)
})

test_that("tidy()/glance() return broom columns matching the source table", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res <- contrast_test(fit, contrasts = "pairwise",
                      adjust = "bonferroni", conf_level = 0.95)

  td <- generics::tidy(res)
  expect_equal(nrow(td), nrow(res))
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper",
                     "statistic", "df", "p_value", "p_adjusted",
                     "conf_level"))
  expect_equal(td$term,      res$contrast)
  expect_equal(td$estimate,  res$estimate)
  expect_equal(td$ci_lower,  res$ci_lower)
  expect_equal(td$ci_upper, res$ci_upper)
  expect_equal(td$statistic, res$t)
  expect_equal(td$df,        res$df)
  expect_equal(td$p_value,   res$p_value)
  expect_equal(td$p_adjusted, res$p_adjusted)
  expect_true(all(td$conf_level == 0.95))

  gl <- generics::glance(res)
  expect_equal(nrow(gl), 1L)
  expect_named(gl, c("n_contrasts", "adjust", "var_equal",
                     "p_adjusted_min", "conf_level"))
  expect_equal(gl$n_contrasts, nrow(res))
  expect_equal(gl$adjust, "bonferroni")
  expect_equal(gl$p_adjusted_min, min(res$p_adjusted))
  expect_equal(gl$conf_level, 0.95)

  # The result is a dmar_tbl and carries the broom-dispatch subclass.
  expect_s3_class(res, "dmar_contrast_test")
  expect_s3_class(res, "dmar_tbl")
})
