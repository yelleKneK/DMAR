## Snapshot tests for the display layer.
##
## The print/format conventions (dmar_tbl rounding, whole-number integers,
## fixed-decimal p-values with the "< 0.0001" floor, fixed_terms for
## information criteria, and the format_p / print_anova / print_summary
## helpers) are exactly what testthat snapshots are for: the full rendered
## output is recorded once and any future change to the display rules shows
## up as a reviewable diff instead of slipping past a regex.

test_that("print.dmar_tbl renders a long table with p_terms and a CI footer", {
  # equivalence_smd: three p-values in a shared value column (p_terms route).
  expect_snapshot(print(equivalence_smd(smd = 0.1, n_1 = 50, n_2 = 50,
                                 delta_lower = 0.5, delta_upper = 0.5)))
  # ci_smd: conf_level footer plus mixed integer / fractional values.
  expect_snapshot(print(ci_smd(smd = 0.5, n_1 = 30, n_2 = 30)))
})

test_that("print.dmar_tbl renders a wide table column by column", {
  expect_snapshot(print(ci_mahalanobis(D2 = 103.2, n_1 = 50, n_2 = 50, p = 4)))
})

test_that("format_p renders fixed decimals with the floor label", {
  expect_snapshot(format_p(c(0.234567, 0.05, 0.0499, 1.010542e-05, 0, 1)))
  expect_snapshot(format_p(c(0.234567, 1e-04), digits_p = 2))
})

test_that("print_anova renders an anova table with DMAR p-value formatting", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_snapshot(print_anova(anova(fit)))
})

test_that("print_summary renders an lm summary with DMAR p-value formatting", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_snapshot(print_summary(summary(fit)))
})
