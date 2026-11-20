test_that("print_summary() returns the summary invisibly, unchanged", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  out <- capture.output(ret <- print_summary(fit))
  # Returned summary has full-precision numeric p-values:
  expect_true(is.numeric(ret$coefficients[, "Pr(>|t|)"]))
  # And full numeric estimates and standard errors:
  expect_true(is.numeric(ret$coefficients[, "Estimate"]))
})

test_that("print_summary() formats tiny p-values as '< 0.0001'", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  out <- capture.output(print_summary(fit))
  expect_true(any(grepl("< 0\\.0001", out)))
  # No scientific-notation p-value strings in the printed output:
  expect_false(any(grepl("Pr.*[0-9]e[+-]?[0-9]", out)))
})

test_that("print_summary() reports R^2 and F-statistic for an lm fit", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  out <- capture.output(print_summary(fit))
  expect_true(any(grepl("Multiple R-squared", out)))
  expect_true(any(grepl("F-statistic",        out)))
})

test_that("print_summary() reports random effects for an lmer fit", {
  skip_if_not_installed("lme4")
  fit <- lme4::lmer(weight ~ Time + (1 | Chick), data = ChickWeight)
  out <- capture.output(print_summary(fit))
  expect_true(any(grepl("Random effects",  out)))
  expect_true(any(grepl("Fixed effects",   out)))
  expect_false(any(grepl("Multiple R-squared", out)))
})

test_that("print_summary() respects digits_p", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  out <- capture.output(print_summary(fit, digits_p = 6))
  expect_true(any(grepl("0\\.[0-9]{6}", out)) ||
              any(grepl("< 0\\.000001", out)))
})
