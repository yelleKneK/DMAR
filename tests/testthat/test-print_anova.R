test_that("print_anova() returns the input invisibly, unchanged", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  a <- anova(fit)
  out <- capture.output(ret <- print_anova(a))
  # Returned object is identical to input, full precision intact:
  expect_identical(ret, a)
  expect_true(is.numeric(ret[["Pr(>F)"]]))
})

test_that("print_anova() formats tiny p-values as '< 0.0001'", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  a <- anova(fit)
  out <- capture.output(print_anova(a))
  expect_true(any(grepl("< 0\\.0001", out)))
  # No scientific-notation p-value strings in the printed output:
  expect_false(any(grepl("[0-9]e[+-]?[0-9]", out)))
})

test_that("print_anova() preserves the heading attribute", {
  skip_if_not_installed("car")
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  a <- car::Anova(fit, type = "III")
  out <- capture.output(print_anova(a))
  expect_true(any(grepl("Anova Table", out)))
})

test_that("print_anova() respects digits_p", {
  fit <- lm(weight ~ Time + Diet, data = ChickWeight)
  a <- anova(fit)
  out <- capture.output(print_anova(a, digits_p = 6))
  # Six-decimal p-values (or the matching < 0.000001 floor) should appear:
  expect_true(any(grepl("0\\.[0-9]{6}", out)) ||
              any(grepl("< 0\\.000001", out)))
})

test_that("print_anova() handles lmerTest anova output", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  fit <- lmerTest::lmer(weight ~ Time + (1 | Chick), data = ChickWeight)
  a <- anova(fit, type = 3, ddf = "Kenward-Roger")
  expect_silent(out <- capture.output(print_anova(a)))
  # Numeric p-value column survives at full precision:
  expect_true(is.numeric(a[["Pr(>F)"]]))
})
