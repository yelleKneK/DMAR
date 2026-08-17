test_that("ss_power_rm_anova() returns a tidy data.frame with the documented terms", {
  res <- ss_power_rm_anova(a = 4, f = 0.25)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_subjects", "a", "effect_df", "error_df",
                    "noncentrality", "actual_power") %in% res$term))
})

test_that("ss_power_rm_anova() actual_power reaches desired_power", {
  res <- ss_power_rm_anova(a = 4, f = 0.25, desired_power = 0.85)
  expect_gte(res$value[res$term == "actual_power"], 0.85 - 1e-6)
})

test_that("ss_power_rm_anova() necessary_n grows when rho shrinks (correlation lifts power)", {
  # With higher within-subject correlation, fewer subjects are needed for
  # the within-subjects test.
  low_rho  <- ss_power_rm_anova(a = 4, f = 0.25, rho = 0.10)$value[1]
  high_rho <- ss_power_rm_anova(a = 4, f = 0.25, rho = 0.70)$value[1]
  expect_lt(high_rho, low_rho)
})

test_that("ss_power_rm_anova() necessary_n grows under sphericity correction (epsilon < 1)", {
  full     <- ss_power_rm_anova(a = 4, f = 0.25, epsilon = 1.00)$value[1]
  adjusted <- ss_power_rm_anova(a = 4, f = 0.25, epsilon = 0.60)$value[1]
  expect_gt(adjusted, full)
})

test_that("ss_power_rm_anova() effect_df = a - 1 and error_df = (n-1)(a-1)", {
  a <- 5
  res <- ss_power_rm_anova(a = a, f = 0.25)
  n   <- res$value[res$term == "necessary_n_subjects"]
  expect_equal(res$value[res$term == "effect_df"], a - 1)
  expect_equal(res$value[res$term == "error_df"], (n - 1) * (a - 1))
})

test_that("ss_power_rm_anova() f and eta_squared (partial) give compatible answers", {
  # eta^2 / (1 - eta^2) = f^2 (Cohen), so f = 0.25 -> eta^2 ~ 0.0588
  via_f   <- ss_power_rm_anova(a = 4, f = 0.25)$value[1]
  via_eta <- ss_power_rm_anova(a = 4, eta_squared = 0.25^2 / (1 + 0.25^2))$value[1]
  expect_equal(via_f, via_eta)  # the two parameterizations are exactly equivalent
})

test_that("ss_power_rm_anova returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_rm_anova(a = 3, f = 0.6, rho = 0.5, n = n)
    r$value[r$term == "actual_power"] }
  target <- pw(2) - 1e-4
  r <- ss_power_rm_anova(a = 3, f = 0.6, rho = 0.5, desired_power = target)
  expect_equal(r$value[r$term == "necessary_n_subjects"], 2)
})
