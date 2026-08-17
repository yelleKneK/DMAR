test_that("ss_power_factorial_anova() returns a tidy data.frame with the documented terms", {
  res <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                  effect_indices = 1, f = 0.25)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_per_cell", "total_N", "effect_df",
                    "error_df", "noncentrality", "actual_power") %in% res$term))
})

test_that("ss_power_factorial_anova() power at returned n_per_cell meets desired_power", {
  res <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                  effect_indices = 1, f = 0.25,
                                  desired_power = 0.85)
  expect_gte(res$value[res$term == "actual_power"], 0.85 - 1e-6)
})

test_that("ss_power_factorial_anova() total_N = n_per_cell * product(factor_levels)", {
  res <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                  effect_indices = 1, f = 0.25)
  n <- res$value[res$term == "necessary_n_per_cell"]
  expect_equal(res$value[res$term == "total_N"], n * 2 * 3)
})

test_that("ss_power_factorial_anova() effect_df = product(level-1) over the requested effect", {
  # Main effect on a 2-level factor: df = 1
  res1 <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                   effect_indices = 1, f = 0.25)
  expect_equal(res1$value[res1$term == "effect_df"], 1)
  # Two-way interaction on a 2 x 3: df = (2-1)(3-1) = 2
  res2 <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                   effect_indices = c(1, 2), f = 0.25)
  expect_equal(res2$value[res2$term == "effect_df"], 2)
})

test_that("ss_power_factorial_anova() smaller f -> larger n_per_cell", {
  big_eff   <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                        effect_indices = 1, f = 0.40)$value[1]
  small_eff <- ss_power_factorial_anova(factor_levels = c(2, 3),
                                        effect_indices = 1, f = 0.15)$value[1]
  expect_gt(small_eff, big_eff)
})
