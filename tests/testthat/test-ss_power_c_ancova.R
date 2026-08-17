test_that("ss_power_c_ancova() returns a tidy data.frame", {
  res <- ss_power_c_ancova(psi = 1.5, c_weights = c(1, -1),
                           sigma = 2, rho = 0.3)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_per_group", "actual_power",
                    "noncentral_t_parm") %in% res$term))
})

test_that("ss_power_c_ancova() reaches desired_power at the returned n", {
  res <- ss_power_c_ancova(psi = 1.5, c_weights = c(1, -1),
                           sigma = 2, rho = 0.3, desired_power = 0.85)
  expect_gte(res$value[res$term == "actual_power"], 0.85 - 1e-6)
})

test_that("ss_power_c_ancova() necessary_n shrinks as rho grows (covariate adjustment helps)", {
  low_rho  <- ss_power_c_ancova(psi = 1.5, c_weights = c(1, -1),
                                sigma = 2, rho = 0.10)$value[1]
  high_rho <- ss_power_c_ancova(psi = 1.5, c_weights = c(1, -1),
                                sigma = 2, rho = 0.70)$value[1]
  expect_lt(high_rho, low_rho)
})

test_that("ss_power_c_ancova() with rho = 0 returns a value at or near ss_power_c()", {
  ancova_rho_0 <- ss_power_c_ancova(psi = 1.5, c_weights = c(1, -1),
                                    sigma = 2, rho = 0)$value[1]
  anova_only   <- ss_power_c(psi = 1.5, c_weights = c(1, -1), sigma = 2)$value[1]
  # ANCOVA at rho = 0 loses one df relative to ANOVA, so the ANCOVA n
  # should be at most one or two larger -- not smaller.
  expect_gte(ancova_rho_0, anova_only - 1)
  expect_lt(ancova_rho_0, anova_only + 5)
})

test_that("ss_power_c_ancova returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_c_ancova(psi = 0.6, c_weights = c(1, -1), sigma = 1, rho = 0.3, n = n)
    r$value[r$term == "actual_power"] }
  target <- pw(2) - 1e-4
  r <- ss_power_c_ancova(psi = 0.6, c_weights = c(1, -1), sigma = 1, rho = 0.3, desired_power = target)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})
