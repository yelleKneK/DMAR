test_that("ss_power_rc() returns a tidy data.frame with the documented columns", {
  res <- ss_power_rc(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_N", "actual_power",
                    "noncentral_t_parm", "effect_size") %in% res$term))
})

test_that("ss_power_rc() actual_power is at or above desired_power", {
  res <- ss_power_rc(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5,
                     desired_power = 0.85)
  expect_gte(res$value[res$term == "actual_power"], 0.85 - 1e-6)
})

test_that("ss_power_rc() effect size = sqrt((R2 - R2_without_j) / (1 - R2))", {
  rho2 <- 0.50; rho2_k <- 0.40
  res  <- ss_power_rc(rho2_Y_X = rho2, rho2_Y_X_without_j = rho2_k, p = 5)
  f2 <- (rho2 - rho2_k) / (1 - rho2)
  expect_equal(res$value[res$term == "effect_size"], sqrt(f2),
               tolerance = 1e-9)
})

test_that("ss_power_rc() necessary_N grows as the unique increment shrinks", {
  small <- ss_power_rc(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.49, p = 5)$value[1]
  big   <- ss_power_rc(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5)$value[1]
  expect_gt(small, big)
})
