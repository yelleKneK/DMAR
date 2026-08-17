test_that("ss_aipe_reg_coef() returns a positive integer sample size", {
  result <- ss_aipe_reg_coef(rho2_Y_X = 0.3, rho2_j_X_without_j = 0.2, p = 3, b_j = 0.4, width = 0.3)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_gt(result$value, 0)
})

test_that("ss_aipe_reg_coef() requires a larger n for a narrower width", {
  n_wide   <- ss_aipe_reg_coef(rho2_Y_X = 0.3, rho2_j_X_without_j = 0.2, p = 3, b_j = 0.4, width = 0.40)$value
  n_narrow <- ss_aipe_reg_coef(rho2_Y_X = 0.3, rho2_j_X_without_j = 0.2, p = 3, b_j = 0.4, width = 0.20)$value
  expect_lt(n_wide, n_narrow)
})
