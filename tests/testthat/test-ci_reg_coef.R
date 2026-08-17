test_that("ci_reg_coef() returns a 3-row tidy data frame with the expected columns", {
  result <- ci_reg_coef(b_j = 0.5, SE_b_j = 0.2, N = 100, p = 3)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("lower_limit", "reg_coef", "upper_limit"))
})

test_that("ci_reg_coef() interval brackets the point estimate", {
  b_j <- 0.5
  result <- ci_reg_coef(b_j = b_j, SE_b_j = 0.2, N = 100, p = 3)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_lt(ll, b_j)
  expect_gt(ul, b_j)
})

test_that("ci_reg_coef() symmetric 95% CI matches t-based formula closely", {
  b_j <- 0.5; SE <- 0.2; N <- 100; p <- 3
  result <- ci_reg_coef(b_j = b_j, SE_b_j = SE, N = N, p = p, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  t_crit <- qt(0.975, df = N - p - 1)
  expect_equal(ll, b_j - t_crit * SE, tolerance = 1e-4)
  expect_equal(ul, b_j + t_crit * SE, tolerance = 1e-4)
})
