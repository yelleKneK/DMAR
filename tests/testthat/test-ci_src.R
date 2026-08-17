test_that("ci_src() returns a 3-row data frame with the estimate between the limits", {
  result <- ci_src(beta_j = 0.4, SE_beta_j = 0.1, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("lower_limit", "src", "upper_limit"))
  expect_equal(result$value[result$term == "src"], 0.4)
})

test_that("ci_src() interval brackets the standardized point estimate", {
  beta_j <- 0.4
  result <- ci_src(beta_j = beta_j, SE_beta_j = 0.1, N = 100, p = 5)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_lt(ll, beta_j)
  expect_gt(ul, beta_j)
})
