test_that("ci_rc() returns a 2-row tidy data frame with the expected columns", {
  result <- ci_rc(b_j = 0.5, SE_b_j = 0.2, s_Y = 1, s_X = 1, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("lower_limit", "upper_limit"))
})

test_that("ci_rc() return contract matches the help page", {
  # The @return block promises a 2-row table with four columns, terms
  # lower_limit and upper_limit, and tail probabilities beside each limit.
  result <- ci_rc(b_j = 0.61319, SE_b_j = 0.16098, N = 30, p = 6, conf_level = 0.95)
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("lower_limit", "upper_limit"))
  expect_equal(nrow(result), 2L)
  expect_true(is.numeric(result$value))
  expect_s3_class(result, "dmar_tbl")
  expect_equal(attr(result, "conf_level"), 0.95)
  # Central approach: the tail probabilities are the nominal alphas.
  expect_equal(result$prob_less, c(0.025, 0.975))
  expect_equal(result$prob_greater, c(0.975, 0.025))

  # The noncentral path returns the same shape, with achieved tail
  # probabilities near the nominal alphas.
  result_nc <- ci_rc(
    b_j = 0.61319, SE_b_j = 0.16098, N = 30, p = 6,
    conf_level = 0.95, noncentral = TRUE
  )
  expect_named(result_nc, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result_nc$term, c("lower_limit", "upper_limit"))
  expect_equal(result_nc$prob_less, c(0.025, 0.975), tolerance = 1e-4)
  expect_equal(result_nc$prob_greater, c(0.975, 0.025), tolerance = 1e-4)
})

test_that("ci_rc() interval brackets the point estimate", {
  b_j <- 0.5
  result <- ci_rc(b_j = b_j, SE_b_j = 0.2, s_Y = 1, s_X = 1, N = 100, p = 5)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_lt(ll, b_j)
  expect_gt(ul, b_j)
})
