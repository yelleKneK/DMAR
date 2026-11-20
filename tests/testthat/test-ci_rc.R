test_that("ci_rc() returns a 3-row tidy data frame with the expected columns", {
  result <- ci_rc(reg_coef_j = 0.5, se_reg_coef_j = 0.2, sd_Y = 1, sd_X = 1, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("ci_lower", "reg_coef", "ci_upper"))
})

test_that("ci_rc() return contract matches the help page", {
  # The @return block promises a 3-row table with four columns, the
  # estimate between its limits, and tail probabilities beside each limit.
  result <- ci_rc(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098, N = 30, p = 6, conf_level = 0.95)
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("ci_lower", "reg_coef", "ci_upper"))
  expect_equal(nrow(result), 3L)
  expect_equal(result$value[2], 0.61319)
  expect_true(is.numeric(result$value))
  expect_s3_class(result, "dmar_tbl")
  expect_equal(attr(result, "conf_level"), 0.95)
  # Central approach: the tail probabilities are the nominal alphas.
  expect_equal(result$prob_less, c(0.025, NA, 0.975))
  expect_equal(result$prob_greater, c(0.975, NA, 0.025))

  # The noncentral path returns the same shape, with achieved tail
  # probabilities near the nominal alphas.
  result_nc <- ci_rc(
    reg_coef_j = 0.61319, se_reg_coef_j = 0.16098, N = 30, p = 6,
    conf_level = 0.95, noncentral = TRUE
  )
  expect_named(result_nc, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result_nc$term, c("ci_lower", "reg_coef", "ci_upper"))
  expect_equal(result_nc$prob_less, c(0.025, NA, 0.975), tolerance = 1e-4)
  expect_equal(result_nc$prob_greater, c(0.975, NA, 0.025), tolerance = 1e-4)
})

test_that("ci_rc() interval brackets the point estimate", {
  b_j <- 0.5
  result <- ci_rc(reg_coef_j = b_j, se_reg_coef_j = 0.2, sd_Y = 1, sd_X = 1, N = 100, p = 5)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  expect_lt(ll, b_j)
  expect_gt(ul, b_j)
})

test_that("ci_rc() takes the coefficient from the squared multiple correlations", {
  # Through the first 1.1.0 builds this route stopped with "arguments imply
  # differing number of rows": the estimate row held the NULL argument.
  N <- 100; p <- 3; R2 <- .30; R2_without <- .20; R2_j <- .10
  sd_Y <- 5; sd_X <- 2
  b <- sqrt((R2 - R2_without) / (1 - R2_j)) * sd_Y / sd_X
  se <- sqrt((1 - R2) / ((1 - R2_j) * (N - p - 1))) * sd_Y / sd_X
  res <- suppressWarnings(ci_rc(sd_Y = sd_Y, sd_X = sd_X, N = N, p = p, R2_Y_X = R2,
                                R2_Y_X_without_j = R2_without, R2_j_X_without_j = R2_j))
  expect_equal(res$term, c("ci_lower", "reg_coef", "ci_upper"))
  expect_equal(res$value, b + c(-1, 0, 1) * qt(.975, N - p - 1) * se, tolerance = 1e-12)
  expect_warning(ci_rc(sd_Y = sd_Y, sd_X = sd_X, N = N, p = p, R2_Y_X = R2,
                       R2_Y_X_without_j = R2_without, R2_j_X_without_j = R2_j),
                 class = "dmar_result_qualified")
})
