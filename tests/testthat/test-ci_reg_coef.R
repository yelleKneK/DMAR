test_that("ci_reg_coef() returns a 3-row tidy data frame with the expected columns", {
  result <- ci_reg_coef(reg_coef_j = 0.5, se_reg_coef_j = 0.2, N = 100, p = 3)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_equal(result$term, c("ci_lower", "reg_coef", "ci_upper"))
})

test_that("ci_reg_coef() interval brackets the point estimate", {
  b_j <- 0.5
  result <- ci_reg_coef(reg_coef_j = b_j, se_reg_coef_j = 0.2, N = 100, p = 3)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  expect_lt(ll, b_j)
  expect_gt(ul, b_j)
})

test_that("ci_reg_coef() symmetric 95% CI matches t-based formula closely", {
  b_j <- 0.5; SE <- 0.2; N <- 100; p <- 3
  result <- ci_reg_coef(reg_coef_j = b_j, se_reg_coef_j = SE, N = N, p = p, conf_level = 0.95)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  t_crit <- qt(0.975, df = N - p - 1)
  expect_equal(ll, b_j - t_crit * SE, tolerance = 1e-4)
  expect_equal(ul, b_j + t_crit * SE, tolerance = 1e-4)
})

test_that("ci_reg_coef() refusals state the requirement on the tail areas", {
  # A tail area outside [0, .5] is refused, and the message states the
  # range the check enforces.
  expect_error(ci_reg_coef(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098,
                           N = 30, p = 6, alpha_lower = 0.6, alpha_upper = 0.025),
               "'alpha_lower' must be between 0 and .5", fixed = TRUE,
               class = "dmar_argument_invalid")
  expect_error(ci_reg_coef(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098,
                           N = 30, p = 6, alpha_lower = 0.025, alpha_upper = -0.1),
               "'alpha_upper' must be between 0 and .5", fixed = TRUE,
               class = "dmar_argument_invalid")
  # The boundary values 0 and .5 are accepted (a one-sided bound at 0).
  expect_no_error(ci_reg_coef(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098,
                              N = 30, p = 6, alpha_lower = 0, alpha_upper = 0.05))
})
