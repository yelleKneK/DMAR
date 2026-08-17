test_that("conf_limits_ncf() returns a 2-row data frame with expected columns", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_equal(result$term, c("lower_limit", "upper_limit"))
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
})

test_that("conf_limits_ncf() lower_limit < upper_limit for significant F", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_lt(ll, ul)
})

test_that("conf_limits_ncf() achieves target tail probabilities to high precision", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  alpha <- 0.025
  expect_equal(result$prob_greater[result$term == "lower_limit"], alpha, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "upper_limit"],    alpha, tolerance = 1e-10)
})

test_that("conf_limits_ncf() one-sided upper interval: lower_limit is 0", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, alpha_lower = 0, alpha_upper = 0.05, conf_level = NULL)
  expect_equal(result$value[result$term == "lower_limit"], 0)
})

test_that("conf_limits_ncf() one-sided lower interval: upper_limit is Inf", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, alpha_lower = 0.05, alpha_upper = 0, conf_level = NULL)
  expect_equal(result$value[result$term == "upper_limit"], Inf)
})

test_that("conf_limits_ncf() warns and clamps lower limit when F is not significant", {
  expect_warning(
    result <- conf_limits_ncf(F_value = 0.5, df_1 = 5, df_2 = 100, conf_level = 0.95),
    "below the alpha_lower critical value"
  )
  expect_equal(result$value[result$term == "lower_limit"], 0)
})

test_that("conf_limits_ncf() errors on bad inputs", {
  expect_error(conf_limits_ncf(),                             "F_value")
  expect_error(conf_limits_ncf(F_value = -1, df_1 = 5, df_2 = 100), "non-negative")
  expect_error(conf_limits_ncf(F_value = 5),                  "degrees of freedom")
  expect_error(conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95, alpha_lower = 0.025), "Specify either")
  expect_error(conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, conf_level = NULL, alpha_lower = 0.5, alpha_upper = 0.6), "non-negative and sum to less than 1")
})

test_that("conf_limits_ncf(verbose = FALSE) drops probability columns", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("conf_limits_ncf() asymmetric alphas use them in the correct direction", {
  result <- conf_limits_ncf(F_value = 5, df_1 = 5, df_2 = 100, alpha_lower = 0.01, alpha_upper = 0.04, conf_level = NULL)
  expect_equal(result$prob_greater[result$term == "lower_limit"], 0.01, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "upper_limit"],    0.04, tolerance = 1e-10)
})

test_that("conf_limits_ncf() clamp warning names the consequence for the interval and carries the dmar_ncf_clamp class", {
  msgs <- capture_warnings(
    result <- conf_limits_ncf(F_value = 0.5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the noncentrality parameter is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_warning(
    conf_limits_ncf(F_value = 0.5, df_1 = 5, df_2 = 100, conf_level = 0.95),
    class = "dmar_ncf_clamp"
  )
  expect_equal(result$value[result$term == "lower_limit"], 0)
})

test_that("conf_limits_ncf() names the inputs and what to try when the root search fails", {
  expect_error(
    conf_limits_ncf(F_value = Inf, df_1 = 5, df_2 = 100, conf_level = 0.95),
    "In conf_limits_ncf\\(\\).*F_value = Inf.*adjust 'tol'"
  )
})
