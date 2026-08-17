test_that("conf_limits_nc_chisq() returns a 2-row data frame with expected columns", {
  result <- conf_limits_nc_chisq(chi_square = 30, df = 15, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_equal(result$term, c("lower_limit", "upper_limit"))
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
})

test_that("conf_limits_nc_chisq() lower_limit < upper_limit for significant Chi", {
  result <- conf_limits_nc_chisq(chi_square = 30, df = 15, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_lt(ll, ul)
})

test_that("conf_limits_nc_chisq() achieves target tail probabilities to high precision", {
  result <- conf_limits_nc_chisq(chi_square = 30, df = 15, conf_level = 0.95)
  alpha <- 0.025
  expect_equal(result$prob_greater[result$term == "lower_limit"], alpha, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "upper_limit"],    alpha, tolerance = 1e-10)
})

test_that("conf_limits_nc_chisq() one-sided upper interval: lower_limit is 0", {
  result <- conf_limits_nc_chisq(chi_square = 30, df = 15, alpha_lower = 0, alpha_upper = 0.05, conf_level = NULL)
  expect_equal(result$value[result$term == "lower_limit"], 0)
})

test_that("conf_limits_nc_chisq() one-sided lower interval: upper_limit is Inf", {
  result <- conf_limits_nc_chisq(chi_square = 30, df = 15, alpha_lower = 0.05, alpha_upper = 0, conf_level = NULL)
  expect_equal(result$value[result$term == "upper_limit"], Inf)
})

test_that("conf_limits_nc_chisq() warns and clamps when chi_square is below central critical value", {
  # chi_square = 5 on 15 df is so small that BOTH tails degenerate: the lower
  # limit clamps to 0 and the upper limit is undefined (NA); assert both warnings.
  expect_warning(
    expect_warning(
      result <- conf_limits_nc_chisq(chi_square = 5, df = 15, conf_level = 0.95),
      "below the alpha_lower critical value"),
    "upper noncentrality limit is undefined"
  )
  expect_equal(result$value[result$term == "lower_limit"], 0)
  expect_true(is.na(result$value[result$term == "upper_limit"]))
})

test_that("conf_limits_nc_chisq() errors on bad inputs", {
  expect_error(conf_limits_nc_chisq(),                              "chi_square")
  expect_error(conf_limits_nc_chisq(chi_square = -1, df = 15),      "non-negative")
  expect_error(conf_limits_nc_chisq(chi_square = 30),               "degrees of freedom")
  expect_error(conf_limits_nc_chisq(chi_square = 30, df = 15, conf_level = 0.95, alpha_lower = 0.025), "Specify either")
  expect_error(conf_limits_nc_chisq(chi_square = 30, df = 15, conf_level = NULL, alpha_lower = 0.5, alpha_upper = 0.6), "non-negative and sum to less than 1")
})

test_that("conf_limits_nc_chisq(verbose = FALSE) drops probability columns", {
  result <- conf_limits_nc_chisq(chi_square = 30, df = 15, conf_level = 0.95, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("conf_limits_nc_chisq() asymmetric alphas use them in the correct direction", {
  # Chi = 50 with df = 15 is highly significant (central upper-tail ~ 1e-5),
  # so a positive lower limit exists at the stricter alpha_lower = 0.01.
  result <- conf_limits_nc_chisq(chi_square = 50, df = 15, alpha_lower = 0.01, alpha_upper = 0.04, conf_level = NULL)
  expect_equal(result$prob_greater[result$term == "lower_limit"], 0.01, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "upper_limit"],    0.04, tolerance = 1e-10)
})
