test_that("ci_nc_F() returns a 2-row data frame with expected columns", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_equal(result$term, c("ci_lower", "ci_upper"))
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
})

test_that("ci_nc_F() ci_lower < ci_upper for significant F", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  expect_lt(ll, ul)
})

test_that("ci_nc_F() achieves target tail probabilities to high precision", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  alpha <- 0.025
  expect_equal(result$prob_greater[result$term == "ci_lower"], alpha, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "ci_upper"],    alpha, tolerance = 1e-10)
})

test_that("ci_nc_F() one-sided upper interval: ci_lower is 0", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, alpha_lower = 0, alpha_upper = 0.05, conf_level = NULL)
  expect_equal(result$value[result$term == "ci_lower"], 0)
})

test_that("ci_nc_F() one-sided lower interval: ci_upper is Inf", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, alpha_lower = 0.05, alpha_upper = 0, conf_level = NULL)
  expect_equal(result$value[result$term == "ci_upper"], Inf)
})

test_that("ci_nc_F() warns and clamps lower limit when F is not significant", {
  expect_warning(
    result <- ci_nc_F(F_value = 0.5, df_1 = 5, df_2 = 100, conf_level = 0.95),
    "below the alpha_lower critical value"
  )
  expect_equal(result$value[result$term == "ci_lower"], 0)
})

test_that("ci_nc_F() errors on bad inputs", {
  expect_error(ci_nc_F(),                             "F_value")
  expect_error(ci_nc_F(F_value = -1, df_1 = 5, df_2 = 100), "non-negative")
  expect_error(ci_nc_F(F_value = 5),                  "degrees of freedom")
  expect_error(ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95, alpha_lower = 0.025), "Specify either")
  expect_error(ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, conf_level = NULL, alpha_lower = 0.5, alpha_upper = 0.6), "non-negative and sum to less than 1")
})

test_that("ci_nc_F(verbose = FALSE) drops probability columns", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, conf_level = 0.95, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("ci_nc_F() asymmetric alphas use them in the correct direction", {
  result <- ci_nc_F(F_value = 5, df_1 = 5, df_2 = 100, alpha_lower = 0.01, alpha_upper = 0.04, conf_level = NULL)
  expect_equal(result$prob_greater[result$term == "ci_lower"], 0.01, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "ci_upper"],    0.04, tolerance = 1e-10)
})

test_that("ci_nc_F() clamp warning names the consequence for the interval and carries the dmar_nc_F_clamp class", {
  msgs <- capture_warnings(
    result <- ci_nc_F(F_value = 0.5, df_1 = 5, df_2 = 100, conf_level = 0.95)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the noncentrality parameter is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_warning(
    ci_nc_F(F_value = 0.5, df_1 = 5, df_2 = 100, conf_level = 0.95),
    class = "dmar_nc_F_clamp"
  )
  expect_equal(result$value[result$term == "ci_lower"], 0)
})

test_that("ci_nc_F() names the inputs and what to try when the root search fails", {
  expect_error(
    ci_nc_F(F_value = Inf, df_1 = 5, df_2 = 100, conf_level = 0.95),
    "In ci_nc_F\\(\\).*F_value = Inf.*adjust 'tol'"
  )
})

test_that("ci_nc_F() reproduces the published limits of Bargmann's (1970) example", {
  # Venables (1975, p. 411): F = 11.2213 on 4 and 50 df gives the .90
  # interval 19.381 to 71.551, which he reports agrees to three decimals
  # with an accurate solution. Steiger (2004, p. 169) gives 19.380 to
  # 71.549 for F as he reports it, 11.221.
  v <- ci_nc_F(F_value = 11.2213, df_1 = 4, df_2 = 50, conf_level = .90)$value
  expect_equal(round(v, 3), c(19.381, 71.551))
  s <- ci_nc_F(F_value = 11.221, df_1 = 4, df_2 = 50, conf_level = .90)$value
  expect_equal(round(s, 3), c(19.380, 71.549))
})
