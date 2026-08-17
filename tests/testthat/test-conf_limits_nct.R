test_that("conf_limits_nct() returns a data frame with expected columns", {
  result <- conf_limits_nct(ncp = 2.83, df = 126, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_true("lower_limit" %in% result$term)
  expect_true("upper_limit" %in% result$term)
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
})

test_that("conf_limits_nct() lower_limit < ncp < upper_limit", {
  ncp <- 2.83
  result <- conf_limits_nct(ncp = ncp, df = 126, conf_level = 0.95)
  ll <- result[result$term == "lower_limit", "value"]
  ul <- result[result$term == "upper_limit", "value"]
  expect_lt(ll, ncp)
  expect_gt(ul, ncp)
})

test_that("conf_limits_nct() one-sided upper interval: lower_limit is -Inf", {
  result <- conf_limits_nct(ncp = 2, df = 30, alpha_lower = 0, alpha_upper = 0.05, conf_level = NULL)
  ll <- result[result$term == "lower_limit", "value"]
  expect_equal(ll, -Inf)
})

test_that("conf_limits_nct() one-sided lower interval: upper_limit is Inf", {
  result <- conf_limits_nct(ncp = 2, df = 30, alpha_lower = 0.05, alpha_upper = 0, conf_level = NULL)
  ul <- result[result$term == "upper_limit", "value"]
  expect_equal(ul, Inf)
})

test_that("conf_limits_nct() errors when df <= 0", {
  expect_error(conf_limits_nct(ncp = 1, df = 0), "degrees of freedom")
})

test_that("conf_limits_nct() t_value alias works identically to ncp", {
  r1 <- conf_limits_nct(ncp = 2.5, df = 50, conf_level = 0.95)
  r2 <- conf_limits_nct(t_value = 2.5, df = 50, conf_level = 0.95)
  expect_equal(r1$value, r2$value)
})

test_that("conf_limits_nct() achieves target tail probabilities to high precision", {
  # P(T > ncp | nu_lower) should equal alpha_lower; P(T < ncp | nu_upper)
  # should equal alpha_upper. These are the achieved tail probabilities of
  # the search.
  result <- conf_limits_nct(ncp = 2.83, df = 126, conf_level = 0.95)
  alpha <- 0.025
  expect_equal(result[result$term == "lower_limit", "prob_greater"], alpha, tolerance = 1e-8)
  expect_equal(result[result$term == "upper_limit", "prob_less"],    alpha, tolerance = 1e-8)
})

test_that("conf_limits_nct() central case is symmetric about zero when ncp = 0", {
  # When ncp = 0 the problem is symmetric: nu_L solves P(T > 0 | nu_L) = alpha
  # and nu_U solves P(T < 0 | nu_U) = alpha, so nu_U = -nu_L by symmetry of
  # the noncentral t about its noncentrality.
  result <- conf_limits_nct(ncp = 0, df = 30, conf_level = 0.95)
  ll <- result[result$term == "lower_limit", "value"]
  ul <- result[result$term == "upper_limit", "value"]
  expect_equal(ul, -ll, tolerance = 1e-10)
})

test_that("conf_limits_nct() asymmetric alphas use them in the correct direction", {
  result <- conf_limits_nct(ncp = 2.83, df = 126, alpha_lower = 0.01, alpha_upper = 0.04, conf_level = NULL)
  expect_equal(result[result$term == "lower_limit", "prob_greater"], 0.01, tolerance = 1e-8)
  expect_equal(result[result$term == "upper_limit", "prob_less"],    0.04, tolerance = 1e-8)
})

test_that("conf_limits_nct(verbose = FALSE) drops the probability columns", {
  result <- conf_limits_nct(ncp = 2.83, df = 126, conf_level = 0.95, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("conf_limits_nct() rejects mixing conf_level with explicit alphas", {
  expect_error(
    conf_limits_nct(ncp = 2.83, df = 126, conf_level = 0.95, alpha_lower = 0.025),
    "Specify either"
  )
})

test_that("conf_limits_nct() rejects alphas that sum to 1 or more", {
  expect_error(
    conf_limits_nct(ncp = 1, df = 30, conf_level = NULL, alpha_lower = 0.6, alpha_upper = 0.5),
    "must be non-negative and sum to less than 1"
  )
})

test_that("conf_limits_nct() handles large ncp without warning at the boundary", {
  # ncp at R's documented stable limit (37.62) should not warn.
  expect_no_warning(conf_limits_nct(ncp = 37.62, df = 200, conf_level = 0.95))
  # Beyond it should warn.
  expect_warning(conf_limits_nct(ncp = 50, df = 200, conf_level = 0.95), "37.62")
})
