test_that("ci_smd_c() returns a 3-row tidy data frame with the expected terms", {
  result <- ci_smd_c(ncp = 2.83, n_C = 60, n_E = 70, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("lower_limit", "smd_c", "upper_limit"))
})

test_that("ci_smd_c() lower_limit < smd_c < upper_limit for a positive ncp", {
  result <- ci_smd_c(ncp = 2.83, n_C = 60, n_E = 70, conf_level = 0.95)
  ll  <- result$value[result$term == "lower_limit"]
  est <- result$value[result$term == "smd_c"]
  ul  <- result$value[result$term == "upper_limit"]
  expect_lt(ll, est)
  expect_lt(est, ul)
})

test_that("ci_smd_c() errors when neither ncp nor smd_c is provided", {
  expect_error(
    ci_smd_c(n_C = 60, n_E = 70, conf_level = 0.95),
    "must specify either"
  )
})

test_that("ci_smd_c() back-transforms ncp -> smd_c using sqrt((n_C + n_E)/(n_C * n_E))", {
  # Regression test: a typo in the back-transformation used n_C * n_C instead
  # of n_C * n_E, which gave the wrong smd_c whenever n_C != n_E.
  ncp <- 3.0
  n_C <- 50
  n_E <- 100
  result <- ci_smd_c(ncp = ncp, n_C = n_C, n_E = n_E)
  expected_smd_c <- ncp * sqrt((n_C + n_E) / (n_C * n_E))
  expect_equal(result$value[result$term == "smd_c"], expected_smd_c)
})

test_that("ci_smd_c() ncp <-> smd_c round-trip is consistent for unequal n", {
  # Supplying smd_c directly and supplying the matching ncp must produce the
  # same data frame.
  n_C <- 50; n_E <- 100
  smd_c_val <- 0.5
  ncp_eq <- smd_c_val * sqrt((n_C * n_E) / (n_C + n_E))
  r1 <- ci_smd_c(smd_c = smd_c_val, n_C = n_C, n_E = n_E)
  r2 <- ci_smd_c(ncp    = ncp_eq,   n_C = n_C, n_E = n_E)
  expect_equal(r1$value, r2$value, tolerance = 1e-9)
})

test_that("ci_smd_c() accepts asymmetric alpha bounds", {
  res <- ci_smd_c(smd_c = 0.5, n_C = 50, n_E = 100,
                  conf_level = NULL, alpha_lower = .01, alpha_upper = .04)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, c("lower_limit", "smd_c", "upper_limit"))
})

test_that("ci_smd_c() rejects mixing conf_level and alphas", {
  expect_error(
    ci_smd_c(smd_c = 0.5, n_C = 50, n_E = 100,
             conf_level = .95, alpha_lower = .025, alpha_upper = .025),
    "cannot mix them"
  )
})
