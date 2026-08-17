test_that("ci_sm() returns the expected 3-row tidy data frame", {
  result <- ci_sm(sm = 0.5, N = 50, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("lower_limit", "std_mean", "upper_limit"))
})

test_that("ci_sm() lower_limit < std_mean < upper_limit", {
  result <- ci_sm(sm = 0.5, N = 50, conf_level = 0.95)
  ll  <- result$value[result$term == "lower_limit"]
  est <- result$value[result$term == "std_mean"]
  ul  <- result$value[result$term == "upper_limit"]
  expect_lt(ll, est)
  expect_lt(est, ul)
})

test_that("ci_sm() accepts (mean, sd) in place of sm and gives consistent results", {
  m <- 5; s <- 2; n <- 50
  r1 <- ci_sm(sm = m / s, N = n, conf_level = 0.95)
  r2 <- ci_sm(mean = m, sd = s, N = n, conf_level = 0.95)
  expect_equal(r1$value, r2$value, tolerance = 1e-8)
})

test_that("ci_sm() labels the confidence level correctly when alpha bounds are supplied", {
  # A supplied pair of tail areas defines the coverage; the default conf_level
  # must not be attached (it would mislabel a 90% interval as 95%).
  a <- ci_sm(sm = .5, N = 50, alpha_lower = .05, alpha_upper = .05)
  expect_true(is.null(attr(a, "conf_level")))
  expect_equal(a$value, ci_sm(sm = .5, N = 50, conf_level = .90)$value,
               tolerance = 1e-9)
  # The conf_level path still labels correctly.
  expect_equal(attr(ci_sm(sm = .5, N = 50, conf_level = .95), "conf_level"), .95)
  # Mixing an explicit conf_level with alpha bounds is rejected.
  expect_error(
    ci_sm(sm = .5, N = 50, conf_level = .99, alpha_lower = .05, alpha_upper = .05),
    "cannot mix")
})
