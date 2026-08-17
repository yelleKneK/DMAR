test_that("ci_cv() returns the expected 4-row tidy data frame", {
  result <- ci_cv(cv = 0.2, n = 30, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_setequal(result$term, c("lower_limit", "upper_limit", "c_of_v", "c_of_v_unbiased"))
})

test_that("ci_cv() lower_limit < c_of_v < upper_limit for a typical positive cv", {
  result <- ci_cv(cv = 0.2, n = 30, conf_level = 0.95)
  ll  <- result$value[result$term == "lower_limit"]
  est <- result$value[result$term == "c_of_v"]
  ul  <- result$value[result$term == "upper_limit"]
  expect_lt(ll, est)
  expect_lt(est, ul)
})

test_that("ci_cv() accepts (mean, sd) in place of cv and gives consistent limits", {
  m <- 10; s <- 2; n <- 30
  r1 <- ci_cv(cv = s / m, n = n, conf_level = 0.95)
  r2 <- ci_cv(mean = m, sd = s, n = n, conf_level = 0.95)
  expect_equal(
    r1$value[r1$term %in% c("lower_limit", "upper_limit")],
    r2$value[r2$term %in% c("lower_limit", "upper_limit")],
    tolerance = 1e-8
  )
})
