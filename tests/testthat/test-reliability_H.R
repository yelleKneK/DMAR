test_that("reliability_H() reproduces the Hancock-Mueller (2001) formula", {
  loadings <- c(0.6, 0.7, 0.8)
  theta <- loadings^2 / (1 - loadings^2)
  H_expected <- sum(theta) / (1 + sum(theta))
  res <- reliability_H(loadings)
  expect_equal(res$value[res$term == "reliability_H"], H_expected,
               tolerance = 1e-12)
})

test_that("reliability_H() H monotone increasing in number of indicators", {
  # With identical loadings, adding more indicators must increase H.
  H_two   <- reliability_H(rep(0.7, 2))$value[1]
  H_three <- reliability_H(rep(0.7, 3))$value[1]
  H_six   <- reliability_H(rep(0.7, 6))$value[1]
  expect_gt(H_three, H_two)
  expect_gt(H_six,   H_three)
})

test_that("reliability_H() delta method CI is returned when SEs are supplied", {
  res <- reliability_H(c(0.6, 0.7, 0.8), se_loadings = c(0.05, 0.04, 0.03))
  expect_true(all(c("lower_limit", "upper_limit", "var_H") %in% res$term))
  expect_gte(res$value[res$term == "lower_limit"], 0)
  expect_lte(res$value[res$term == "upper_limit"], 1)
})

test_that("reliability_H() rejects |loading| >= 1", {
  expect_error(reliability_H(c(0.5, 1.0)), "in \\(-1, 1\\)")
})
