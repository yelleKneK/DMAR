test_that("var_smd() returns exact and approximate rows", {
  res <- var_smd(delta = 0.5, n_1 = 20)
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("var_smd_exact", "var_smd_approx"))
  expect_true(all(res$value > 0))
})

test_that("var_smd() approximation converges to exact as n grows", {
  small <- var_smd(0.5, n_1 = 10)
  large <- var_smd(0.5, n_1 = 1000)
  # Drift between exact and approx is bigger at small n:
  rel_small <- abs(small$value[1] - small$value[2]) / small$value[1]
  rel_large <- abs(large$value[1] - large$value[2]) / large$value[1]
  expect_gt(rel_small, rel_large)
})

test_that("var_smd() unbiased = TRUE rescales by J^2", {
  res_d <- var_smd(0.5, n_1 = 10)
  res_g <- var_smd(0.5, n_1 = 10, unbiased = TRUE)
  df <- 18
  J  <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  expect_equal(res_g$value[1] / res_d$value[1], J^2, tolerance = 1e-8)
})

test_that("var_smd() rejects bad inputs", {
  expect_error(var_smd("a", 10), "numeric")
  expect_error(var_smd(0.5, 2),  ">=")
})
