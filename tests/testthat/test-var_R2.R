## var_R2() -- asymptotic variance of the squared multiple correlation.

test_that("var_R2() returns a tidy data.frame with term var_R2", {
  res <- var_R2(population_R2 = 0.3, N = 100, p = 3)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "var_R2")
  expect_gt(res$value, 0)
})

test_that("var_R2() shrinks toward zero as N grows", {
  small <- var_R2(0.3, N = 50,   p = 3)$value
  large <- var_R2(0.3, N = 5000, p = 3)$value
  expect_gt(small, large)
})

test_that("var_R2() approaches zero near the R^2 = 1 boundary", {
  near1 <- var_R2(0.999, N = 100, p = 3)$value
  middle <- var_R2(0.5, N = 100, p = 3)$value
  expect_lt(near1, middle)
})

test_that("var_R2() responds smoothly to p (monotone in N at fixed p)", {
  v_n200 <- var_R2(0.3, N = 200, p = 5)$value
  v_n400 <- var_R2(0.3, N = 400, p = 5)$value
  expect_gt(v_n200, v_n400)
})
