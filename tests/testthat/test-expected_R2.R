test_that("expected_R2() returns a 1-row data frame with the expected term", {
  result <- expected_R2(population_R2 = 0.3, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_equal(result$term, "expected_value_population_R2")
})

test_that("expected_R2() expected value is at least the population R2 (positive bias)", {
  pop <- 0.3
  result <- expected_R2(population_R2 = pop, N = 100, p = 5)
  expect_gte(result$value, pop)
})

test_that("expected_R2() bias shrinks toward zero as N grows", {
  pop <- 0.3
  bias_small <- expected_R2(population_R2 = pop, N = 50,    p = 5)$value - pop
  bias_large <- expected_R2(population_R2 = pop, N = 5000,  p = 5)$value - pop
  expect_gt(bias_small, bias_large)
  expect_lt(bias_large, 0.01)
})
