test_that("simulate_anova_data() returns a long-format data frame with the expected columns", {
  result <- simulate_anova_data(mu = c(10, 12, 14), sigma = 2, a = 3, n = 5, seed = 1)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("group", "y"))
})

test_that("simulate_anova_data() returns a*n total rows with each group of size n", {
  a <- 3; n <- 5
  result <- simulate_anova_data(mu = c(10, 12, 14), sigma = 2, a = a, n = n, seed = 1)
  expect_equal(nrow(result), a * n)
  expect_true(all(table(result$group) == n))
})

test_that("simulate_anova_data() seed is reproducible", {
  r1 <- simulate_anova_data(mu = c(10, 12, 14), sigma = 2, a = 3, n = 5, seed = 42)
  r2 <- simulate_anova_data(mu = c(10, 12, 14), sigma = 2, a = 3, n = 5, seed = 42)
  expect_identical(r1, r2)
})
