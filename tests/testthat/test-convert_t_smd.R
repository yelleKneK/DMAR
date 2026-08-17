test_that("convert_delta_lambda() and convert_lambda_delta() are inverses", {
  delta <- 0.5; n_1 <- 60; n_2 <- 70
  lambda <- convert_delta_lambda(delta = delta, n_1 = n_1, n_2 = n_2)$value
  delta_back <- convert_lambda_delta(lambda = lambda, n_1 = n_1, n_2 = n_2)$value
  expect_equal(delta, delta_back, tolerance = 1e-12)
})

test_that("convert_delta_lambda() returns the expected tidy data frame", {
  result <- convert_delta_lambda(delta = 0.5, n_1 = 60, n_2 = 70)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "delta_lambda")
})

test_that("convert_lambda_delta() matches the textbook formula", {
  lambda <- 2.83; n_1 <- 60; n_2 <- 70
  result <- convert_lambda_delta(lambda = lambda, n_1 = n_1, n_2 = n_2)
  expect_equal(result$value, lambda * sqrt((n_1 + n_2) / (n_1 * n_2)), tolerance = 1e-12)
})
