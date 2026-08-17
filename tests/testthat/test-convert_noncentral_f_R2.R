test_that("convert_f_R2() and convert_R2_f() are inverses", {
  F_in <- 5; df_1 <- 5; df_2 <- 100
  r2 <- convert_f_R2(F_value = F_in, df_1 = df_1, df_2 = df_2)$value
  F_back <- convert_R2_f(R2 = r2, df_1 = df_1, df_2 = df_2)$value
  expect_equal(F_in, F_back, tolerance = 1e-10)
})

test_that("convert_lambda_R2() and convert_R2_lambda() are inverses", {
  lambda_in <- 5; N <- 100
  r2 <- convert_lambda_R2(lambda = lambda_in, N = N)$value
  lambda_back <- convert_R2_lambda(R2 = r2, N = N)$value
  expect_equal(lambda_in, lambda_back, tolerance = 1e-10)
})

test_that("convert_lambda_R2(0, N) is 0 and convert_R2_lambda(0, N) is 0", {
  expect_equal(convert_lambda_R2(lambda = 0, N = 100)$value, 0, tolerance = 1e-12)
  expect_equal(convert_R2_lambda(R2 = 0, N = 100)$value, 0, tolerance = 1e-12)
})

test_that("convert_*() functions return tidy data frames", {
  expect_named(convert_f_R2(F_value = 5, df_1 = 5, df_2 = 100), c("term", "value"))
  expect_named(convert_R2_f(R2 = 0.2, df_1 = 5, df_2 = 100), c("term", "value"))
  expect_named(convert_lambda_R2(lambda = 5, N = 100), c("term", "value"))
  expect_named(convert_R2_lambda(R2 = 0.2, N = 100), c("term", "value"))
})
