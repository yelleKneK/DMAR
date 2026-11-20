test_that("convert_cor_cov() returns a matrix of the right dimensions", {
  R <- diag(3)
  sds <- c(1, 2, 3)
  result <- convert_cor_cov(cor_mat = R, sd = sds)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3L, 3L))
})

test_that("convert_cor_cov() diagonal equals sd^2 for identity correlation", {
  sds <- c(2, 3, 5)
  result <- convert_cor_cov(cor_mat = diag(3), sd = sds)
  expect_equal(diag(result), sds^2)
})

test_that("convert_cor_cov() produces symmetric covariance matrix", {
  R <- matrix(c(1, 0.5, 0.3,
                0.5, 1, 0.4,
                0.3, 0.4, 1), nrow = 3)
  sds <- c(1, 2, 3)
  result <- convert_cor_cov(cor_mat = R, sd = sds)
  expect_true(isSymmetric(result))
})

test_that("convert_cor_cov() off-diagonal equals rho * sd_i * sd_j", {
  R <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
  sds <- c(2, 3)
  result <- convert_cor_cov(cor_mat = R, sd = sds)
  expect_equal(result[1, 2], 0.6 * 2 * 3)
})

test_that("convert_cor_cov() errors on non-square matrix", {
  expect_error(convert_cor_cov(cor_mat = matrix(1:6, nrow = 2), sd = c(1, 2)), "square")
})
