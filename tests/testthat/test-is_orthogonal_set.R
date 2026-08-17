test_that("is_orthogonal_set() detects orthogonal contrasts under equal n", {
  cmat <- cbind(
    c_linear = c(-3, -1,  1,  3),
    c_quad   = c( 1, -1, -1,  1)
  )
  res <- is_orthogonal_set(cmat)
  expect_equal(res$value[res$term == "all_orthogonal"], 1)
  expect_equal(res$value[res$term == "all_contrasts_sum_to_zero"], 1)
})

test_that("is_orthogonal_set() catches non-orthogonal pairs", {
  cmat <- cbind(
    c_1 = c(1, -1,  0,  0),
    c_2 = c(1,  0, -1,  0)
  )
  res <- is_orthogonal_set(cmat)
  expect_equal(res$value[res$term == "all_orthogonal"], 0)
})

test_that("is_orthogonal_set() flags non-zero column sums", {
  cmat <- cbind(c_bad = c(1, 1, 0, 0), c_ok = c(1, -1, 1, -1))
  res <- is_orthogonal_set(cmat)
  expect_equal(res$value[res$term == "all_contrasts_sum_to_zero"], 0)
})

test_that("is_orthogonal_set() with unequal n uses c_i c_j / n_i", {
  cmat <- cbind(c_1 = c(-1, 0, 1), c_2 = c(1, -2, 1))
  expect_equal(is_orthogonal_set(cmat)$value[1], 1)  # equal-n orthogonal
  res_un <- is_orthogonal_set(cmat, n = c(10, 5, 10))
  expect_true(is.finite(res_un$value[res_un$term == "all_orthogonal"]))
})
