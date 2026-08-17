test_that("common_method_marker() applies the Lindell-Whitney adjustment", {
  R <- matrix(c(1, .5, .4, .5, 1, .45, .4, .45, 1), 3, 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  res <- common_method_marker(R, marker_r = 0.10)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(res$value[res$term == "marker_correlation"], 0.10)

  adj <- attr(res, "adjusted")
  expect_equal(adj["a", "b"], (0.5 - 0.1) / (1 - 0.1))
  expect_equal(diag(adj), c(a = 1, b = 1, c = 1))
  expect_equal(res$value[res$term == "mean_abs_r_adjusted"],
               mean(abs(adj[upper.tri(adj)])))
})

test_that("common_method_marker() uses the smallest positive r as a proxy", {
  R <- matrix(c(1, .5, .2, .5, 1, .45, .2, .45, 1), 3, 3)
  res <- common_method_marker(R)
  expect_equal(res$value[res$term == "marker_correlation"], 0.2)
})

test_that("common_method_marker() validates input", {
  expect_error(common_method_marker(matrix(2, 2, 2)), "correlation matrix")
  expect_error(common_method_marker(diag(3), marker_r = 1.2), "\\(-1, 1\\)")
})
