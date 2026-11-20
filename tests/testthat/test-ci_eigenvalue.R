test_that("ci_eigenvalue() returns the largest eigenvalue and its CI", {
  set.seed(113)
  X <- data.frame(matrix(rnorm(200), nrow = 50))
  res <- ci_eigenvalue(X, which_eigenvalue = 1)
  expect_setequal(res$term, c("eigenvalue", "ci_lower", "ci_upper"))
  expect_lt(res$value[res$term == "ci_lower"],
            res$value[res$term == "eigenvalue"])
  expect_gt(res$value[res$term == "ci_upper"],
            res$value[res$term == "eigenvalue"])
})

test_that("ci_eigenvalue() accepts data.frame or matrix", {
  set.seed(113)
  X <- data.frame(matrix(rnorm(200), nrow = 50))
  res_df  <- ci_eigenvalue(X, which_eigenvalue = 1)
  res_mat <- ci_eigenvalue(cov(X), N = 50, which_eigenvalue = 1)
  expect_equal(res_df$value, res_mat$value, tolerance = 1e-12)
})

test_that("ci_eigenvalue() second eigenvalue is smaller than first", {
  set.seed(113)
  X <- data.frame(matrix(rnorm(200), nrow = 50))
  ev1 <- ci_eigenvalue(X, which_eigenvalue = 1)$value[1]
  ev2 <- ci_eigenvalue(X, which_eigenvalue = 2)$value[1]
  expect_gt(ev1, ev2)
})
