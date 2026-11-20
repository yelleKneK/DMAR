
test_that("epsilon estimation declines when n - 1 is below the contrast rank", {
  # Same ruling as anova_within_two_way(): with n <= k - 1 the contrast
  # covariance matrix is singular and the sample epsilon is a rank
  # artifact bounded by (n - 1)/(k - 1), not an estimate.
  set.seed(113)
  Y <- matrix(rnorm(3 * 5), nrow = 3)   # n = 3 subjects, k = 5 levels
  expect_warning(res <- epsilon_corrections(Y), "cannot be")
  eps <- setNames(res$epsilon, res$epsilon_method)
  expect_true(is.na(eps[["Greenhouse-Geisser"]]))
  expect_true(is.na(eps[["Huynh-Feldt"]]))
  expect_equal(unname(eps[["lower_bound"]]), 1 / 4)
  # And the well-posed regime still estimates.
  Y2 <- matrix(rnorm(12 * 4), nrow = 12)
  expect_no_warning(res2 <- epsilon_corrections(Y2))
  expect_true(all(!is.na(res2$epsilon)))
})
