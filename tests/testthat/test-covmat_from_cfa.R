test_that("covmat_from_cfa() returns a list containing population_cov", {
  result <- covmat_from_cfa(lambda = matrix(c(0.7, 0.8, 0.6), ncol = 1), psi_squared = c(0.3, 0.2, 0.4))
  expect_type(result, "list")
  expect_true("population_cov" %in% names(result))
  expect_true(is.matrix(result$population_cov))
})

test_that("covmat_from_cfa() reconstructs Sigma = Lambda Lambda' + diag(Psi)", {
  Lambda <- matrix(c(0.7, 0.8, 0.6), ncol = 1)
  Psi    <- c(0.3, 0.2, 0.4)
  expected <- Lambda %*% t(Lambda) + diag(Psi)
  result <- covmat_from_cfa(lambda = Lambda, psi_squared = Psi)
  expect_equal(unname(result$population_cov), unname(expected), tolerance = 1e-12)
})

test_that("covmat_from_cfa() errors are uncorrelated (Psi is strictly diagonal)", {
  # The implied error covariance matrix Sigma - lambda lambda' is diagonal:
  # there is no way to specify a residual covariance between two indicators.
  Lambda <- matrix(c(0.7, 0.8, 0.6), ncol = 1)
  Psi    <- c(0.3, 0.2, 0.4)
  S <- covmat_from_cfa(lambda = Lambda, psi_squared = Psi)$population_cov
  residual_cov <- S - tcrossprod(as.numeric(Lambda))
  expect_equal(residual_cov, diag(Psi), tolerance = 1e-12)
  expect_equal(residual_cov[upper.tri(residual_cov)], rep(0, 3), tolerance = 1e-12)
})

test_that("covmat_from_cfa() returned matrix is symmetric and positive definite", {
  Lambda <- matrix(c(0.7, 0.8, 0.6), ncol = 1)
  Psi    <- c(0.3, 0.2, 0.4)
  S <- covmat_from_cfa(lambda = Lambda, psi_squared = Psi)$population_cov
  expect_equal(S, t(S), tolerance = 1e-12)
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("covmat_from_cfa() recycles a scalar psi_squared across indicators", {
  res <- covmat_from_cfa(lambda = rep(0.7, 4), psi_squared = 0.51)
  expect_equal(dim(res$population_cov), c(4L, 4L))
  expect_equal(diag(res$population_cov), rep(0.7^2 + 0.51, 4), tolerance = 1e-12)
})

test_that("covmat_from_cfa() rejects length-mismatched psi_squared", {
  expect_error(covmat_from_cfa(lambda = c(0.5, 0.6, 0.7),
                               psi_squared = c(0.4, 0.4)),
               "same length as 'lambda'")
})

test_that("covmat_from_cfa() requires at least two indicators", {
  expect_error(covmat_from_cfa(lambda = 0.7, psi_squared = 0.51),
               "at least two")
})

test_that("covmat_from_cfa() warns when a typo is passed through `...`", {
  expect_warning(
    covmat_from_cfa(lambda = c(0.5, 0.6, 0.7),
                    psi_squared = c(0.75, 0.64, 0.51),
                    tol_dt = 1e-3),
    "Unrecognized argument"
  )
})

test_that("covmat_from_cfa() honors tol_det passed through `...`", {
  # An obviously well-conditioned matrix should not warn under tol_det = 0.01.
  expect_silent(
    covmat_from_cfa(lambda = c(0.5, 0.6, 0.7),
                    psi_squared = c(0.75, 0.64, 0.51),
                    tol_det = 0.01)
  )
  # A pathological model with tiny psi_squared and unit loadings is near
  # singular, so a generous tol_det should still flag it.
  expect_warning(
    covmat_from_cfa(lambda = c(1, 1, 1, 1),
                    psi_squared = 1e-10,
                    tol_det = 1),
    "may not be positive-definite"
  )
})

