# Oracles: Lambda Phi Lambda' + Theta written out with matrix products and
# element by element, and the model implied covariance matrix that lavaan
# reports for a fitted one-factor model.

test_that("covmat_from_cfa() returns a list whose one element is model_implied_cov", {
  result <- covmat_from_cfa(factor_loading = matrix(c(0.7, 0.8, 0.6), ncol = 1),
                            residual_var = c(0.3, 0.2, 0.4))
  expect_type(result, "list")
  expect_identical(names(result), "model_implied_cov")
  expect_true(is.matrix(result$model_implied_cov))
  expect_null(result$population_cov)
})

test_that("covmat_from_cfa() takes factor_loading and residual_var, and refuses the earlier names", {
  expect_identical(names(formals(covmat_from_cfa)),
                   c("factor_loading", "residual_var", "tol_det"))
  # The DMAR 1.0.x names stop with R's "unused argument" error.
  expect_error(
    covmat_from_cfa(lambda = c(.5, .6, .7), psi_squared = c(.75, .64, .51)),
    "unused argument")
  err <- expect_error(covmat_from_cfa(c(.5, .6, .7)),
                      class = "dmar_argument_invalid")
  expect_identical(err$arg, "residual_var")
})

test_that("covmat_from_cfa() equals Lambda Phi Lambda' + Theta with the factor variance fixed at 1", {
  Lambda <- matrix(c(0.7, 0.8, 0.6, 0.45), ncol = 1)
  Phi    <- matrix(1, 1, 1)
  theta  <- c(0.3, 0.2, 0.4, 0.62)
  Theta  <- diag(theta)
  expected <- Lambda %*% Phi %*% t(Lambda) + Theta
  result <- covmat_from_cfa(factor_loading = Lambda, residual_var = theta)
  expect_equal(unname(result$model_implied_cov), unname(expected), tolerance = 1e-14)
  # Element by element: lambda_i lambda_j off the diagonal, lambda_i^2 +
  # theta_i on it.
  S <- result$model_implied_cov
  for (i in 1:4) for (j in 1:4) {
    expect_equal(S[i, j], Lambda[i] * Lambda[j] + (i == j) * theta[i],
                 tolerance = 1e-14)
  }
  # A one-row matrix of loadings is read the same way.
  expect_equal(covmat_from_cfa(factor_loading = t(Lambda),
                               residual_var = theta)$model_implied_cov,
               S, tolerance = 1e-14)
})

test_that("covmat_from_cfa() from fitted estimates reproduces lavaan's model implied covariance matrix", {
  skip_if_not_installed("lavaan")
  items <- c("t5_general_information", "t6_paragraph_comprehension",
             "t7_sentence", "t9_word_meaning")
  model <- paste("f =~", paste(items, collapse = " + "))
  # Factor variance fixed at 1
  fit <- lavaan::cfa(model, data = holzinger_swineford, std.lv = TRUE)
  est <- lavaan::lavInspect(fit, "est")
  implied <- lavaan::lavInspect(fit, "implied")$cov
  result <- covmat_from_cfa(factor_loading = est$lambda,
                            residual_var = diag(est$theta))$model_implied_cov
  expect_equal(unname(result), unname(unclass(implied)), tolerance = 1e-10)

  # A marker loading fixed at 1 with the factor variance estimated: the
  # loadings times the factor standard deviation give the same matrix.
  fit_marker <- lavaan::cfa(model, data = holzinger_swineford)
  est_marker <- lavaan::lavInspect(fit_marker, "est")
  rescaled <- covmat_from_cfa(
    factor_loading = est_marker$lambda * sqrt(est_marker$psi[1, 1]),
    residual_var = diag(est_marker$theta))$model_implied_cov
  expect_equal(unname(rescaled),
               unname(unclass(lavaan::lavInspect(fit_marker, "implied")$cov)),
               tolerance = 1e-10)
  expect_equal(unname(rescaled), unname(result), tolerance = 1e-6)

  # The estimates cfa_1() reports, as the help page example uses them
  fit_dmar <- cfa_1(data = holzinger_swineford, items = items)
  from_dmar <- covmat_from_cfa(
    factor_loading = fit_dmar$estimate[grepl("=~", fit_dmar$syntax)],
    residual_var = fit_dmar$estimate[fit_dmar$syntax %in%
                                       paste(items, "~~", items)])$model_implied_cov
  expect_equal(unname(from_dmar), unname(result), tolerance = 1e-6)
})

test_that("covmat_from_cfa() residuals are uncorrelated (Theta is strictly diagonal)", {
  # The implied residual covariance matrix Sigma - lambda lambda' is
  # diagonal: there is no way to specify a residual covariance between two
  # indicators.
  Lambda <- matrix(c(0.7, 0.8, 0.6), ncol = 1)
  theta  <- c(0.3, 0.2, 0.4)
  S <- covmat_from_cfa(factor_loading = Lambda, residual_var = theta)$model_implied_cov
  residual_cov <- S - tcrossprod(as.numeric(Lambda))
  expect_equal(residual_cov, diag(theta), tolerance = 1e-12)
  expect_equal(residual_cov[upper.tri(residual_cov)], rep(0, 3), tolerance = 1e-12)
})

test_that("covmat_from_cfa() returned matrix is symmetric and positive definite", {
  Lambda <- matrix(c(0.7, 0.8, 0.6), ncol = 1)
  theta  <- c(0.3, 0.2, 0.4)
  S <- covmat_from_cfa(factor_loading = Lambda, residual_var = theta)$model_implied_cov
  expect_equal(S, t(S), tolerance = 1e-12)
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("covmat_from_cfa() recycles a scalar residual_var across indicators", {
  res <- covmat_from_cfa(factor_loading = rep(0.7, 4), residual_var = 0.51)
  expect_equal(dim(res$model_implied_cov), c(4L, 4L))
  expect_equal(diag(res$model_implied_cov), rep(0.7^2 + 0.51, 4), tolerance = 1e-12)
})

test_that("covmat_from_cfa() rejects length-mismatched residual_var", {
  expect_error(covmat_from_cfa(factor_loading = c(0.5, 0.6, 0.7),
                               residual_var = c(0.4, 0.4)),
               "same length as 'factor_loading'",
               class = "dmar_argument_invalid")
})

test_that("covmat_from_cfa() requires at least two indicators and finite values", {
  invalid <- "dmar_argument_invalid"
  expect_error(covmat_from_cfa(factor_loading = 0.7, residual_var = 0.51),
               "at least two", class = invalid)
  expect_error(covmat_from_cfa(factor_loading = c(0.7, NA), residual_var = 0.51),
               "finite", class = invalid)
  expect_error(covmat_from_cfa(factor_loading = c(0.7, 0.6), residual_var = Inf),
               "finite", class = invalid)
  expect_error(covmat_from_cfa(factor_loading = "0.7", residual_var = 0.51),
               "not missing", class = invalid)
  expect_error(covmat_from_cfa(factor_loading = c(0.7, 0.6), residual_var = 0.51,
                               tol_det = "a"),
               "'tol_det' must be a single number", class = invalid)
})

test_that("covmat_from_cfa() refuses a misspelled argument", {
  expect_error(
    covmat_from_cfa(factor_loading = c(0.5, 0.6, 0.7),
                    residual_var = c(0.75, 0.64, 0.51),
                    tol_dt = 1e-3),
    "unused argument")
})

test_that("covmat_from_cfa() honors tol_det", {
  # An obviously well-conditioned matrix should not warn under tol_det = 0.01.
  expect_silent(
    covmat_from_cfa(factor_loading = c(0.5, 0.6, 0.7),
                    residual_var = c(0.75, 0.64, 0.51),
                    tol_det = 0.01)
  )
  # A pathological model with tiny residual variances and unit loadings is
  # near singular, so a generous tol_det should still flag it.
  expect_warning(
    covmat_from_cfa(factor_loading = c(1, 1, 1, 1),
                    residual_var = 1e-10,
                    tol_det = 1),
    "may not be positive-definite",
    class = "dmar_result_qualified"
  )
  # A large negative residual variance makes the determinant negative, so
  # the check warns; the value is still used as given.
  expect_warning(
    S <- covmat_from_cfa(factor_loading = c(.9, .8, 1.1),
                         residual_var = c(.2, .3, -.3))$model_implied_cov,
    "may not be positive-definite")
  expect_equal(diag(S), c(.81 + .2, .64 + .3, 1.21 - .3), tolerance = 1e-14)
})

test_that("covmat_from_cfa() uses a Heywood residual variance as given and need not warn", {
  # The help's example: loadings .9, .8, 1.1 with residual variances .2, .3,
  # and -.01 (or 0). The sign is not checked, the implied matrix is positive
  # definite, its determinant is far above tol_det, and the call is silent.
  lambda <- c(.9, .8, 1.1)
  for (theta in list(c(.2, .3, -.01), c(.2, .3, 0), c(.36, .51, -.05))) {
    S <- expect_silent(covmat_from_cfa(factor_loading = lambda,
                                       residual_var = theta))$model_implied_cov
    oracle <- outer(lambda, lambda)
    diag(oracle) <- lambda^2 + theta
    expect_equal(S, oracle, tolerance = 1e-15)
    expect_gt(det(S), 1e-05)
    expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  # By the determinant formula for lambda lambda' + Theta with Theta
  # diagonal and nonsingular, det = prod(theta) (1 + sum(lambda^2 / theta)).
  theta <- c(.2, .3, -.01)
  S <- covmat_from_cfa(factor_loading = lambda, residual_var = theta)$model_implied_cov
  expect_equal(det(S), prod(theta) * (1 + sum(lambda^2 / theta)), tolerance = 1e-12)
  expect_equal(det(S), .06829, tolerance = 1e-12)
})

test_that("covmat_from_cfa() determinant check is silent with two negative eigenvalues", {
  # Two negative residual variances with small loadings give two negative
  # eigenvalues and a positive determinant (.99), so the determinant check,
  # as documented, does not establish positive definiteness.
  lambda <- c(.1, .1, .1)
  theta <- c(-1, -1, 1)
  S <- expect_silent(covmat_from_cfa(factor_loading = lambda,
                                     residual_var = theta))$model_implied_cov
  expect_equal(det(S), prod(theta) * (1 + sum(lambda^2 / theta)), tolerance = 1e-12)
  expect_equal(det(S), .99, tolerance = 1e-12)
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  expect_identical(sum(eig < 0), 2L)
})
