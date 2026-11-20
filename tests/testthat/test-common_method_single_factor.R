test_that("common_method_single_factor() reproduces the factanal one-factor solution", {
  set.seed(1)
  f <- rnorm(200)
  d <- data.frame(
    x1 = f + rnorm(200), x2 = f + rnorm(200), x3 = f + rnorm(200),
    x4 = rnorm(200),     x5 = rnorm(200),     x6 = rnorm(200))
  res <- common_method_single_factor(d)
  expect_s3_class(res, "dmar_tbl")

  # Anchor: the maximum likelihood loadings from factanal() reproduce the
  # reported proportion of common variance to 1e-6.
  R <- cor(d)
  fit <- factanal(covmat = R, factors = 1)
  Lambda <- fit$loadings[, 1]
  expect_equal(res$value[res$term == "variance_explained"],
               sum(Lambda^2) / 6, tolerance = 1e-6)
  expect_equal(res$value[res$term == "variance_explained"],
               0.2248317622, tolerance = 1e-6)
  expect_equal(res$value[res$term == "n_items"], 6)

  # The factor-analytic proportion excludes unique variance, so it sits
  # below the first-eigenvalue version of the screen.
  ev1 <- eigen(R, symmetric = TRUE, only.values = TRUE)$values[1]
  expect_lt(res$value[res$term == "variance_explained"], ev1 / 6)

  # Same answer from the correlation matrix directly.
  res2 <- common_method_single_factor(R = R)
  expect_equal(res2$value, res$value)
})

test_that("common_method_single_factor() separates a true one-factor structure from noise", {
  # Monte Carlo sanity check: data generated from a one-factor model with
  # communalities of .7 must show a high proportion of common variance, and
  # independent noise a low one.
  skip_on_cran()
  set.seed(113)
  strong <- replicate(20, {
    g <- rnorm(500)
    X <- sapply(1:6, function(i) sqrt(0.7) * g + sqrt(0.3) * rnorm(500))
    res <- common_method_single_factor(X)
    res$value[res$term == "variance_explained"]
  })
  noise <- replicate(20, {
    X <- matrix(rnorm(500 * 6), 500, 6)
    res <- common_method_single_factor(X)
    res$value[res$term == "variance_explained"]
  })
  expect_gt(mean(strong), 0.65)
  expect_lt(mean(noise), 0.15)
  expect_gt(min(strong), max(noise))
})

test_that("common_method_single_factor() accepts a covariance matrix", {
  set.seed(113)
  f <- rnorm(200)
  d <- data.frame(
    x1 = f + rnorm(200), x2 = f + rnorm(200), x3 = f + rnorm(200),
    x4 = rnorm(200),     x5 = rnorm(200),     x6 = rnorm(200))
  S <- cov(d)
  res_cov <- common_method_single_factor(S = S)
  expect_s3_class(res_cov, "dmar_tbl")

  # A covariance matrix is standardized to a correlation matrix, so the
  # screen returns the same scale-free quantity as the correlation input.
  res_cor <- common_method_single_factor(R = cov2cor(S))
  expect_equal(res_cov$value, res_cor$value)
  expect_equal(res_cov$value[res_cov$term == "n_items"], 6)

  expect_error(
    common_method_single_factor(S = matrix(c(-1, 0.5, 0.5, 2), 2, 2)),
    "covariance matrix")
})

test_that("common_method_single_factor() validates input", {
  expect_error(common_method_single_factor(), "exactly one")
  expect_error(
    common_method_single_factor(data = matrix(rnorm(18), 6, 3), R = diag(3)),
    "exactly one")
  expect_error(common_method_single_factor(data = matrix(1:6, 6, 1)),
               "three or more")
  expect_error(common_method_single_factor(data = matrix(rnorm(12), 6, 2)),
               "three or more")
  expect_error(common_method_single_factor(R = matrix(2, 2, 2)),
               "correlation matrix")
  expect_error(common_method_single_factor(R = diag(2)),
               "three or more")
})
