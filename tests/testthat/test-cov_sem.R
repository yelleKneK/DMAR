test_that("cov_sem() returns the lavaan model implied covariance of the observed variables", {
  skip_if_not_installed("lavaan")

  pop_model <- "
    f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
    f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
    f2 ~ 0.5*f1
    f1 ~~ 1*f1
    f2 ~~ 0.75*f2
    y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
    y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
  "

  res   <- cov_sem(pop_model)
  Sigma <- res$sigma_theta
  ov    <- c("y1", "y2", "y3", "y4", "y5", "y6")

  expect_named(res, c("sigma_theta", "mu_theta", "observed_vars"))
  expect_identical(res$observed_vars, ov)

  # No mean structure: the model implied means are zero for every observed
  # variable.
  expect_identical(res$mu_theta, stats::setNames(rep(0, 6L), ov))

  expect_true(is.matrix(Sigma))
  expect_true(is.numeric(Sigma))
  expect_equal(dim(Sigma), c(6L, 6L))
  expect_identical(dimnames(Sigma), list(ov, ov))
  expect_false(anyNA(Sigma))

  # Symmetric.
  expect_true(isSymmetric(unname(Sigma)))

  # Positive definite: all eigenvalues strictly greater than zero.
  expect_true(all(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values > 0))

  # Hand-computed implied (co)variances from the fixed population values.
  # var(f1) = 1; var(f2) = beta^2 var(f1) + psi22 = 0.5^2 * 1 + 0.75 = 1.
  # var(yj)  = lambda_j^2 var(factor) + residual variance.
  # cov within a factor: lambda_j lambda_k var(factor).
  # cov across factors:  lambda_j lambda_k cov(f1, f2), cov(f1, f2) = 0.5.
  expect_equal(diag(Sigma),
               c(y1 = 1.5, y2 = 1.14, y3 = 1.14,
                 y4 = 1.5, y5 = 1.14, y6 = 1.14))
  expect_equal(Sigma["y1", "y2"], 0.8)        # 1 * 0.8 * var(f1)
  expect_equal(Sigma["y2", "y3"], 0.64)       # 0.8 * 0.8 * var(f1)
  expect_equal(Sigma["y1", "y4"], 0.5)        # 1 * 1 * cov(f1, f2)
  expect_equal(Sigma["y2", "y5"], 0.32)       # 0.8 * 0.8 * cov(f1, f2)
})

test_that("cov_sem() errors clearly on bad input", {
  skip_if_not_installed("lavaan")

  expect_error(cov_sem(123),
               "single character string", fixed = TRUE)
  expect_error(cov_sem(c("f1 =~ 1*y1", "f1 ~~ 1*f1")),
               "single character string", fixed = TRUE)
  expect_error(cov_sem(NA_character_),
               "single character string", fixed = TRUE)
})

test_that("cov_sem() returns the model implied means of a mean structure", {
  skip_if_not_installed("lavaan")

  # A linear latent growth curve over four waves. The implied wave means are
  # mean(i) + loading * mean(s): 5.0, 5.3, 5.6, 5.9. The implied variances
  # and covariances follow from var(i) = 1, var(s) = 0.2, cov(i, s) = -0.15,
  # and residual variance 0.5: var(t1) = 1.5, cov(t1, t2) = 0.85.
  pop_lgm <- "
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
    i ~~ 1*i
    s ~~ 0.2*s
    i ~~ -0.15*s
    t1 ~~ 0.5*t1; t2 ~~ 0.5*t2; t3 ~~ 0.5*t3; t4 ~~ 0.5*t4
    t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
    i ~ 5*1
    s ~ 0.3*1
  "
  res <- cov_sem(pop_lgm)
  waves <- paste0("t", 1:4)
  expect_equal(res$mu_theta,
               stats::setNames(c(5.0, 5.3, 5.6, 5.9), waves))
  expect_equal(res$sigma_theta["t1", "t1"], 1.5)
  expect_equal(res$sigma_theta["t1", "t2"], 0.85)

  # A mean structure whose intercepts or latent means are left free is
  # refused, the same way a free variance is.
  free_mean <- "
    i =~ 1*t1 + 1*t2
    i ~~ 1*i
    t1 ~~ 0.5*t1; t2 ~~ 0.5*t2
    t1 ~ 0*1; t2 ~ 0*1
    i ~ 1
  "
  expect_error(cov_sem(free_mean), "parameters are free")
})

test_that("cov_sem() refuses a model whose variances are not all stated", {
  # lavaanify() without auto.var = TRUE adds an unstated variance fixed at
  # zero with free = 0, so the fully-specified guard used to pass it
  # through and return a silently wrong Sigma. Each of these leaves one
  # variance implied rather than stated.
  expect_error(
    cov_sem("x ~~ 1*x\n m ~ 0.4*x\n y ~ 0.35*m + 0.15*x\n y ~~ 0.813*y"),
    "these parameters are free: m ~~ m")
  expect_error(
    cov_sem("f =~ 1*y1 + 0.8*y2 + 0.8*y3\n y1 ~~ 0.5*y1\n y2 ~~ 0.5*y2\n y3 ~~ 0.5*y3"),
    "these parameters are free: f ~~ f")
  # A fully specified model is unaffected, and its Sigma is the closed form.
  S <- cov_sem(paste("x ~~ 1*x", "m ~ 0.4*x", "m ~~ 0.84*m",
                     "y ~ 0.35*m + 0.15*x", "y ~~ 0.813*y",
                     sep = "\n"))$sigma_theta
  expect_equal(unname(S["x", "x"]), 1, tolerance = 1e-10)
  expect_equal(unname(S["m", "m"]), 1, tolerance = 1e-10)
  expect_equal(unname(S["x", "m"]), 0.4, tolerance = 1e-10)
  expect_equal(unname(S["y", "y"]), 1, tolerance = 1e-10)
})
