test_that("returns a 1-row tidy table with the method-specific term", {
  op <- unbiased_R2(R2 = 0.50, N = 50, p = 5)
  expect_s3_class(op, "data.frame")
  expect_equal(nrow(op), 1L)
  expect_equal(op$term, "unbiased_population_R2")

  ez <- unbiased_R2(R2 = 0.50, N = 50, p = 5, method = "ezekiel")
  expect_equal(ez$term, "adjusted_population_R2")
  expect_type(op$value, "double")
  expect_type(ez$value, "double")
})

test_that("the Ezekiel option reproduces summary(lm)$adj.r.squared exactly", {
  set.seed(113)
  d   <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
  fit <- lm(V1 ~ ., data = d)
  s   <- summary(fit)
  ez  <- unbiased_R2(R2 = s$r.squared, N = 50, p = 5, method = "ezekiel")
  expect_equal(ez$value, s$adj.r.squared)
})

test_that("the closed form matches a direct Ezekiel computation", {
  R2 <- 0.42; N <- 80; p <- 4
  expect_equal(
    unbiased_R2(R2, N, p, method = "ezekiel")$value,
    1 - ((N - 1) / (N - p - 1)) * (1 - R2)
  )
})

test_that("both corrections shrink a typical R2 downward, and Olkin-Pratt over Ezekiel", {
  R2 <- 0.50; N <- 50; p <- 5
  op <- unbiased_R2(R2, N, p)$value
  ez <- unbiased_R2(R2, N, p, method = "ezekiel")$value
  expect_lt(op, R2)
  expect_lt(ez, R2)
  expect_gt(op, ez)         # Ezekiel over-shrinks relative to Olkin-Pratt
})

test_that("estimates are allowed to fall below zero for very small R2 (not truncated)", {
  ez <- unbiased_R2(R2 = 0.001, N = 30, p = 8, method = "ezekiel")$value
  op <- unbiased_R2(R2 = 0.001, N = 30, p = 8)$value
  expect_lt(ez, 0)
  expect_lt(op, 0)
})

test_that("R2 = 0 is handled (nudged off the 2F1 boundary) and stays finite", {
  op <- unbiased_R2(R2 = 0, N = 50, p = 5)$value
  expect_true(is.finite(op))
})

test_that("the correction shrinks toward the observed R2 as N grows", {
  R2 <- 0.50; p <- 5
  small <- unbiased_R2(R2, 50,  p)$value
  large <- unbiased_R2(R2, 500, p)$value
  expect_lt(abs(large - R2), abs(small - R2))
})

test_that("invalid arguments error informatively", {
  expect_error(unbiased_R2(R2 = 1.2, N = 50, p = 5), "\\[0, 1\\]")
  expect_error(unbiased_R2(R2 = -0.1, N = 50, p = 5), "\\[0, 1\\]")
  expect_error(unbiased_R2(R2 = 0.5, N = 50, p = 0), "positive integer")
  expect_error(unbiased_R2(R2 = 0.5, N = 50, p = 2.5), "positive integer")
  expect_error(unbiased_R2(R2 = 0.5, N = 6, p = 5), "greater than 'p' \\+ 1")
})

test_that("the Olkin-Pratt estimator is essentially unbiased by Monte Carlo", {
  skip_on_cran()
  set.seed(113)
  N <- 40; p <- 3; rho2 <- 0.25; G <- 3000
  b <- sqrt(rho2 / p)                 # independent predictors => population R^2 = rho2
  raw <- op <- ez <- numeric(G)
  for (g in seq_len(G)) {
    X <- matrix(rnorm(N * p), N, p)
    y <- as.numeric(X %*% rep(b, p)) + sqrt(1 - rho2) * rnorm(N)
    R2 <- summary(lm(y ~ X))$r.squared
    raw[g] <- R2
    op[g]  <- unbiased_R2(R2, N, p)$value
    ez[g]  <- unbiased_R2(R2, N, p, method = "ezekiel")$value
  }
  # Olkin-Pratt is essentially unbiased; raw R^2 is biased upward.
  expect_lt(abs(mean(op) - rho2), 0.01)
  expect_gt(mean(raw) - rho2, 0.02)
  # Olkin-Pratt has smaller absolute bias than the raw estimator.
  expect_lt(abs(mean(op) - rho2), abs(mean(raw) - rho2))
})
