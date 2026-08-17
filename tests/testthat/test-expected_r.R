test_that("expected_r() returns a tidy data.frame with the documented columns", {
  res <- expected_r(rho = 0.5, n = 10)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("rho", "n", "expected_r", "bias", "relative_bias"))
})

test_that("expected_r() returns a dmar_tbl", {
  expect_s3_class(expected_r(rho = 0.5, n = 10), "dmar_tbl")
})

test_that("expected_r() handles the boundary cases exactly", {
  # rho = 0 => E[r] = 0 exactly, regardless of n
  expect_equal(expected_r(rho = 0, n = 10)$expected_r, 0)
  expect_equal(expected_r(rho = 0, n = 100)$expected_r, 0)

  # rho = +/- 1 => E[r] = +/- 1 exactly
  expect_equal(expected_r(rho =  1, n = 10)$expected_r,  1)
  expect_equal(expected_r(rho = -1, n = 10)$expected_r, -1)
})

test_that("expected_r() is downwardly biased: |E[r]| < |rho| for 0 < |rho| < 1", {
  cases <- expand.grid(rho = c(0.2, 0.5, 0.8),
                       n   = c(5, 10, 20, 50))
  out <- expected_r(rho = cases$rho, n = cases$n)
  expect_true(all(out$expected_r < out$rho))
  expect_true(all(out$bias > 0))
})

test_that("expected_r() bias shrinks toward zero as n grows", {
  out <- expected_r(rho = 0.5, n = c(5, 10, 20, 50, 100, 500))
  # Bias is monotone decreasing in n at fixed rho.
  expect_true(all(diff(out$bias) < 0))
  # And approaches zero.
  expect_lt(out$bias[length(out$bias)], 1e-3)
})

test_that("expected_r() matches the leading-order asymptotic ~ rho * (1 - (1-rho^2)/(2(n-1))) for large n", {
  # For large n, the exact formula and the first-order asymptotic agree.
  rho <- 0.5
  n   <- 200
  exact <- expected_r(rho = rho, n = n)$expected_r
  asymp <- rho - rho * (1 - rho^2) / (2 * (n - 1))
  expect_equal(exact, asymp, tolerance = 1e-3)
})

test_that("expected_r() matches Monte Carlo simulation (deterministic seed) within MC noise", {
  skip_on_cran()  # 20,000 Monte Carlo replications; the exact boundary and asymptotic anchors above run on CRAN
  set.seed(113)
  M   <- 20000
  rho <- 0.4
  n   <- 30
  rs  <- replicate(M, {
    x <- rnorm(n)
    y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    cor(x, y)
  })
  mc_mean <- mean(rs)
  mc_se   <- sd(rs) / sqrt(M)
  exact   <- expected_r(rho = rho, n = n)$expected_r
  # Within 3 MC standard errors.
  expect_lt(abs(exact - mc_mean), 3 * mc_se)
})

test_that("expected_r() is symmetric: E[r | -rho, n] = -E[r | rho, n]", {
  pos <- expected_r(rho =  0.5, n = 15)$expected_r
  neg <- expected_r(rho = -0.5, n = 15)$expected_r
  expect_equal(pos, -neg, tolerance = 1e-12)
})

test_that("expected_r() rejects malformed inputs", {
  expect_error(expected_r(rho = 1.5, n = 10), "lie in")
  expect_error(expected_r(rho = -2,  n = 10), "lie in")
  expect_error(expected_r(rho = 0.5, n = 3),  "n.+>=.+4")
  expect_error(expected_r(rho = "a", n = 10), "numeric")
  expect_error(expected_r(rho = 0.5, n = "10"), "numeric")
})

test_that("expected_r() vectorizes over rho and n", {
  out <- expected_r(rho = c(0.1, 0.3, 0.5), n = c(10, 20, 30))
  expect_equal(nrow(out), 3L)
  expect_equal(out$rho, c(0.1, 0.3, 0.5))
  expect_equal(out$n,   c(10, 20, 30))
})

test_that("expected_r() relative_bias is NA at rho = 0 and reasonable elsewhere", {
  out <- expected_r(rho = c(0, 0.5), n = c(10, 10))
  expect_true(is.na(out$relative_bias[1]))
  expect_false(is.na(out$relative_bias[2]))
  expect_gt(out$relative_bias[2], 0)
})

test_that("Olkin-Pratt unbiased estimator (worked in @examples) is approximately unbiased over MC", {
  skip_on_cran()  # 20,000 Monte Carlo replications, each summing a hypergeometric series; the exact anchors above run on CRAN
  # If r is observed and we apply the Olkin-Pratt correction, the average
  # should be unbiased for rho. This is the inverse use case of expected_r().
  op_unbiased <- function(r, n, tol = 1e-15, max_iter = 5000) {
    z <- 1 - r^2; c_par <- (n - 2) / 2
    s <- 1; term <- 1
    for (k in seq_len(max_iter)) {
      term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
      s <- s + term
      if (abs(term) < tol * abs(s)) break
    }
    r * s
  }
  set.seed(113)
  M   <- 20000
  rho <- 0.4
  n   <- 15
  unbiased_estimates <- replicate(M, {
    x <- rnorm(n)
    y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    op_unbiased(cor(x, y), n)
  })
  mc_mean <- mean(unbiased_estimates)
  mc_se   <- sd(unbiased_estimates) / sqrt(M)
  # OP estimator should be unbiased (within MC noise) for rho:
  expect_lt(abs(mc_mean - rho), 3 * mc_se)
})

test_that("expected_r() matches the documented closed form computed independently", {
  # The help page prints the Olkin-Pratt (Hotelling, 1953) expectation. These
  # anchors evaluate that expression with an independent hypergeometric and
  # check the package agrees. The points are small n with large rho on
  # purpose: that is where a misstated parameter shows. An earlier version of
  # the help page carried the series at c = (n-1)/2 with the gamma ratio
  # written one step down, and it agreed to three decimals by n = 25 and to
  # five by n = 60, so only small n discriminates. The asymptotic and Monte
  # Carlo checks above cannot see it.
  # Pinned from gsl::hyperg_2F1 (gsl 2.1.8, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  grid <- list(c(0.9, 5), c(0.8, 6), c(0.5, 8), c(0.3, 10), c(0.95, 7))
  documented <- c(0.8687129498217154, 0.7646354230230249, 0.4724694274071386,
                  0.2849994536105499, 0.9393395498646683)
  for (i in seq_along(grid)) {
    rho <- grid[[i]][1]; n <- grid[[i]][2]
    expect_equal(expected_r(rho = rho, n = n)$expected_r,
                 documented[i], tolerance = 1e-10,
                 info = sprintf("rho = %.2f, n = %d", rho, n))
  }
})

test_that("expected_r() bias magnitudes match the values quoted in the help page", {
  # The page quotes +0.021 (rho = 0.5, n = 10) and +0.0065 (n = 30);
  # the exact values are 0.02134122 and 0.00652862. Both are positive:
  # bias is rho - E[r], and the sample r underestimates a positive rho.
  b10 <- expected_r(rho = 0.5, n = 10)$bias
  b30 <- expected_r(rho = 0.5, n = 30)$bias
  expect_gt(b10, 0)
  expect_gt(b30, 0)
  expect_equal(round(b10, 3), 0.021)
  expect_equal(round(b30, 4), 0.0065)
})
