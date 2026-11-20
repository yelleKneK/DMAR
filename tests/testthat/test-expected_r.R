test_that("expected_r() returns a tidy data.frame with the documented columns", {
  res <- expected_r(true_cor = 0.5, N = 10)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("true_cor", "N", "expected_r", "bias", "relative_bias"))
})

test_that("expected_r() returns a dmar_tbl", {
  expect_s3_class(expected_r(true_cor = 0.5, N = 10), "dmar_tbl")
})

test_that("expected_r() handles the boundary cases exactly", {
  # rho = 0 => E[r] = 0 exactly, regardless of n
  expect_equal(expected_r(true_cor = 0, N = 10)$expected_r, 0)
  expect_equal(expected_r(true_cor = 0, N = 100)$expected_r, 0)

  # rho = +/- 1 => E[r] = +/- 1 exactly
  expect_equal(expected_r(true_cor =  1, N = 10)$expected_r,  1)
  expect_equal(expected_r(true_cor = -1, N = 10)$expected_r, -1)
})

test_that("expected_r() is downwardly biased: |E[r]| < |rho| for 0 < |rho| < 1", {
  cases <- expand.grid(rho = c(0.2, 0.5, 0.8),
                       n   = c(5, 10, 20, 50))
  out <- expected_r(true_cor = cases$rho, N = cases$n)
  expect_true(all(out$expected_r < out$true_cor))
  expect_true(all(out$bias > 0))
})

test_that("expected_r() bias shrinks toward zero as n grows", {
  out <- expected_r(true_cor = 0.5, N = c(5, 10, 20, 50, 100, 500))
  # Bias is monotone decreasing in n at fixed rho.
  expect_true(all(diff(out$bias) < 0))
  # And approaches zero.
  expect_lt(out$bias[length(out$bias)], 1e-3)
})

test_that("expected_r() matches the leading-order asymptotic ~ rho * (1 - (1-rho^2)/(2(n-1))) for large n", {
  # For large n, the exact formula and the first-order asymptotic agree.
  rho <- 0.5
  n   <- 200
  exact <- expected_r(true_cor = rho, N = n)$expected_r
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
  exact   <- expected_r(true_cor = rho, N = n)$expected_r
  # Within 3 MC standard errors.
  expect_lt(abs(exact - mc_mean), 3 * mc_se)
})

test_that("expected_r() is symmetric: E[r | -rho, n] = -E[r | rho, n]", {
  pos <- expected_r(true_cor =  0.5, N = 15)$expected_r
  neg <- expected_r(true_cor = -0.5, N = 15)$expected_r
  expect_equal(pos, -neg, tolerance = 1e-12)
})

test_that("expected_r() rejects malformed inputs", {
  expect_error(expected_r(true_cor = 1.5, N = 10), "lie in")
  expect_error(expected_r(true_cor = -2,  N = 10), "lie in")
  expect_error(expected_r(true_cor = 0.5, N = 3),  "'N' must be >= 4")
  expect_error(expected_r(true_cor = "a", N = 10), "numeric")
  expect_error(expected_r(true_cor = 0.5, N = "10"), "numeric")
})

test_that("expected_r() vectorizes over rho and N", {
  out <- expected_r(true_cor = c(0.1, 0.3, 0.5), N = c(10, 20, 30))
  expect_equal(nrow(out), 3L)
  expect_equal(out$true_cor, c(0.1, 0.3, 0.5))
  expect_equal(out$N,   c(10, 20, 30))
})

test_that("expected_r() relative_bias is NA at rho = 0 and reasonable elsewhere", {
  out <- expected_r(true_cor = c(0, 0.5), N = c(10, 10))
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
    expect_equal(expected_r(true_cor = rho, N = n)$expected_r,
                 documented[i], tolerance = 1e-10,
                 info = sprintf("rho = %.2f, n = %d", rho, n))
  }
})

test_that("expected_r() bias magnitudes match the values quoted in the help page", {
  # The page quotes +0.021 (rho = 0.5, n = 10) and +0.0065 (n = 30);
  # the exact values are 0.02134122 and 0.00652862. Both are positive:
  # bias is rho - E[r], and the sample r underestimates a positive rho.
  b10 <- expected_r(true_cor = 0.5, N = 10)$bias
  b30 <- expected_r(true_cor = 0.5, N = 30)$bias
  expect_gt(b10, 0)
  expect_gt(b30, 0)
  expect_equal(round(b10, 3), 0.021)
  expect_equal(round(b30, 4), 0.0065)
})

# The gamma ratio Gamma(n/2)^2 / {Gamma((n - 1)/2) Gamma((n + 1)/2)} is
# c4(n)^2, the square of E(S) / sigma on n - 1 df (R/c4_internals.R). Through
# 1.0.x it was the sum of three lgamma() values, which nearly cancel at a
# large n: off by 1e-7 at n = 1e8, larger than the relative bias itself.
# The bias was rho - E[r], which keeps only the digits the two do not share;
# it is now rho w, with the relative bias w = 1 - 2F1 c4(n)^2 formed from
# logarithms. The pins were computed with Rmpfr at 320 bits, summing the
# 2F1 series to 1e-40 and taking c4(n)^2 = exp(log(2 / (n - 1)) +
# 2 lgamma(n / 2) - 2 lgamma((n - 1) / 2)), as w, bias = rho w, and
# E[r] = rho (1 - w).
er_pin <- data.frame(
  rho = rep(c(0.5, 0.9), 6),
  n = rep(c(10, 100, 1000, 1e6, 1e8, 1e9), each = 2),
  w = c(0.042682442326507693, 0.012609995873731108,
        0.0037995057459232097, 0.00097501127177053959,
        0.00037549247723497236, 9.5244950772162376e-5,
        3.7500049218779004e-7, 9.5000244388062233e-8,
        3.7500000492187503e-9, 9.5000002443875036e-10,
        3.7500000049218750e-10, 9.5000000244387481e-11),
  bias = c(0.021341221163253847, 0.011348996286357998,
           0.0018997528729616049, 0.00087751014459348566,
           0.00018774623861748618, 8.5720455694946141e-5,
           1.8750024609389502e-7, 8.5500219949256011e-8,
           1.8750000246093751e-9, 8.5500002199487535e-10,
           1.8750000024609375e-10, 8.5500000219948735e-11),
  expected_r = c(0.47865877883674615, 0.88865100371364202,
                 0.49810024712703840, 0.89912248985540654,
                 0.49981225376138251, 0.89991427954430508,
                 0.49999981249975391, 0.89999991449978007,
                 0.49999999812499998, 0.89999999914500000,
                 0.49999999981250000, 0.89999999991450002)
)

test_that("expected_r() keeps a positive bias of the right size at a very large n", {
  # At rho = .5 the 1.0.x bias was 26 percent too large at n = 1e7, 31
  # times its value at n = 1e8, and negative at n = 3e8, where E[r] came out
  # above rho.
  rho <- 0.5
  for (n in c(1e7, 1e8, 3e8)) {
    res <- expected_r(true_cor = rho, N = n)
    expect_gt(res$bias, 0)
    expect_lt(res$expected_r, rho)
    expect_equal(res$bias, rho * (1 - rho^2) / (2 * (n - 1)),
                 tolerance = 2e-6, info = paste("n =", n))
  }
})

test_that("expected_r() bias keeps its full relative precision at every n", {
  skip_on_cran()
  # Formed as rho - E[r], with the gamma ratio right to a few units in the
  # last place, the bias at rho = .5 is off by 3e-9 at n = 1e6, 3e-7 at
  # n = 1e8, and 7e-6 at n = 1e9, relative. Formed directly it is right to
  # about 1e-15 there. At n = 10 and rho = .9 the stopping rule of the 2F1
  # series leaves about 1e-14.
  res <- expected_r(true_cor = er_pin$rho, N = er_pin$n)
  for (i in seq_len(nrow(er_pin))) {
    info <- sprintf("rho = %.1f, n = %g", er_pin$rho[i], er_pin$n[i])
    tol <- if (er_pin$n[i] < 100) 5e-14 else 1e-14
    expect_equal(res$relative_bias[i], er_pin$w[i], tolerance = tol, info = info)
    expect_equal(res$bias[i], er_pin$bias[i], tolerance = tol, info = info)
    expect_equal(res$expected_r[i], er_pin$expected_r[i], tolerance = 1e-14,
                 info = info)
  }
  # The bias does not depend on the sign of rho except through rho itself.
  neg <- expected_r(true_cor = -er_pin$rho, N = er_pin$n)
  expect_identical(neg$bias, -res$bias)
  expect_identical(neg$relative_bias, res$relative_bias)
})

test_that("expected_r() is unchanged at ordinary n by the stable gamma ratio", {
  # The 1.0.x expression, with the 2F1 series summed here independently.
  # At n of 200 or less its lgamma() difference loses fewer than four
  # digits, and the two agree to 1e-12.
  old_form <- function(rho, n) {
    z <- rho^2; c_par <- (n + 1) / 2
    s <- 1; term <- 1; k <- 0
    repeat {
      k <- k + 1
      term <- term * (k - 0.5)^2 / ((c_par + k - 1) * k) * z
      s <- s + term
      if (abs(term) < 1e-17 * abs(s)) break
    }
    rho * s * exp(2 * lgamma(n / 2) - lgamma((n - 1) / 2) -
                    lgamma((n + 1) / 2))
  }
  grid <- expand.grid(rho = c(-0.95, -0.5, 0.05, 0.3, 0.6, 0.9, 0.99),
                      n = c(4:60, seq(70, 200, by = 10)))
  res <- expected_r(true_cor = grid$rho, N = grid$n)
  ref <- mapply(old_form, grid$rho, grid$n)
  expect_lt(max(abs(res$expected_r / ref - 1)), 1e-12)
  # The bias rho - E[r] of the 1.0.x expression was right to 1e-12 only
  # where the bias is large, at small n and moderate rho (its error is
  # under 6e-13 for n of 30 or less and |rho| of .6 or less); there the
  # direct bias agrees with it to 1e-12.
  small <- grid$n <= 30 & abs(grid$rho) <= 0.6
  expect_lt(max(abs(res$bias[small] / (grid$rho - ref)[small] - 1)), 1e-12)
})

test_that("expected_r() warns about the series where 1.0.x did", {
  # The series now runs on until a term is small against the terms after
  # the leading 1, which is what the bias needs, but the warning still
  # turns on the whole sum. At rho = .999 and n = 5 the series met the
  # 1.0.x rule within 5000 terms but not the new one, and gives no warning;
  # at n = 4 neither rule is met within 5000 terms.
  expect_no_warning(expected_r(true_cor = 0.999, N = 5))
  expect_no_warning(expected_r(true_cor = 0.9995, N = 6))
  expect_warning(expected_r(true_cor = 0.999, N = 4), "did not converge")
})

test_that("expected_r() is right when rho^2 underflows to 0", {
  # Below about 1.6e-162, rho^2 is 0 in double precision. Through 1.0.x
  # that gave E[r] = 0 and a relative bias of 1; the relative bias is
  # 1 - c4(n)^2 there, its limit as rho goes to 0, which is Var(S) / sigma^2
  # on n - 1 df.
  res <- expected_r(true_cor = c(1e-170, 1e-8), N = 10)
  expect_gt(res$expected_r[1], 0)
  expect_equal(res$relative_bias[1], res$relative_bias[2], tolerance = 1e-15)
  expect_equal(res$relative_bias[1], DMAR:::.chi_var_factor(9),
               tolerance = 1e-15)
})

test_that("the size of the bias relative to sampling variability stated on ?expected_r", {
  skip_on_cran()
  # The help page says the bias of a single r is about rho / (2 sqrt(N - 1))
  # of its standard deviation (exact ratio 0.09 at rho = .8, N = 20) and,
  # for the average of K = 50 such correlations, about
  # rho sqrt(K) / (2 sqrt(N - 1)) of its standard error (exact ratio 0.63).
  # Oracle, with no package code: the exact density of r under bivariate
  # normality in Fisher's (1915) integral form,
  #   f(r) = (N - 2) / pi (1 - rho^2)^((N - 1) / 2) (1 - r^2)^((N - 4) / 2)
  #          int_0^Inf (cosh w - rho r)^(-(N - 1)) dw,
  # integrated numerically for the mean and standard deviation of r.
  dens_r <- function(r, rho, N) {
    inner <- sapply(r, function(ri)
      integrate(function(w) (cosh(w) - rho * ri)^(-(N - 1)), 0, Inf,
                rel.tol = 1e-12)$value)
    (N - 2) / pi * (1 - rho^2)^((N - 1) / 2) * (1 - r^2)^((N - 4) / 2) * inner
  }
  moment <- function(k, rho, N)
    integrate(function(r) r^k * dens_r(r, rho, N), -1, 1, rel.tol = 1e-10)$value
  rho <- 0.8; N <- 20; K <- 50
  expect_equal(moment(0, rho, N), 1, tolerance = 1e-8)
  mean_r <- moment(1, rho, N)
  sd_r <- sqrt(moment(2, rho, N) - mean_r^2)
  bias <- expected_r(true_cor = rho, N = N)$bias
  expect_equal(bias, rho - mean_r, tolerance = 1e-8)
  expect_equal(round(bias / sd_r, 2), 0.09)
  expect_equal(round(bias / (sd_r / sqrt(K)), 2), 0.63)
  # The approximations quoted on the page are within 5 percent of these.
  expect_equal(bias / sd_r, rho / (2 * sqrt(N - 1)), tolerance = 0.05)
  expect_equal(bias / (sd_r / sqrt(K)), rho * sqrt(K) / (2 * sqrt(N - 1)),
               tolerance = 0.05)
})
