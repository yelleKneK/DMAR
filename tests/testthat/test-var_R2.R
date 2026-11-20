## var_R2() -- variance of the squared multiple correlation.

test_that("var_R2() returns a tidy data.frame with term var_R2", {
  res <- var_R2(true_R2 = 0.3, N = 100, p = 3)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "var_R2")
  expect_gt(res$value, 0)
})

test_that("var_R2() shrinks toward zero as N grows", {
  small <- var_R2(0.3, N = 50,   p = 3)$value
  large <- var_R2(0.3, N = 5000, p = 3)$value
  expect_gt(small, large)
})

test_that("var_R2() approaches zero near the R^2 = 1 boundary", {
  near1 <- var_R2(0.999, N = 100, p = 3)$value
  middle <- var_R2(0.5, N = 100, p = 3)$value
  expect_lt(near1, middle)
})

test_that("var_R2() responds smoothly to p (monotone in N at fixed p)", {
  v_n200 <- var_R2(0.3, N = 200, p = 5)$value
  v_n400 <- var_R2(0.3, N = 400, p = 5)$value
  expect_gt(v_n200, v_n400)
})

test_that("var_R2() takes true_R2; the DMAR 1.0.x name population_R2 stops", {
  expect_identical(names(formals(var_R2)), c("true_R2", "N", "p"))
  expect_error(var_R2(population_R2 = 0.3, N = 100, p = 3), "unused argument")
})

test_that("var_R2() is the variance of the exact distribution of R^2", {
  # Oracle, with no package code: under multivariate normality the sample
  # R^2 is a mixture of Beta(p/2 + j, (N - p - 1)/2) distributions with
  # negative binomial weights Gamma(m + j) / (Gamma(m) j!) rho2^j
  # (1 - rho2)^m, m = (N - 1)/2, so its variance is the mean of the Beta
  # variances plus the variance of the Beta means over the weights (taken
  # in that form, since E[R^4] - E[R^2]^2 loses four digits at rho2 = .9).
  # The lane I5 evidence (work/round6/evidence/i5_names_bp) checks these
  # moments against R^2 from raw multivariate normal data.
  mixture_moments <- function(rho2, N, p) {
    m <- (N - 1) / 2
    b <- (N - p - 1) / 2
    j <- 0:20000
    w <- exp(lgamma(m + j) - lgamma(m) - lgamma(j + 1) + j * log(rho2) +
               m * log1p(-rho2))
    w <- w / sum(w)
    a <- p / 2 + j
    mu <- a / (a + b)
    e1 <- sum(w * mu)
    c(mean = e1, var = sum(w * (a * b / ((a + b)^2 * (a + b + 1)) + (mu - e1)^2)))
  }
  for (cell in list(c(0.3, 100, 3), c(0.05, 20, 2), c(0.6, 40, 6),
                    c(0.9, 500, 10), c(0.99, 50, 3))) {
    o <- mixture_moments(cell[1], cell[2], cell[3])
    expect_equal(var_R2(true_R2 = cell[1], N = cell[2], p = cell[3])$value,
                 o[["var"]], tolerance = 1e-10, info = paste(cell, collapse = " "))
  }
})
