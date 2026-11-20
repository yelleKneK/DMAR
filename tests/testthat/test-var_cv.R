# var_cv(): the first-order variance of the sample coefficient of variation
# k = s / M under normality, kappa^2 (1/2 + kappa^2) / (N - 1), in one row.
#
# The reference variance used below does not use that formula. Under
# normality s and M are independent, so the variance of k is
# E[s^2] E[1/M^2] - (E[s] E[1/M])^2. The moments of 1/M are not finite over
# the whole line, so M is kept within `z_max` of its standard errors of mu
# and at least `floor_mean` times mu (mu = 1 and sigma = kappa), and they are
# computed with base R integrate().
restricted_var_k <- function(kappa, N, z_max = 8, floor_mean = 0) {
  se <- kappa / sqrt(N)
  lo <- max(-z_max, (floor_mean - 1) / se)
  stopifnot(1 + lo * se > 0)
  mass <- stats::pnorm(z_max) - stats::pnorm(lo)
  m1 <- stats::integrate(function(z) stats::dnorm(z) / (1 + se * z),
                         lo, z_max, rel.tol = 1e-12,
                         subdivisions = 2000L)$value / mass
  m2 <- stats::integrate(function(z) stats::dnorm(z) / (1 + se * z)^2,
                         lo, z_max, rel.tol = 1e-12,
                         subdivisions = 2000L)$value / mass
  c4 <- sqrt(2 / (N - 1)) * exp(lgamma(N / 2) - lgamma((N - 1) / 2))
  kappa^2 * m2 - (kappa * c4 * m1)^2
}

test_that("var_cv() returns one row, var_cv_normal, with its closed form", {
  res <- var_cv(true_cv = 0.20, N = 30)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(names(res), c("term", "value"))
  expect_identical(res$term, "var_cv_normal")
  expect_type(res$value, "double")
  # 0.04 * 0.54 / 29
  expect_equal(res$value, 0.000744827586206897, tolerance = 1e-12)
  for (kappa in c(0.05, 0.3, 0.8, 1.5)) {
    for (N in c(3, 10, 30.5, 400)) {
      expect_equal(var_cv(kappa, N)$value,
                   (kappa^2 / 2 + kappa^4) / (N - 1), tolerance = 1e-14,
                   info = paste(kappa, N))
    }
  }
  # The arguments by position are the coefficient and then the sample size.
  expect_identical(var_cv(0.20, 30)$value, res$value)
  expect_identical(names(formals(var_cv)), c("true_cv", "N"))
})

test_that("var_cv() is within 3 percent of an independently computed variance", {
  # kappa = .3, N = 10, the mean within 8 standard errors of mu. var_cv()
  # is 1.4 percent above the reference. With N in place of N - 1 the ratio
  # would be 0.913, without the kappa^4 term 0.859, and with kappa^2
  # multiplied by (N + 1)/(N - 1) in the second term 1.049.
  ref <- restricted_var_k(0.3, 10)
  expect_equal(ref, 0.005817623423, tolerance = 1e-8)
  ratio <- var_cv(true_cv = 0.3, N = 10)$value / ref
  expect_gt(ratio, 1)
  expect_lt(ratio, 1.03)
  expect_equal(ratio, 1.01416, tolerance = 1e-5)
})

test_that("the accuracy figures on ?var_cv hold", {
  # Percent by which var_cv() differs from the variance of k with the mean
  # within 8 standard errors of mu and at least mu / 50.
  pct <- function(kappa, N, z_max = 8, floor_mean = 1 / 50)
    100 * (var_cv(kappa, N)$value /
             restricted_var_k(kappa, N, z_max = z_max,
                              floor_mean = floor_mean) - 1)
  # Too large for kappa <= .3, most for a small kappa. As kappa goes to 0
  # the excess goes to that of sigma^2 / (2 (N - 1)) over the exact variance
  # of s, sigma^2 (1 - c4^2): 3.0067, 1.3676, 0.8843, and 0.5180 percent at
  # N = 10, 20, 30, and 50, which the page gives to one decimal place.
  c4 <- function(N) sqrt(2 / (N - 1)) * exp(lgamma(N / 2) - lgamma((N - 1) / 2))
  sizes <- c(10, 20, 30, 50)
  limit <- 100 * ((1 / (2 * (sizes - 1))) / (1 - c4(sizes)^2) - 1)
  expect_equal(limit, c(3.0067, 1.3676, 0.8843, 0.5180), tolerance = 1e-4)
  expect_equal(round(limit, 1), c(3.0, 1.4, 0.9, 0.5))
  small <- vapply(sizes, function(N) pct(0.01, N), numeric(1))
  expect_equal(small, limit, tolerance = 1e-3)
  expect_true(all(small < limit))
  for (N in c(10, 20, 30, 50, 100, 500)) {
    over <- vapply(c(0.01, 0.05, 0.1, 0.2, 0.3), function(k) pct(k, N),
                   numeric(1))
    expect_true(all(over > 0), info = N)
    expect_true(all(diff(over) < 0), info = N)
  }
  # Too small at kappa = .5: 5.7 percent at N = 10 and 1.7 percent at 30.
  expect_equal(round(c(pct(0.5, 10), pct(0.5, 30)), 1), c(-5.7, -1.7))
  # Too small from kappa = .4 on, and more so the larger kappa is.
  for (N in c(5, 10, 30, 100)) {
    under <- vapply(c(0.4, 0.45, 0.5), function(k) pct(k, N), numeric(1))
    expect_true(all(under < 0), info = N)
    expect_true(all(diff(under) < 0), info = N)
  }
  # The figures of the page are the same for a lower bound on the mean from
  # mu / 1000 to mu / 5, with the mean within 6, 8, or 10 standard errors.
  # With a lower bound of mu / 10000 the figure at kappa = .5 and N = 10 is
  # 5.8 percent, so the page states the range.
  page <- list(c(0.01, 10, 3.0), c(0.01, 20, 1.4), c(0.01, 30, 0.9),
               c(0.01, 50, 0.5), c(0.5, 10, -5.7), c(0.5, 30, -1.7))
  for (cell in page) {
    for (floor_mean in c(1 / 1000, 1 / 100, 1 / 10, 1 / 5)) {
      for (z_max in c(6, 8, 10)) {
        expect_equal(round(pct(cell[1], cell[2], z_max, floor_mean), 1),
                     cell[3],
                     info = paste(cell[1], cell[2], z_max, floor_mean))
      }
    }
  }
  expect_equal(round(pct(0.5, 10, floor_mean = 1e-4), 1), -5.8)
})

test_that("the comparison with a population that is not normal on ?var_cv holds", {
  # The large-sample variance of k is xi^2 / N, with xi^2 written here in
  # the central moments of Kelley, Darku, and Chattopadhyay (2018, Equation
  # 30); N var_cv() approaches kappa^2 (1/2 + kappa^2).
  xi2 <- function(mu, sigma, mu3, mu4)
    mu4 / (4 * mu^2 * sigma^2) - sigma^2 / (4 * mu^2) - mu3 / mu^3 +
      sigma^4 / mu^4
  normal_theory <- function(kappa) 1e6 * var_cv(kappa, 1e6 + 1)$value
  # A normal population: the two agree.
  expect_equal(xi2(1, 0.3, 0, 3 * 0.3^4), normal_theory(0.3),
               tolerance = 1e-12)
  # A gamma population with shape a has mean a, variance a, third central
  # moment 2a, and fourth 3a^2 + 6a; kappa = .8 is a = 1 / .64.
  a <- 1 / 0.8^2
  gamma_xi2 <- xi2(a, sqrt(a), 2 * a, 3 * a^2 + 6 * a)
  expect_equal(gamma_xi2, 0.5248, tolerance = 1e-12)
  expect_equal(round(100 * (normal_theory(0.8) / gamma_xi2 - 1)), 39)
  # A Laplace population with mean 1 and standard deviation .3 has third
  # central moment 0 and fourth 6 sd^4.
  laplace_xi2 <- xi2(1, 0.3, 0, 6 * 0.3^4)
  expect_equal(laplace_xi2, 0.1206, tolerance = 1e-12)
  expect_equal(round(100 * (1 - normal_theory(0.3) / laplace_xi2)), 56)
})

test_that("var_cv() is near a simulated variance", {
  skip_on_cran()
  # 40,000 samples of 40 from a normal population with kappa = .15; the
  # first-order value is within a few percent of the simulated variance.
  set.seed(20261110)
  k <- replicate(40000, {
    y <- stats::rnorm(40, 10, 1.5)
    stats::sd(y) / mean(y)
  })
  expect_lt(abs(stats::var(k) / var_cv(.15, 40)$value - 1), 0.05)
})

test_that("var_cv() refuses a coefficient that is not one finite positive number", {
  bad <- list(NA_real_, NaN, NA, Inf, -Inf, 0, -0.1, c(0.2, 0.3), numeric(0),
              "0.2", NULL, TRUE, list(0.2))
  for (x in bad) {
    err <- expect_error(var_cv(true_cv = x, N = 30),
                        class = "dmar_argument_invalid",
                        info = paste(deparse(x), collapse = ""))
    expect_identical(err$arg, "true_cv")
    expect_match(conditionMessage(err),
                 "'true_cv' must be a single finite number greater than 0",
                 fixed = TRUE)
  }
  expect_match(conditionMessage(expect_error(var_cv(NA_real_, 30))),
               "NA was given", fixed = TRUE)
  expect_match(conditionMessage(expect_error(var_cv(Inf, 30))),
               "Inf was given", fixed = TRUE)
})

test_that("var_cv() refuses a sample size that is not one finite number of at least 3", {
  bad <- list(NA_real_, NaN, NA, NA_integer_, Inf, -Inf, 2, 2.999, 0, -30,
              c(30, 40), numeric(0), "30", NULL, TRUE)
  for (x in bad) {
    err <- expect_error(var_cv(true_cv = 0.2, N = x),
                        class = "dmar_argument_invalid",
                        info = paste(deparse(x), collapse = ""))
    expect_identical(err$arg, "N")
    expect_match(conditionMessage(err),
                 "'N' must be a single finite sample size of at least 3",
                 fixed = TRUE)
  }
  expect_match(conditionMessage(expect_error(var_cv(0.2, Inf))),
               "Inf was given", fixed = TRUE)
  # The smallest sample size, a fractional one, and an integer are accepted.
  expect_equal(var_cv(0.2, 3)$value, 0.04 * 0.54 / 2, tolerance = 1e-14)
  expect_equal(var_cv(0.2, 30.5)$value, 0.04 * 0.54 / 29.5, tolerance = 1e-14)
  expect_identical(var_cv(0.2, 30L)$value, var_cv(0.2, 30)$value)
})

test_that("var_cv() does not accept the names its arguments had in DMAR 1.0.x", {
  # The names are data here, so that a rename sweep does not change them.
  expect_error(do.call(var_cv, list(cv = 0.2, N = 30)), "unused argument")
  expect_error(do.call(var_cv, list(true_cv = 0.2, n = 30)), "unused argument")
})
