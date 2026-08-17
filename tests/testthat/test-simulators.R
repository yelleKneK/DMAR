# =========================================================================
# simulate_anova_data
# =========================================================================
test_that("simulate_anova_data() returns a long-format data.frame", {
  set.seed(113)
  d <- simulate_anova_data(mu = c(50, 55, 60), sigma = 8, a = 3, n = 30)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 90L)
  expect_named(d, c("group", "y"))
  expect_s3_class(d$group, "factor")
})

test_that("simulate_anova_data() supports unequal n per group", {
  set.seed(113)
  d <- simulate_anova_data(mu = c(50, 55, 60), sigma = 8, a = 3,
                           n = c(40, 30, 20))
  expect_equal(as.integer(table(d$group)), c(40L, 30L, 20L))
})

test_that("simulate_anova_data() recovers approximate group means with large n", {
  set.seed(113)
  mu <- c(50, 55, 60, 65)
  d  <- simulate_anova_data(mu = mu, sigma = 8, a = 4, n = 5000)
  obs_mu <- as.numeric(by(d$y, d$group, mean))
  expect_equal(obs_mu, mu, tolerance = 0.3)
})

test_that("simulate_anova_data() supports per-group sigma (heteroscedastic)", {
  set.seed(113)
  d <- simulate_anova_data(mu = c(50, 50), sigma = c(2, 10),
                           a = 2, n = 5000)
  obs_sd <- as.numeric(by(d$y, d$group, sd))
  expect_equal(obs_sd, c(2, 10), tolerance = 0.3)
})

test_that("simulate_anova_data() rejects bad inputs", {
  expect_error(simulate_anova_data(mu = c(1, 2), sigma = 1, a = 1, n = 10),
               "a.*>= 2")
  expect_error(simulate_anova_data(mu = c(1, 2, 3), sigma = 1, a = 2, n = 10),
               "mu.*length")
  expect_error(simulate_anova_data(mu = c(1, 2), sigma = c(1, 2, 3),
                                   a = 2, n = 10),
               "sigma")
  expect_error(simulate_anova_data(mu = c(1, 2), sigma = -1, a = 2, n = 10),
               "single positive number")
})

test_that("simulate_anova_data() seed gives reproducible output", {
  d1 <- simulate_anova_data(mu = c(1, 2), sigma = 1, a = 2, n = 30, seed = 113)
  d2 <- simulate_anova_data(mu = c(1, 2), sigma = 1, a = 2, n = 30, seed = 113)
  expect_equal(d1, d2)
})


# =========================================================================
# simulate_regression_data
# =========================================================================
test_that("simulate_regression_data() returns a data.frame with the documented shape", {
  set.seed(113)
  d <- simulate_regression_data(N = 200, p = 5, rho_YX = rep(0.30, 5))
  expect_s3_class(d, "data.frame")
  expect_equal(dim(d), c(200L, 6L))
  expect_named(d, c("y", paste0("x", 1:5)))
})

test_that("simulate_regression_data() recovers rho_YX with large N", {
  set.seed(113)
  rho_target <- c(.50, .40, .30, .20, .10)
  rho_XX <- diag(5)
  d <- simulate_regression_data(N = 50000, p = 5, rho_YX = rho_target,
                                rho_XX = rho_XX)
  rho_obs <- as.numeric(cor(d$y, d[, paste0("x", 1:5)]))
  expect_equal(rho_obs, rho_target, tolerance = 0.02)
})

test_that("simulate_regression_data() recovers an exchangeable rho_XX", {
  set.seed(113)
  RHO <- matrix(0.5, 5, 5); diag(RHO) <- 1
  d <- simulate_regression_data(N = 20000, p = 5, rho_YX = rep(0.3, 5),
                                rho_XX = RHO)
  XX <- d[, paste0("x", 1:5)]
  off <- cor(XX)[lower.tri(diag(5))]
  expect_equal(mean(off), 0.5, tolerance = 0.02)
})

test_that("simulate_regression_data() rejects an inconsistent joint correlation", {
  # Y correlates .9 with each X, but the X's correlate -.9 with each other.
  # Implied joint correlation matrix is not positive definite.
  bad_RHO <- matrix(c( 1, -0.9,
                      -0.9,  1), nrow = 2)
  expect_error(simulate_regression_data(N = 100, p = 2,
                                        rho_YX = c(0.9, 0.9),
                                        rho_XX = bad_RHO),
               "positive definite")
})

test_that("simulate_regression_data() rejects bad inputs", {
  # N must be at least p + 2 (so residual df is positive after fitting).
  expect_error(simulate_regression_data(N = 6, p = 5, rho_YX = rep(.3, 5)),
               "N.*> p \\+ 1")
  expect_error(simulate_regression_data(N = 100, p = 5, rho_YX = c(.3, .3)),
               "rho_YX.*length")
  expect_error(simulate_regression_data(N = 100, p = 2, rho_YX = c(1.5, .3)),
               "rho_YX.*\\(-1, 1\\)")
})


# =========================================================================
# simulate_ancova_factorial_data
# =========================================================================
test_that("simulate_ancova_factorial_data() 1-factor case matches simulate_ancova_data()", {
  set.seed(113)
  d_simple <- simulate_ancova_data(mu_y = c(50, 55), mu_x = 10,
                                   sigma_y = 8, sigma_x = 3,
                                   rho = 0.4, a = 2, n = 30)
  set.seed(113)
  d_fact <- simulate_ancova_factorial_data(
    a = 2, mu_y = c(50, 55),
    mu_x = matrix(10, nrow = 2, ncol = 1),
    sigma_y = 8, sigma_x = 3, rho_y_x = 0.4, n = 30
  )
  # Same seed and same parameters -> identical Y values (column orders differ)
  expect_equal(d_simple$y, d_fact$y, tolerance = 1e-12)
})

test_that("simulate_ancova_factorial_data() 2x2 design produces the right cell counts", {
  set.seed(113)
  d <- simulate_ancova_factorial_data(
    a = 2, b = 2, mu_y = c(50, 60, 55, 65),
    mu_x = matrix(10, nrow = 4, ncol = 1),
    sigma_y = 8, sigma_x = 3, rho_y_x = 0.4, n = 30
  )
  expect_named(d, c("A", "B", "x1", "y"))
  expect_equal(nrow(d), 120L)
  expect_equal(as.integer(table(d$A, d$B)), c(30L, 30L, 30L, 30L))
})

test_that("simulate_ancova_factorial_data() 4-way 2x2x2x2 design works", {
  set.seed(113)
  mu_y <- 50 + (0:15) * 0.5
  d <- simulate_ancova_factorial_data(
    a = 2, b = 2, c = 2, d = 2,
    mu_y = mu_y,
    mu_x = matrix(10, nrow = 16, ncol = 1),
    sigma_y = 8, sigma_x = 3, rho_y_x = 0.3, n = 10
  )
  expect_equal(nrow(d), 160L)
  expect_named(d, c("A", "B", "C", "D", "x1", "y"))
  expect_equal(nlevels(d$A), 2L)
  expect_equal(nlevels(d$D), 2L)
})

test_that("simulate_ancova_factorial_data() supports multiple covariates", {
  set.seed(113)
  d <- simulate_ancova_factorial_data(
    a = 3, b = 2, n_covariates = 3,
    mu_y = c(50, 55, 60, 52, 58, 64),
    mu_x = matrix(c(10, 11, 12, 9, 10, 11,
                    5,  6,  7, 4,  5,  6,
                    1,  2,  3, 0,  1,  2),
                  nrow = 6, ncol = 3),
    sigma_y = 8, sigma_x = c(3, 2, 1),
    rho_y_x = c(0.4, 0.3, 0.2),
    rho_x_x = diag(3),
    n = 30, randomized = FALSE
  )
  expect_named(d, c("A", "B", "x1", "x2", "x3", "y"))
  expect_equal(nrow(d), 6 * 30L)
})

test_that("simulate_ancova_factorial_data() errors when randomized = TRUE and mu_x varies", {
  expect_error(
    simulate_ancova_factorial_data(
      a = 2, b = 2, mu_y = c(50, 55, 60, 65),
      mu_x = matrix(c(10, 11, 12, 13), nrow = 4, ncol = 1),
      sigma_y = 8, sigma_x = 3, rho_y_x = 0.4, n = 30,
      randomized = TRUE),
    "randomized = TRUE.*every column.*constant"
  )
})

test_that("simulate_ancova_factorial_data() recovers cell means with large n", {
  set.seed(113)
  cell_mu_y <- c(50, 55, 60, 65)
  d <- simulate_ancova_factorial_data(
    a = 2, b = 2, mu_y = cell_mu_y,
    mu_x = matrix(10, nrow = 4, ncol = 1),
    sigma_y = 8, sigma_x = 3, rho_y_x = 0.4, n = 2000
  )
  obs <- aggregate(y ~ A + B, data = d, FUN = mean)
  obs_ordered <- obs$y[order(obs$A, obs$B)]
  target_ordered <- cell_mu_y[order(rep(1:2, 2), rep(1:2, each = 2))]
  expect_equal(obs_ordered, target_ordered, tolerance = 0.3)
})

test_that("simulate_ancova_factorial_data() supports unequal n per cell", {
  set.seed(113)
  d <- simulate_ancova_factorial_data(
    a = 2, b = 2, mu_y = c(50, 55, 60, 65),
    mu_x = matrix(10, nrow = 4, ncol = 1),
    sigma_y = 8, sigma_x = 3, rho_y_x = 0.4,
    n = c(20, 25, 30, 35)
  )
  expect_equal(nrow(d), 110L)
})

test_that("simulate_ancova_factorial_data() rejects inconsistent rho_y_x and rho_x_x", {
  # Y near-perfectly correlated with both x1 and x2, but x1 nearly uncorrelated
  # with x2 -> implied joint matrix not positive-definite.
  expect_error(
    simulate_ancova_factorial_data(
      a = 2, n_covariates = 2,
      mu_y = c(50, 55),
      mu_x = matrix(c(10, 10, 5, 5), nrow = 2, ncol = 2),
      sigma_y = 8, sigma_x = c(3, 2),
      rho_y_x = c(0.95, 0.95),
      rho_x_x = matrix(c(1, 0, 0, 1), 2, 2),
      n = 30),
    "positive definite"
  )
})
