test_that("simulate_ancova_data() returns a long-format data.frame", {
  set.seed(113)
  d <- simulate_ancova_data(mu_y = c(3, 5), mu_x = 10,
                            sigma_y = 1, sigma_x = 2,
                            rho = 0.8, a = 2, n = 20)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 40L)
  expect_named(d, c("group", "y", "x"))
  expect_s3_class(d$group, "factor")
  expect_equal(levels(d$group), c("1", "2"))
  expect_type(d$y, "double")
  expect_type(d$x, "double")
})

test_that("simulate_ancova_data() supports unequal sample sizes", {
  set.seed(113)
  d <- simulate_ancova_data(mu_y = c(50, 55, 60, 65),
                            mu_x = c(10, 12, 11, 13),
                            sigma_y = 8, sigma_x = 3,
                            rho = 0.3, a = 4,
                            n = c(40, 35, 45, 30),
                            randomized = FALSE)
  expect_equal(as.integer(table(d$group)), c(40L, 35L, 45L, 30L))
  expect_equal(nrow(d), 150L)
})

test_that("simulate_ancova_data() honors per-group correlations under randomized = FALSE", {
  set.seed(113)
  rho_target <- c(0.30, 0.30, 0.60, 0.60)
  d <- simulate_ancova_data(mu_y = c(50, 55, 60, 65),
                            mu_x = c(10, 12, 11, 13),
                            sigma_y = 8, sigma_x = 3,
                            rho = rho_target, a = 4,
                            n = rep(2000, 4),    # large n for tight verification
                            randomized = FALSE)
  rho_obs <- sapply(split(d, d$group), function(g) cor(g$y, g$x))
  expect_equal(unname(rho_obs), rho_target, tolerance = 0.05)
})

test_that("simulate_ancova_data() errors when randomized = TRUE and rho is a vector", {
  expect_error(
    simulate_ancova_data(mu_y = c(3, 5), mu_x = 10, sigma_y = 1, sigma_x = 2,
                         rho = c(0.3, 0.7), a = 2, n = 20,
                         randomized = TRUE),
    "randomized = TRUE.*single number"
  )
})

test_that("simulate_ancova_data() errors when randomized = TRUE and mu_x is a vector", {
  expect_error(
    simulate_ancova_data(mu_y = c(3, 5), mu_x = c(10, 12),
                         sigma_y = 1, sigma_x = 2,
                         rho = 0.5, a = 2, n = 20,
                         randomized = TRUE),
    "randomized = TRUE.*single number"
  )
})

test_that("simulate_ancova_data() validates argument shapes", {
  expect_error(simulate_ancova_data(mu_y = c(3, 5, 7), mu_x = 10,
                                    sigma_y = 1, sigma_x = 2,
                                    rho = 0.5, a = 2, n = 20),
               "mu_y.*length")
  expect_error(simulate_ancova_data(mu_y = c(3, 5), mu_x = c(10, 11, 12),
                                    sigma_y = 1, sigma_x = 2,
                                    rho = 0.5, a = 2, n = 20,
                                    randomized = FALSE),
               "mu_x.*length")
  expect_error(simulate_ancova_data(mu_y = c(3, 5), mu_x = c(10, 12),
                                    sigma_y = 1, sigma_x = 2,
                                    rho = c(0.3, 0.5, 0.7),
                                    a = 2, n = 20,
                                    randomized = FALSE),
               "single number.*vector of length")
  expect_error(simulate_ancova_data(mu_y = c(3, 5), mu_x = 10,
                                    sigma_y = 1, sigma_x = 2,
                                    rho = 1.2, a = 2, n = 20),
               "between -1 and 1")
  expect_error(simulate_ancova_data(mu_y = c(3, 5), mu_x = 10,
                                    sigma_y = -1, sigma_x = 2,
                                    rho = 0.5, a = 2, n = 20),
               "sigma_y.*positive")
  expect_error(simulate_ancova_data(mu_y = c(3, 5), mu_x = 10,
                                    sigma_y = 1, sigma_x = 2,
                                    rho = 0.5, a = 1, n = 20),
               "a.*>= 2")
  expect_error(simulate_ancova_data(mu_y = c(3, 5), mu_x = 10,
                                    sigma_y = 1, sigma_x = 2,
                                    rho = 0.5, a = 2, n = c(10, 20, 30)),
               "single number.*vector of length")
})

test_that("simulate_ancova_data() recovers approximate population means with large n", {
  set.seed(113)
  mu_y <- c(50, 60)
  mu_x <- c(10, 12)
  d <- simulate_ancova_data(mu_y = mu_y, mu_x = mu_x,
                            sigma_y = 8, sigma_x = 3,
                            rho = 0.4, a = 2, n = c(5000, 5000),
                            randomized = FALSE)
  obs_my <- as.numeric(by(d$y, d$group, mean))
  obs_mx <- as.numeric(by(d$x, d$group, mean))
  expect_equal(obs_my, mu_y, tolerance = 0.2)
  expect_equal(obs_mx, mu_x, tolerance = 0.2)
})

test_that("seed reproducibility: identical seed -> identical data", {
  set.seed(113)
  d1 <- simulate_ancova_data(mu_y = c(3, 5), mu_x = 10, sigma_y = 1,
                             sigma_x = 2, rho = 0.8, a = 2, n = 20)
  set.seed(113)
  d2 <- simulate_ancova_data(mu_y = c(3, 5), mu_x = 10, sigma_y = 1,
                             sigma_x = 2, rho = 0.8, a = 2, n = 20)
  expect_equal(d1, d2)
})
