test_that("simulate_ancova_data() returns a long-format data.frame", {
  set.seed(113)
  d <- simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10,
                            true_sd_Y = 1, true_sd_X = 2,
                            true_cor = 0.8, a = 2, n_per_group = 20)
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
  d <- simulate_ancova_data(true_mean_Y = c(50, 55, 60, 65),
                            true_mean_X = c(10, 12, 11, 13),
                            true_sd_Y = 8, true_sd_X = 3,
                            true_cor = 0.3, a = 4,
                            n_per_group = c(40, 35, 45, 30),
                            randomized = FALSE)
  expect_equal(as.integer(table(d$group)), c(40L, 35L, 45L, 30L))
  expect_equal(nrow(d), 150L)
})

test_that("simulate_ancova_data() honors per-group correlations under randomized = FALSE", {
  set.seed(113)
  rho_target <- c(0.30, 0.30, 0.60, 0.60)
  d <- simulate_ancova_data(true_mean_Y = c(50, 55, 60, 65),
                            true_mean_X = c(10, 12, 11, 13),
                            true_sd_Y = 8, true_sd_X = 3,
                            true_cor = rho_target, a = 4,
                            n_per_group = rep(2000, 4),    # large n for tight verification
                            randomized = FALSE)
  rho_obs <- sapply(split(d, d$group), function(g) cor(g$y, g$x))
  expect_equal(unname(rho_obs), rho_target, tolerance = 0.05)
})

test_that("simulate_ancova_data() errors when randomized = TRUE and true_cor is a vector", {
  expect_error(
    simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10, true_sd_Y = 1, true_sd_X = 2,
                         true_cor = c(0.3, 0.7), a = 2, n_per_group = 20,
                         randomized = TRUE),
    "randomized = TRUE.*single number"
  )
})

test_that("simulate_ancova_data() errors when randomized = TRUE and true_mean_X is a vector", {
  expect_error(
    simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = c(10, 12),
                         true_sd_Y = 1, true_sd_X = 2,
                         true_cor = 0.5, a = 2, n_per_group = 20,
                         randomized = TRUE),
    "randomized = TRUE.*single number"
  )
})

test_that("simulate_ancova_data() validates argument shapes", {
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5, 7), true_mean_X = 10,
                                    true_sd_Y = 1, true_sd_X = 2,
                                    true_cor = 0.5, a = 2, n_per_group = 20),
               "true_mean_Y.*length")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = c(10, 11, 12),
                                    true_sd_Y = 1, true_sd_X = 2,
                                    true_cor = 0.5, a = 2, n_per_group = 20,
                                    randomized = FALSE),
               "true_mean_X.*length")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = c(10, 12),
                                    true_sd_Y = 1, true_sd_X = 2,
                                    true_cor = c(0.3, 0.5, 0.7),
                                    a = 2, n_per_group = 20,
                                    randomized = FALSE),
               "single number.*vector of length")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10,
                                    true_sd_Y = 1, true_sd_X = 2,
                                    true_cor = 1.2, a = 2, n_per_group = 20),
               "between -1 and 1")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10,
                                    true_sd_Y = -1, true_sd_X = 2,
                                    true_cor = 0.5, a = 2, n_per_group = 20),
               "true_sd_Y.*positive")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10,
                                    true_sd_Y = 1, true_sd_X = 2,
                                    true_cor = 0.5, a = 1, n_per_group = 20),
               "a.*>= 2")
  expect_error(simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10,
                                    true_sd_Y = 1, true_sd_X = 2,
                                    true_cor = 0.5, a = 2, n_per_group = c(10, 20, 30)),
               "single number.*vector of length")
})

test_that("simulate_ancova_data() recovers approximate population means with large n", {
  set.seed(113)
  mu_y <- c(50, 60)
  mu_x <- c(10, 12)
  d <- simulate_ancova_data(true_mean_Y = mu_y, true_mean_X = mu_x,
                            true_sd_Y = 8, true_sd_X = 3,
                            true_cor = 0.4, a = 2, n_per_group = c(5000, 5000),
                            randomized = FALSE)
  obs_my <- as.numeric(by(d$y, d$group, mean))
  obs_mx <- as.numeric(by(d$x, d$group, mean))
  expect_equal(obs_my, mu_y, tolerance = 0.2)
  expect_equal(obs_mx, mu_x, tolerance = 0.2)
})

test_that("seed reproducibility: identical seed -> identical data", {
  set.seed(113)
  d1 <- simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10, true_sd_Y = 1,
                             true_sd_X = 2, true_cor = 0.8, a = 2, n_per_group = 20)
  set.seed(113)
  d2 <- simulate_ancova_data(true_mean_Y = c(3, 5), true_mean_X = 10, true_sd_Y = 1,
                             true_sd_X = 2, true_cor = 0.8, a = 2, n_per_group = 20)
  expect_equal(d1, d2)
})
