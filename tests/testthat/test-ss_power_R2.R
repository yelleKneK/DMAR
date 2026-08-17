test_that("ss_power_R2() solves for the necessary sample size", {
  result <- ss_power_R2(population_R2 = 0.3, p = 5, alpha_level = 0.05, desired_power = 0.8)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_true("necessary_N" %in% result$term)
  expect_true("actual_power" %in% result$term)
  expect_gt(result$value[result$term == "necessary_N"], 0)
  expect_gte(result$value[result$term == "actual_power"], 0.8)
})

test_that("ss_power_R2() with specified_N returns realized power", {
  result <- ss_power_R2(population_R2 = 0.3, p = 5, alpha_level = 0.05, specified_N = 100)
  expect_true("specified_N" %in% result$term)
  expect_equal(result$value[result$term == "specified_N"], 100)
  expect_gt(result$value[result$term == "actual_power"], 0.99)
})

test_that("ss_power_R2() requires fewer subjects for a larger effect", {
  small <- ss_power_R2(population_R2 = 0.10, p = 5, alpha_level = 0.05, desired_power = 0.8)$value[1]
  large <- ss_power_R2(population_R2 = 0.30, p = 5, alpha_level = 0.05, desired_power = 0.8)$value[1]
  expect_gt(small, large)
})

test_that("ss_power_R2() random predictors requires at least as many subjects as fixed predictors", {
  # Cohen (1988) noncentral-F overstates power under random predictors and so
  # understates required N. The Lee (1971) random-predictor path corrects this,
  # so N_random >= N_fixed across plausible scenarios.
  scenarios <- expand.grid(rho2 = c(0.10, 0.25, 0.50),
                           p    = c(3, 5),
                           pow  = c(0.80, 0.90))
  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    n_fixed  <- ss_power_R2(population_R2 = s$rho2, p = s$p, desired_power = s$pow,
                            random_predictors = FALSE)$value[1]
    n_random <- ss_power_R2(population_R2 = s$rho2, p = s$p, desired_power = s$pow,
                            random_predictors = TRUE)$value[1]
    expect_gte(n_random, n_fixed)
  }
})

test_that("ss_power_R2() random predictors at specified N has noncentral_f_parm NA", {
  res <- ss_power_R2(population_R2 = 0.3, p = 5, specified_N = 50,
                     random_predictors = TRUE)
  expect_true(is.na(res$value[res$term == "noncentral_f_parm"]))
})

test_that("ss_power_R2() respects random_predictors", {
  # The fixed-predictors and random-predictors paths give different
  # noncentralities at the same population R^2 and N; verify the
  # function actually branches.
  r_fixed  <- ss_power_R2(population_R2 = 0.3, p = 5, specified_N = 50,
                          random_predictors = FALSE)$value[2]
  r_random <- ss_power_R2(population_R2 = 0.3, p = 5, specified_N = 50,
                          random_predictors = TRUE)$value[2]
  expect_false(isTRUE(all.equal(r_fixed, r_random, tolerance = 1e-3)))
})

test_that("ss_power_R2() cohen_f2 back-derivation respects null_R2", {
  # Math: cohen_f2 = (population_R2 - null_R2) / (1 - population_R2)
  # so    population_R2 = (cohen_f2 + null_R2) / (1 + cohen_f2).
  # The two entry points must agree under both fixed and random predictors.
  null_R2  <- 0.10
  true_R2  <- 0.30
  cohen_f2 <- (true_R2 - null_R2) / (1 - true_R2)

  for (rp in c(FALSE, TRUE)) {
    n_via_R2 <- ss_power_R2(population_R2 = true_R2, null_R2 = null_R2, p = 5,
                            desired_power = 0.80, random_predictors = rp)$value[1]
    n_via_f2 <- ss_power_R2(cohen_f2 = cohen_f2, null_R2 = null_R2, p = 5,
                            desired_power = 0.80, random_predictors = rp)$value[1]
    expect_equal(n_via_R2, n_via_f2)
  }
})

test_that("ss_power_R2() power at population_R2 <= null_R2 equals alpha_level (both paths)", {
  # The right-tail F-test of H0: rho^2 = null_R2 vs H1: rho^2 > null_R2
  # rejects at rate alpha_level when the truth is at or below the null.
  # Previously the fixed path produced NaNs because pf() got a negative ncp;
  # the random path was fine. The shared guard now hoists this above the
  # branch.
  for (rp in c(FALSE, TRUE)) {
    res_at <- ss_power_R2(population_R2 = 0.20, null_R2 = 0.20, p = 5,
                          specified_N = 200, random_predictors = rp)
    expect_equal(res_at$value[res_at$term == "actual_power"], 0.05,
                 tolerance = 1e-6)
    res_below <- ss_power_R2(population_R2 = 0.05, null_R2 = 0.20, p = 5,
                             specified_N = 200, random_predictors = rp)
    expect_equal(res_below$value[res_below$term == "actual_power"], 0.05,
                 tolerance = 1e-6)
  }
})

test_that("ss_power_R2() validates desired_power on the cohen_f2 entry point", {
  # Previously, desired_power validation was inside the population_R2 branch,
  # so an invalid desired_power supplied alongside cohen_f2 alone would drive
  # the iterative search into an infinite loop.
  expect_error(
    ss_power_R2(cohen_f2 = 1, desired_power = 1.5, p = 5,
                random_predictors = FALSE),
    "desired_power"
  )
  expect_error(
    ss_power_R2(cohen_f2 = 1, desired_power = -0.5, p = 5,
                random_predictors = TRUE),
    "desired_power"
  )
})

test_that("ss_power_R2() rejects a specified_N with no residual degrees of freedom (MEDIUM-03)", {
  # df_2 = N - p - 1 = 0 at specified_N = p + 1: the F test is undefined, so the
  # function errors rather than returning a numeric-looking power of 0.
  expect_error(ss_power_R2(population_R2 = 0.30, p = 5, specified_N = 6,
                           random_predictors = TRUE), "whole number")
  expect_error(ss_power_R2(population_R2 = 0.30, p = 5, specified_N = 6,
                           random_predictors = FALSE), "whole number")
})

test_that("ss_power_R2() random-predictor power matches a Monte Carlo simulation", {
  skip_on_cran()
  # The random-predictor path uses the Lee (1971) approximation to the sampling
  # distribution of R^2 when the predictors are themselves a draw from a joint
  # multivariate normal. This pins the analytic power to ground truth: generate
  # data with the predictors RE-DRAWN each replication (the random-predictor
  # model), reject by the usual overall F-test, and compare the empirical
  # rejection rate to ss_power_R2()'s analytic power.
  set.seed(113)
  p <- 3; rho2 <- 0.20; N <- 60; alpha <- 0.05; G <- 4000

  analytic <- ss_power_R2(population_R2 = rho2, p = p, specified_N = N,
                          random_predictors = TRUE,
                          alpha_level = alpha)$value[2]

  # Equal coefficients on independent standard-normal predictors give a
  # population R^2 of exactly rho2 when the error variance is 1.
  b <- sqrt(rho2 / (p * (1 - rho2)))
  reject <- logical(G)
  for (g in seq_len(G)) {
    X    <- matrix(rnorm(N * p), N, p)              # predictors random each rep
    y    <- as.numeric(X %*% rep(b, p)) + rnorm(N)
    fobs <- summary(lm(y ~ X))$fstatistic
    reject[g] <- pf(fobs[1], fobs[2], fobs[3], lower.tail = FALSE) < alpha
  }
  empirical <- mean(reject)

  # The Lee approximation is excellent here; the Monte Carlo standard error is
  # about 0.005, so a 0.02 band is comfortable and still a real test.
  expect_lt(abs(analytic - empirical), 0.02)

  # And the fixed-predictor (Cohen, 1988) formula over-states the random-design
  # power, the discrepancy that motivates the random-predictor default.
  fixed <- ss_power_R2(population_R2 = rho2, p = p, specified_N = N,
                       random_predictors = FALSE, alpha_level = alpha)$value[2]
  expect_gt(fixed, analytic)
})

test_that("ss_power_R2 returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_R2(population_R2 = 0.3, p = 5, specified_N = n)
    r$value[r$term == "actual_power"] }
  target <- pw(7) - 1e-4                        # minimum admissible is p + 2 = 7 (df_2 >= 1)
  r <- ss_power_R2(population_R2 = 0.3, p = 5, desired_power = target)
  expect_equal(r$value[r$term == "necessary_N"], 7)
})

test_that("ss_power_R2 rejects invalid fixed sample sizes (MEDIUM-03)", {
  expect_error(ss_power_R2(population_R2 = 0.3, p = 5, specified_N = 5), "whole number")   # no residual df
  expect_error(ss_power_R2(population_R2 = 0.3, p = 2, specified_N = 20.5), "whole number") # fractional
  expect_silent(ss_power_R2(population_R2 = 0.3, p = 2, specified_N = 20))
})
