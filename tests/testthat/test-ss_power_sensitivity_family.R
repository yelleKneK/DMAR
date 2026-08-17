## Monte Carlo sensitivity-analysis tests for the ss_power_* family. These
## are inherently slow; gated by skip_on_cran() to keep the CRAN test budget
## in check. The closed-form planners (test-ss_power_R2.R,
## test-ss_power_reg_coef.R) cover the same functional surface that runs on
## CRAN.

test_that("ss_power_R2_sensitivity() returns a tidy data.frame with the expected rows", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_power_R2_sensitivity(true_R2 = 0.30, estimated_R2 = 0.30,
                            desired_power = 0.80, p = 5,
                            random_predictors = TRUE, G = 100,
                            print_iter = FALSE)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("total_N", "empirical_power", "analytic_power",
                    "mean_R2", "F_crit") %in% res$term))
  # The planning inputs are echoed under their own names.
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v["true_R2"]), 0.30)
  expect_equal(unname(v["estimated_R2"]), 0.30)
  expect_equal(unname(v["p"]), 5)
  expect_equal(unname(v["desired_power"]), 0.80)
  expect_equal(unname(v["alpha_level"]), 0.05)
  # Empirical power should be within a generous MC tolerance of the analytic
  # power (G = 100, SE ~ 0.04 at power 0.8).
  emp <- res$value[res$term == "empirical_power"]
  ana <- res$value[res$term == "analytic_power"]
  expect_lt(abs(emp - ana), 0.15)
})

test_that("ss_power_R2_sensitivity() honors specified_N", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(
    ss_power_R2_sensitivity(true_R2 = 0.30, specified_N = 50,
                            p = 5, G = 100, print_iter = FALSE)
  )
  expect_equal(res$value[res$term == "total_N"], 50)
  # With a directly specified size there is no planning value to echo.
  v <- stats::setNames(res$value, res$term)
  expect_true(is.na(v["estimated_R2"]))
  expect_true(is.na(v["desired_power"]))
})

test_that("ss_power_R2_sensitivity() crossing fixed planning with random data realizes lower empirical power than analytic", {
  skip_on_cran()
  set.seed(113)
  # When planning under fixed predictors but data are random, the realized
  # power should fall short of the analytic fixed-predictor power because
  # Cohen's noncentral-F overstates random-predictor power.
  res <- suppressMessages(
    ss_power_R2_sensitivity(true_R2 = 0.30, estimated_R2 = 0.30,
                            desired_power = 0.80, p = 5,
                            random_predictors = FALSE,
                            generate_random_predictors = TRUE,
                            G = 1500, print_iter = FALSE)
  )
  emp <- res$value[res$term == "empirical_power"]
  ana <- res$value[res$term == "analytic_power"]
  expect_lt(emp, ana)
})

test_that("ss_power_reg_coef_sensitivity() returns a tidy summary frame", {
  skip_on_cran()
  set.seed(113)
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  res <- suppressMessages(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 1, desired_power = 0.80,
      G = 200, print_iter = FALSE
    )
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("total_N", "empirical_power", "analytic_power",
                    "mean_b_j", "mean_R2", "t_crit", "true_b_j") %in% res$term))
  # The planning inputs are echoed under their own names.
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v["p"]), 2)
  expect_equal(unname(v["which_predictor"]), 1)
  expect_equal(unname(v["desired_power"]), 0.80)
  expect_equal(unname(v["alpha_level"]), 0.05)
  expect_equal(unname(v["estimated_b_j"]), unname(v["true_b_j"]),
               tolerance = 1e-10)
})

test_that("ss_power_reg_coef_sensitivity() empirical power agrees with analytic at G = 1000", {
  skip_on_cran()
  set.seed(113)
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  res <- suppressMessages(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 1, desired_power = 0.80,
      G = 1000, print_iter = FALSE
    )
  )
  emp <- res$value[res$term == "empirical_power"]
  ana <- res$value[res$term == "analytic_power"]
  expect_lt(abs(emp - ana), 0.05)
})

test_that("ss_power_reg_coef_sensitivity() validates which_predictor against p", {
  # Out-of-range which_predictor previously fell through to ss_power_reg_coef
  # and triggered an infinite loop (negative-index past length silently
  # makes rho2_Y_X_without_j equal rho2_Y_X, so f = 0 and the search never
  # converges). The validation now fails fast.
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  expect_error(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 5, G = 20, print_iter = FALSE
    ),
    "which_predictor"
  )
  expect_error(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 0, G = 20, print_iter = FALSE
    ),
    "which_predictor"
  )
  expect_error(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 1.5, G = 20, print_iter = FALSE
    ),
    "which_predictor"
  )
})

test_that("ss_power_R2_sensitivity() and ss_power_reg_coef_sensitivity() reject bad scalar inputs", {
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX  <- c(0.4, 0.3)
  expect_error(
    ss_power_R2_sensitivity(true_R2 = 0.3, estimated_R2 = 0.3, p = 5,
                            desired_power = 1.5, G = 5, print_iter = FALSE),
    "desired_power"
  )
  expect_error(
    ss_power_R2_sensitivity(true_R2 = 0.3, estimated_R2 = 0.3, p = 5,
                            alpha_level = 1.5, G = 5, print_iter = FALSE),
    "alpha_level"
  )
  expect_error(
    ss_power_R2_sensitivity(true_R2 = 0.3, estimated_R2 = 0.3, p = 5,
                            G = -1, print_iter = FALSE),
    "G"
  )
  expect_error(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
      which_predictor = 1, desired_power = 2, G = 20, print_iter = FALSE
    ),
    "desired_power"
  )
})

test_that("ss_power_reg_coef_sensitivity() directional test handles negative coefficients", {
  # If True_b_j < 0, the rejection rule must use the left tail. Otherwise
  # empirical power collapses to 0 even when the analytic power is large.
  skip_on_cran()
  set.seed(113)
  Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX_neg <- c(-0.4, 0.0)
  res <- suppressMessages(
    ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = cov_YX_neg, true_cov_XX = Sigma_X,
      which_predictor = 1, desired_power = 0.80,
      directional = TRUE, G = 1000, print_iter = FALSE
    )
  )
  emp <- res$value[res$term == "empirical_power"]
  ana <- res$value[res$term == "analytic_power"]
  true_b <- res$value[res$term == "true_b_j"]
  expect_lt(true_b, 0)
  expect_lt(abs(emp - ana), 0.05)
})
