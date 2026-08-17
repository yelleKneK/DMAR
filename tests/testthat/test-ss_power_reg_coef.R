test_that("ss_power_reg_coef() returns a tidy data.frame with the documented columns", {
  res <- ss_power_reg_coef(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_N", "actual_power") %in% res$term))
})

test_that("ss_power_reg_coef() actual_power at returned N reaches desired_power", {
  res <- ss_power_reg_coef(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5,
                           desired_power = 0.85)
  expect_gte(res$value[res$term == "actual_power"], 0.85 - 1e-6)
})

test_that("ss_power_reg_coef() effect size is sqrt((R2 - R2_without_j)/(1 - R2))", {
  rho2 <- 0.50; rho2_j <- 0.40
  res  <- ss_power_reg_coef(rho2_Y_X = rho2, rho2_Y_X_without_j = rho2_j, p = 5)
  f2 <- (rho2 - rho2_j) / (1 - rho2)
  expect_equal(res$value[res$term == "effect_size"], sqrt(f2),
               tolerance = 1e-9)
})

test_that("ss_power_reg_coef() rc and reg_coef agree when the parameterizations align", {
  rc_res    <- ss_power_rc(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5)
  coef_res  <- ss_power_reg_coef(rho2_Y_X = 0.50, rho2_Y_X_without_j = 0.40, p = 5)
  # Both functions plan for the same noncentrality, so the n's should match.
  expect_equal(rc_res$value[rc_res$term == "necessary_N"],
               coef_res$value[coef_res$term == "necessary_N"])
})

test_that("ss_power_reg_coef() Cohen's f^2 path matches the rho2 path", {
  rho2 <- 0.50; rho2_j <- 0.40
  f2 <- (rho2 - rho2_j) / (1 - rho2)
  via_cohen <- ss_power_reg_coef(cohen_f2 = f2, p = 5)
  via_rho   <- ss_power_reg_coef(rho2_Y_X = rho2, rho2_Y_X_without_j = rho2_j, p = 5)
  expect_equal(via_cohen$value[via_cohen$term == "necessary_N"],
               via_rho$value[via_rho$term == "necessary_N"])
})

test_that("ss_power_reg_coef() rejects inconsistent rho specifications", {
  expect_error(
    ss_power_reg_coef(rho2_Y_X = 0.5, rho2_Y_X_without_j = 0.4, p = 5,
                      rho_YX = 1:5, rho_XX = diag(5)),
    "do not specify"
  )
})

test_that("ss_power_reg_coef two-sided power counts both tails (HIGH-02)", {
  # A nondirectional level-alpha test must reject with probability alpha at a
  # null effect, and match the base-R noncentral-t two-tail expression at any
  # effect. The bug summed only the upper tail, returning alpha/2 at the null.
  two_tail <- function(ncp, df, a = 0.05) {
    cv <- stats::qt(1 - a / 2, df)
    stats::pt(-cv, df, ncp) + stats::pt(cv, df, ncp, lower.tail = FALSE)
  }
  for (f2 in c(0, 0.001, 0.05, 0.2)) {
    r <- ss_power_reg_coef(cohen_f2 = f2, p = 2, specified_N = 20)
    got <- r$value[r$term == "actual_power"]
    expect_equal(got, two_tail(sqrt(20) * sqrt(f2), 20 - 2 - 1), tolerance = 1e-10)
  }
  # At the null the power is exactly the nominal alpha, not alpha/2.
  r0 <- ss_power_reg_coef(cohen_f2 = 0, p = 2, specified_N = 20)
  expect_equal(r0$value[r0$term == "actual_power"], 0.05, tolerance = 1e-10)
})

test_that("ss_power_reg_coef returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_reg_coef(cohen_f2 = 0.6, p = 5, specified_N = n)
    r$value[r$term == "actual_power"] }
  target <- pw(7) - 1e-4                        # minimum admissible is p + 2 = 7 (df >= 1)
  r <- ss_power_reg_coef(cohen_f2 = 0.6, p = 5, desired_power = target)
  expect_equal(r$value[r$term == "necessary_N"], 7)
})

test_that("ss_power_reg_coef rejects invalid fixed sample sizes (MEDIUM-03)", {
  expect_error(ss_power_reg_coef(cohen_f2 = 0.2, p = 5, specified_N = 5), "whole number")
  expect_error(ss_power_reg_coef(cohen_f2 = 0.2, p = 2, specified_N = 12.5), "whole number")
})
