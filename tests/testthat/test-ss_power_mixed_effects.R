## Tests for ss_power_mixed_effects(), the two-level (random-intercept) power
## planner for a treatment effect at the higher level.

test_that("ss_power_mixed_effects() necessary_J grows with ICC", {
  low_icc  <- ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.01, desired_power = 0.80)
  high_icc <- ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.20, desired_power = 0.80)
  expect_gt(high_icc[high_icc$term == "necessary_J_per_arm", "value"],
            low_icc[low_icc$term   == "necessary_J_per_arm", "value"])
})

test_that("ss_power_mixed_effects() necessary_J shrinks with larger d", {
  small <- ss_power_mixed_effects(d = 0.20, n = 30, rho = 0.05, desired_power = 0.80)
  large <- ss_power_mixed_effects(d = 0.50, n = 30, rho = 0.05, desired_power = 0.80)
  expect_gt(small[small$term == "necessary_J_per_arm", "value"],
            large[large$term == "necessary_J_per_arm", "value"])
})

test_that("ss_power_mixed_effects() with rho = 0 reduces toward an unclustered two-sample t-test", {
  # With rho = 0, the design effect 1 + (n-1)*rho = 1, so each cluster contributes n
  # independent observations. The total sample per arm is J*n, and the test reduces to
  # a standard two-sample t-test on J*n per arm.
  cluster_res  <- ss_power_mixed_effects(d = 0.30, n = 5, rho = 0, desired_power = 0.80)
  J_required <- cluster_res[cluster_res$term == "necessary_J_per_arm", "value"]
  total_n_per_arm <- J_required * 5
  smd_res <- ss_power_smd(smd = 0.30, n_1 = total_n_per_arm)
  realized_power <- smd_res[smd_res$term == "actual_power", "value"]
  # The cluster-level test uses df = 2J - 2 instead of 2*Jn - 2, so it's a bit more conservative.
  # Realized power from the smd unclustered test should at least equal the desired_power level.
  expect_gte(realized_power, 0.80)
})

test_that("ss_power_mixed_effects() power at returned J equals or exceeds desired_power", {
  res <- ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, desired_power = 0.80)
  J   <- res[res$term == "necessary_J_per_arm", "value"]
  pwr <- ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, J = J)
  expect_gte(pwr[pwr$term == "actual_power", "value"], 0.80)
})

test_that("ss_power_mixed_effects() validates inputs", {
  expect_error(ss_power_mixed_effects(d = 0,    n = 30, rho = 0.05, desired_power = 0.80),
               "'d' must be a single finite nonzero numeric")
  expect_error(ss_power_mixed_effects(d = 0.30, n = 30, rho = 1.0,  desired_power = 0.80),
               "'rho' must be a single numeric value in")
  expect_error(ss_power_mixed_effects(d = 0.30, n = 30, rho = -0.1, desired_power = 0.80),
               "'rho' must be a single numeric value in")
})
