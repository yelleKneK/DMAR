test_that("ss_power_equivalence_c() returns documented rows", {
  res <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                delta_upper = 5, desired_power = 0.90)
  expect_setequal(res$term, c("necessary_n_per_group", "total_N", "actual_power"))
  n <- res$value[res$term == "necessary_n_per_group"]
  expect_equal(res$value[res$term == "total_N"], 2 * n)
})

test_that("ss_power_equivalence_c() recommendation is minimal", {
  res <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                delta_upper = 5, desired_power = 0.90)
  n <- res$value[res$term == "necessary_n_per_group"]
  at <- function(n) power_equivalence_c(c_weights = c(1, -1), n = n,
                                        sigma = 15.67,
                                        delta_upper = 5)$value
  expect_gte(at(n), 0.90)
  expect_lt(at(n - 1), 0.90)
  expect_equal(res$value[res$term == "actual_power"], at(n), tolerance = 1e-10)
})

test_that("ss_power_equivalence_c() noninferiority needs fewer observations
           than equivalence", {
  n_eq <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                 delta_upper = 5, desired_power = 0.90
  )$value[1]
  n_ni <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                 delta_upper = 5, desired_power = 0.90,
                                 side = "noninferiority")$value[1]
  expect_lt(n_ni, n_eq)
})

test_that("ss_power_equivalence_c() requirement grows as true_psi approaches
           a bound", {
  n_0 <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                delta_upper = 5, true_psi = 0,
                                desired_power = 0.90)$value[1]
  n_2 <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                delta_upper = 5, true_psi = 2,
                                desired_power = 0.90)$value[1]
  expect_gt(n_2, n_0)
})

test_that("ss_power_equivalence_c() is coherent with the half-width design
           rule", {
  # Targeting a 90% CI half-width of delta / 2 gives a declaration
  # probability of about .90 at psi = 0, so the power-based plan at
  # desired_power = 0.90 should land near the AIPE plan for width = delta.
  n_power <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                    delta_upper = 5,
                                    desired_power = 0.90)$value[1]
  n_aipe <- ss_aipe_c(error_variance = 15.67^2, c_weights = c(1, -1),
                      width = 5, conf_level = 0.90)$value[1]
  expect_lt(abs(n_power - n_aipe) / n_aipe, 0.05)
})

test_that("ss_power_equivalence_c() rejects an infeasible true_psi", {
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), sigma = 10,
                                      delta_upper = 5, true_psi = 5),
               "strictly inside")
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), sigma = 10,
                                      delta_upper = 5, true_psi = -6,
                                      side = "noninferiority"),
               "must exceed")
})

test_that("ss_power_equivalence_c() validates its other inputs", {
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67),
               "delta_upper")
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                                      delta_upper = 5, alpha_level = 0.6),
               "alpha")
  expect_error(ss_power_equivalence_c(c_weights = c(1, 1), sigma = 15.67,
                                      delta_upper = 5),
               "sum of the contrast weights")
  expect_error(ss_power_equivalence_c(c_weights = c(2, -2), sigma = 15.67,
                                      delta_upper = 5),
               "positive weights must sum to 1")
})

test_that("ss_power_equivalence_c() plans the same N regardless of scale", {
  # CRITICAL-01 regression, planner path. The inherited endpoint defect made
  # the resolved sample size depend on the units of the response: the same
  # design planned n = 24 per group on one scale and n = 198 on a rescaled
  # copy. The planned N and the realized power must be scale invariant.
  plan_at_scale <- function(s) {
    r <- ss_power_equivalence_c(c_weights = c(1, -1), sigma = 0.20 * s,
                                delta_upper = 0.2 * s, true_psi = 0.05 * s,
                                desired_power = 0.80)
    c(n = r$value[r$term == "necessary_n_per_group"],
      power = r$value[r$term == "actual_power"])
  }
  base <- plan_at_scale(1)
  for (s in c(1e-6, 1e-3, 1e3, 1e6)) {
    got <- plan_at_scale(s)
    expect_identical(got[["n"]], base[["n"]])
    expect_equal(got[["power"]], base[["power"]], tolerance = 1e-10)
  }
})
