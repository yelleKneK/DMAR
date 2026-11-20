test_that("ss_power_equivalence_c() returns documented rows", {
  res <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                contrast_upper = 5, desired_power = 0.90)
  expect_setequal(res$term, c("necessary_n_per_group", "total_N", "actual_power"))
  n <- res$value[res$term == "necessary_n_per_group"]
  expect_equal(res$value[res$term == "total_N"], 2 * n)
})

test_that("ss_power_equivalence_c() recommendation is minimal", {
  res <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                contrast_upper = 5, desired_power = 0.90)
  n <- res$value[res$term == "necessary_n_per_group"]
  at <- function(n) power_equivalence_c(c_weights = c(1, -1), n_per_group = n,
                                        true_sd = 15.67,
                                        contrast_upper = 5)$value
  expect_gte(at(n), 0.90)
  expect_lt(at(n - 1), 0.90)
  expect_equal(res$value[res$term == "actual_power"], at(n), tolerance = 1e-10)
})

test_that("ss_power_equivalence_c() noninferiority needs fewer observations
           than equivalence", {
  n_eq <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                 contrast_upper = 5, desired_power = 0.90
  )$value[1]
  n_ni <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                 contrast_upper = 5, desired_power = 0.90,
                                 side = "noninferiority")$value[1]
  expect_lt(n_ni, n_eq)
})

test_that("ss_power_equivalence_c() requirement grows as true_contrast approaches
           a bound", {
  n_0 <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                contrast_upper = 5, true_contrast = 0,
                                desired_power = 0.90)$value[1]
  n_2 <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                contrast_upper = 5, true_contrast = 2,
                                desired_power = 0.90)$value[1]
  expect_gt(n_2, n_0)
})

test_that("ss_power_equivalence_c() is coherent with the half-width design
           rule", {
  # Targeting a 90% CI half-width of delta / 2 gives a declaration
  # probability of about .90 at true_contrast = 0, so the power-based plan at
  # desired_power = 0.90 should land near the AIPE plan for width = delta.
  n_power <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                    contrast_upper = 5,
                                    desired_power = 0.90)$value[1]
  n_aipe <- ss_aipe_c(true_var_error = 15.67^2, c_weights = c(1, -1),
                      width = 5, conf_level = 0.90)$value[1]
  expect_lt(abs(n_power - n_aipe) / n_aipe, 0.05)
})

test_that("ss_power_equivalence_c() rejects an infeasible true_contrast", {
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 10,
                                      contrast_upper = 5, true_contrast = 5),
               "strictly inside")
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 10,
                                      contrast_upper = 5, true_contrast = -6,
                                      side = "noninferiority"),
               "must exceed")
})

test_that("ss_power_equivalence_c() validates its other inputs", {
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67),
               "contrast_upper")
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                      contrast_upper = 5, alpha_level = 0.6),
               "alpha")
  expect_error(ss_power_equivalence_c(c_weights = c(1, 1), true_sd = 15.67,
                                      contrast_upper = 5),
               "sum of the contrast weights")
  expect_error(ss_power_equivalence_c(c_weights = c(2, -2), true_sd = 15.67,
                                      contrast_upper = 5),
               "positive weights must sum to 1")
})

test_that("ss_power_equivalence_c() plans the same N regardless of scale", {
  # CRITICAL-01 regression, planner path. The inherited endpoint defect made
  # the resolved sample size depend on the units of the response: the same
  # design planned n = 24 per group on one scale and n = 198 on a rescaled
  # copy. The planned N and the realized power must be scale invariant.
  plan_at_scale <- function(s) {
    r <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 0.20 * s,
                                contrast_upper = 0.2 * s, true_contrast = 0.05 * s,
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

test_that("ss_power_equivalence_c() plans noninferiority without contrast_upper", {
  # The call DMAR 1.0.x refused for want of the upper bound.
  lower_only <- ss_power_equivalence_c(c_weights = c(-1, 1), true_sd = 1,
                                       contrast_lower = .5, desired_power = .90,
                                       side = "noninferiority")
  # Given alone, as 1.0.x required, the upper bound is the margin.
  upper_only <- ss_power_equivalence_c(c_weights = c(-1, 1), true_sd = 1,
                                       contrast_upper = .5, desired_power = .90,
                                       side = "noninferiority")
  expect_identical(lower_only$value, upper_only$value)
  expect_identical(lower_only$term,
                   c("necessary_n_per_group", "total_N", "actual_power"))
  n <- lower_only$value[lower_only$term == "necessary_n_per_group"]
  at <- function(n) power_equivalence_c(c_weights = c(-1, 1), n_per_group = n,
                                        true_sd = 1, contrast_lower = .5,
                                        side = "noninferiority")$value
  expect_gte(at(n), .90)
  expect_lt(at(n - 1), .90)
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 1,
                                      side = "noninferiority"),
               "'contrast_lower' must be specified")
})

test_that("ss_power_equivalence_c() plans a paired design in units", {
  res <- ss_power_equivalence_c(true_sd = 1, contrast_lower = .5,
                                desired_power = .90, side = "noninferiority",
                                design = "paired")
  # Each condition is observed on every unit, so the answer is the number
  # of units, the total sample size: the row is necessary_N, and there is
  # no per-group size and no second total beside it.
  expect_identical(res$term, c("necessary_N", "actual_power"))
  expect_type(res$value, "double")
  expect_s3_class(res, "dmar_ss_power")
  N <- res$value[res$term == "necessary_N"]

  # Checked against the noncentral t directly: 36 units reach .90 and 35
  # do not.
  pw <- function(N) pt(qt(.95, N - 1), N - 1, ncp = .5 * sqrt(N),
                       lower.tail = FALSE)
  expect_equal(N, 36)
  expect_gte(pw(36), .90)
  expect_lt(pw(35), .90)
  expect_equal(res$value[res$term == "actual_power"], pw(36),
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "actual_power"],
               power_equivalence_c(N = N, true_sd = 1, contrast_lower = .5,
                                   side = "noninferiority",
                                   design = "paired")$value,
               tolerance = 1e-12)

  # tidy() and glance() report the planned number of units.
  expect_equal(generics::tidy(res)$estimate, 36)
  expect_equal(generics::tidy(res)$power,
               res$value[res$term == "actual_power"])
  g <- generics::glance(res)
  expect_equal(g$estimate, 36)
  expect_identical(names(g), c("term", "estimate", "power"))

  # The paired equivalence plan is minimal as well.
  eq <- ss_power_equivalence_c(true_sd = 1, contrast_upper = .5, true_contrast = .1,
                               desired_power = .80, design = "paired")
  expect_identical(eq$term, c("necessary_N", "actual_power"))
  N_eq <- eq$value[eq$term == "necessary_N"]
  at <- function(N) power_equivalence_c(N = N, true_sd = 1, contrast_upper = .5,
                                        true_contrast = .1,
                                        design = "paired")$value
  expect_gte(at(N_eq), .80)
  expect_lt(at(N_eq - 1), .80)
})

test_that("ss_power_equivalence_c() paired plan is coherent with the independent one", {
  # Uncorrelated measures with SD 15.67 give difference scores with SD
  # 15.67 * sqrt(2), the standard error of two independent groups of the
  # same size, so the paired plan in units lands at the independent plan
  # per group, moved up by at most a unit or two for its n - 1 degrees of
  # freedom against 2n - 2.
  n_ind <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                                  contrast_upper = 5, desired_power = .90)$value[1]
  N_pair <- ss_power_equivalence_c(true_sd = 15.67 * sqrt(2), contrast_upper = 5,
                                   desired_power = .90,
                                   design = "paired")$value[1]
  expect_gte(N_pair, n_ind)
  expect_lte(N_pair, n_ind + 2)
  # A positive correlation shrinks the difference SD and the plan with it.
  N_rho <- ss_power_equivalence_c(true_sd = 15.67 * sqrt(2 * (1 - .6)),
                                  contrast_upper = 5, desired_power = .90,
                                  design = "paired")$value[1]
  expect_lt(N_rho, N_pair)
})

test_that("ss_power_equivalence_c() names each size row for the quantity it holds", {
  # Independent groups are sized per group, with the total beside the
  # answer; a paired design is sized in units, which is the total sample
  # size N. The class, the columns, and the shape of tidy() are the same
  # for both designs, and each answer is the value power_equivalence_c()
  # takes under the same name.
  for (side in c("equivalence", "noninferiority")) {
    ind <- ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 1,
                                  contrast_upper = .5, desired_power = .80,
                                  side = side)
    pai <- ss_power_equivalence_c(true_sd = 1, contrast_upper = .5,
                                  desired_power = .80, side = side,
                                  design = "paired")
    expect_identical(ind$term,
                     c("necessary_n_per_group", "total_N", "actual_power"))
    expect_identical(pai$term, c("necessary_N", "actual_power"))
    expect_identical(class(pai), class(ind))
    expect_identical(names(pai), names(ind))
    expect_identical(names(generics::tidy(pai)), names(generics::tidy(ind)))
    expect_equal(ind$value[ind$term == "total_N"],
                 2 * ind$value[ind$term == "necessary_n_per_group"])
    expect_equal(
      power_equivalence_c(c_weights = c(1, -1),
                          n_per_group = ind$value[ind$term == "necessary_n_per_group"],
                          true_sd = 1, contrast_upper = .5, side = side)$value,
      ind$value[ind$term == "actual_power"], tolerance = 1e-12)
    expect_equal(
      power_equivalence_c(N = pai$value[pai$term == "necessary_N"],
                          true_sd = 1, contrast_upper = .5, side = side,
                          design = "paired")$value,
      pai$value[pai$term == "actual_power"], tolerance = 1e-12)
  }
})

test_that("ss_power_equivalence_c() validates the design", {
  expect_error(ss_power_equivalence_c(true_sd = 1, contrast_upper = 1),
               "'c_weights' must be supplied")
  expect_error(ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 1,
                                      contrast_upper = 1, design = "paired"),
               "Leave 'c_weights' NULL")
  # Naming the default design changes nothing.
  expect_identical(
    ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                           contrast_upper = 5, desired_power = 0.90,
                           design = "independent")$value,
    ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
                           contrast_upper = 5, desired_power = 0.90)$value)
})
