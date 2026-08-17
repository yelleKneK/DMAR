test_that("reproduces the four-group example: ncp = sqrt(8.4375), power ~ 0.818 at n = 20", {
  res <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    n_per_group   = 20
  )
  ncp_value   <- res$value[res$term == "noncentral_t_parm"]
  power_value <- res$value[res$term == "actual_power"]
  total_value <- res$value[res$term == "total_N"]
  expect_equal(ncp_value,   sqrt(8.4375), tolerance = 1e-6)
  expect_equal(total_value, 80)
  expect_equal(power_value, 0.818031,    tolerance = 1e-4)
})

test_that("matches the F-distribution-based reference power exactly", {
  # 1 - pf(qf(.95, 1, 76), 1, 76, ncp = 8.4375)
  expected <- 1 - pf(qf(.95, 1, 76), 1, 76, ncp = 8.4375)
  res <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    n_per_group   = 20
  )
  power_value <- res$value[res$term == "actual_power"]
  expect_equal(power_value, expected, tolerance = 1e-6)
})

test_that("solving for n returns the smallest n that meets desired_power", {
  res <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    desired_power = 0.90
  )
  n_solved <- res$value[res$term == "necessary_n_per_group"]
  power_at_n <- res$value[res$term == "actual_power"]
  expect_gte(power_at_n, 0.90)
  # Power at (n - 1) should be strictly less than .90.
  res_below <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    n_per_group   = n_solved - 1
  )
  expect_lt(res_below$value[res_below$term == "actual_power"], 0.90)
})

test_that("supplying psi directly is equivalent to supplying mu when psi = c'mu", {
  via_mu <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    n_per_group   = 20
  )
  via_Psi <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    psi           = 9,
    sigma_squared = 144,
    n_per_group   = 20
  )
  expect_equal(via_mu, via_Psi)
})

test_that("Cohen's f matches |ncp| / sqrt(N)", {
  res <- ss_power_contrast(
    c_weights     = c(1/3, 1/3, 1/3, -1),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    n_per_group   = 20
  )
  ncp <- res$value[res$term == "noncentral_t_parm"]
  N   <- res$value[res$term == "total_N"]
  f   <- res$value[res$term == "effect_size_f"]
  expect_equal(f, ncp / sqrt(N), tolerance = 1e-12)
})

test_that("n_per_group as a vector accepts unequal n and reports NA per-group", {
  res <- ss_power_contrast(
    c_weights     = c(0.5, 0.5, -0.5, -0.5),
    mu            = c(90, 92, 88, 81),
    sigma_squared = 144,
    n_per_group   = c(15, 25, 25, 15)
  )
  expect_true(is.na(res$value[res$term == "specified_n_per_group"]))
  expect_equal(res$value[res$term == "total_N"], 80)
})

test_that("a one-sided test gives strictly higher power than two-sided in the same direction", {
  two <- ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                           mu = c(90, 92, 88, 81), sigma_squared = 144,
                           n_per_group = 20, directional = FALSE)
  one <- ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                           mu = c(90, 92, 88, 81), sigma_squared = 144,
                           n_per_group = 20, directional = TRUE)
  expect_gt(one$value[one$term == "actual_power"],
            two$value[two$term == "actual_power"])
})

test_that("c_weights validation: sum to zero, positives sum to 1, negatives sum to -1", {
  bad_sum <- c(1, -0.5, -0.4)        # sums to 0.1, fails sum-to-zero
  expect_error(ss_power_contrast(c_weights = bad_sum, mu = c(1, 2, 3),
                                 sigma_squared = 1, n_per_group = 10),
               "sum to zero")

  bad_pos <- c(2, -1, -1)            # positives sum to 2, not 1
  expect_error(ss_power_contrast(c_weights = bad_pos, mu = c(1, 2, 3),
                                 sigma_squared = 1, n_per_group = 10),
               "positive entries")

  bad_neg <- c(0.5, 0.5, -0.5, -0.5) # OK
  res <- ss_power_contrast(c_weights = bad_neg, mu = c(1, 2, 3, 4),
                           sigma_squared = 1, n_per_group = 10)
  expect_s3_class(res, "data.frame")
})

test_that("requires sigma_squared and either mu or psi", {
  expect_error(ss_power_contrast(c_weights = c(1, -1)), "sigma_squared")
  expect_error(ss_power_contrast(c_weights = c(1, -1), sigma_squared = 1,
                                 n_per_group = 10),
               "psi.*mu|mu.*psi")
  expect_error(ss_power_contrast(c_weights = c(1, -1), sigma_squared = 1,
                                 mu = c(0, 1), psi = 1, n_per_group = 10),
               "either")
})

test_that("rejects n_per_group with wrong length", {
  expect_error(
    ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                      mu = c(90, 92, 88, 81), sigma_squared = 144,
                      n_per_group = c(20, 20, 20)),  # length 3, should be 4
    "n_per_group"
  )
})

test_that("psi == 0 raises a clear error when solving for n", {
  expect_error(
    ss_power_contrast(c_weights = c(1, -1), mu = c(5, 5),
                      sigma_squared = 1, desired_power = 0.80),
    "zero"
  )
})

test_that("returned data.frame has the expected term names and column types", {
  res_solve <- ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                                 mu = c(90, 92, 88, 81), sigma_squared = 144,
                                 desired_power = 0.80)
  expect_named(res_solve, c("term", "value"))
  expect_setequal(res_solve$term,
                  c("necessary_n_per_group", "total_N", "actual_power",
                    "noncentral_t_parm", "effect_size_f"))

  res_pwr <- ss_power_contrast(c_weights = c(1/3, 1/3, 1/3, -1),
                               mu = c(90, 92, 88, 81), sigma_squared = 144,
                               n_per_group = 20)
  expect_setequal(res_pwr$term,
                  c("specified_n_per_group", "total_N", "actual_power",
                    "noncentral_t_parm", "effect_size_f"))
})
