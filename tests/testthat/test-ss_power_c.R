
test_that("ss_power_c() power is exact past a noncentrality of 37.62", {
  # Two per group, a contrast 60 standard deviations from 0 (ncp 60, 2 df),
  # alpha .001: pt() gave .9794981672 in DMAR through 3338679.
  r <- ss_power_c(true_contrast = 60, c_weights = c(1, -1), true_sd = 1,
                  specified_n_per_group = 2, alpha_level = .001)
  crit <- qt(1 - .001 / 2, 2)
  expect_equal(r$value[r$term == "actual_power"], 0.9726544235, tolerance = 1e-9)
  expect_equal(r$value[r$term == "actual_power"],
               .oracle_nct(crit, 2, 60, lower_tail = FALSE) + .oracle_nct(-crit, 2, 60),
               tolerance = 1e-12)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_power_c returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_c(true_contrast = 0.6, c_weights = c(1, -1), true_sd = 1, specified_n_per_group = n)
    r$value[r$term == "actual_power"] }
  target <- pw(2) - 1e-4                       # a target the minimum n = 2 already attains
  r <- ss_power_c(true_contrast = 0.6, c_weights = c(1, -1), true_sd = 1, desired_power = target)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})
