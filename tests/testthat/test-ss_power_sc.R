
test_that("ss_power_sc returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_sc(true_std_contrast = 0.6, c_weights = c(1, -1), specified_n_per_group = n)
    r$value[r$term == "actual_power"] }
  target <- pw(2) - 1e-4
  r <- ss_power_sc(true_std_contrast = 0.6, c_weights = c(1, -1), desired_power = target)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})
