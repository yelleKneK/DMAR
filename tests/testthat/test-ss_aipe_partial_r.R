test_that("ss_aipe_partial_r() returns documented rows", {
  res <- ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20)
  expect_true(all(c("necessary_N", "expected_width", "rho", "J",
                    "width_target", "conf_level") %in% res$term))
  n_val <- res$value[res$term == "necessary_N"]
  expect_gt(n_val, 3)   # sanity: more than k + 1 = 3
})

test_that("ss_aipe_partial_r() yields expected width <= target", {
  res <- ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20)
  expect_lte(res$value[res$term == "expected_width"], 0.20 + 1e-6)
})

test_that("ss_aipe_partial_r() Fisher's Z route gives a similar sample size", {
  res_raw <- ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20,
                                fisher_z = FALSE)
  res_fz  <- ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20,
                                fisher_z = TRUE)
  # Two routes should be in the same ballpark (within a factor of 2):
  ratio <- res_fz$value[1] / res_raw$value[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2)
})

test_that("ss_aipe_partial_r() assurance inflates the sample size", {
  res_50  <- ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20)
  res_80  <- ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20,
                                assurance = 0.80)
  expect_gt(res_80$value[1], res_50$value[1])
})

test_that("ss_aipe_partial_r() rejects bad inputs", {
  expect_error(ss_aipe_partial_r(rho = 1.5, J = 2, width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_partial_r(rho = 0.3, J = 0, width = 0.2), ">=")
  expect_error(ss_aipe_partial_r(rho = 0.3, J = 2, width = 0),   "positive")
})
