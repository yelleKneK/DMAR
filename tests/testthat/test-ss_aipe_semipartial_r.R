test_that("ss_aipe_semipartial_r() returns documented rows and width <= target", {
  res <- ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15)
  expect_true(all(c("necessary_N", "expected_width", "r_sp", "J",
                    "width_target", "conf_level") %in% res$term))
  expect_lte(res$value[res$term == "expected_width"], 0.15 + 1e-6)
})

test_that("ss_aipe_semipartial_r() assurance inflates the sample size", {
  res_50 <- ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15)
  res_80 <- ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15,
                                   assurance = 0.80)
  expect_gt(res_80$value[1], res_50$value[1])
})

test_that("ss_aipe_semipartial_r() rejects bad inputs", {
  expect_error(ss_aipe_semipartial_r(r_sp = 1.5, J = 3, width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_semipartial_r(r_sp = 0.3, J = 0, width = 0.2), ">=")
})
