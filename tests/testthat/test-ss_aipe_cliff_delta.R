test_that("ss_aipe_cliff_delta() returns documented rows and width <= target", {
  res <- ss_aipe_cliff_delta(delta = 0.30, width = 0.20)
  expect_setequal(res$term,
                   c("n_1", "n_2", "necessary_N", "expected_width",
                     "delta", "ratio", "width_target", "conf_level"))
  expect_lte(res$value[res$term == "expected_width"], 0.20 + 1e-6)
})

test_that("ss_aipe_cliff_delta() unbalanced allocation needs larger total N", {
  bal   <- ss_aipe_cliff_delta(0.30, 0.20, ratio = 1)
  unbal <- ss_aipe_cliff_delta(0.30, 0.20, ratio = 3)
  expect_gt(unbal$value[unbal$term == "necessary_N"],
            bal$value[bal$term == "necessary_N"])
})

test_that("ss_aipe_cliff_delta() assurance inflates the sample size", {
  res_50 <- ss_aipe_cliff_delta(0.30, 0.20)
  res_80 <- ss_aipe_cliff_delta(0.30, 0.20, assurance = 0.80)
  expect_gt(res_80$value[res_80$term == "necessary_N"],
            res_50$value[res_50$term == "necessary_N"])
})

test_that("ss_aipe_cliff_delta() rejects bad inputs", {
  expect_error(ss_aipe_cliff_delta(1.2, 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_cliff_delta(0.3, 0),   "positive")
})
