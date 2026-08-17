## ss_aipe_rc() / ss_aipe_src() -- sample size for unstandardized / standardized
## regression coefficient AIPE planning.

test_that("ss_aipe_rc() returns a tidy data.frame with necessary_N", {
  res <- ss_aipe_rc(rho2_Y_X = 0.5, Rho2_j_X_without_j = 0.3,
                    p = 3, b_j = 0.4, width = 0.20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_rc() requires more N for tighter widths", {
  wide  <- ss_aipe_rc(0.5, 0.3, p = 3, b_j = 0.4, width = 0.30)$value[1]
  tight <- ss_aipe_rc(0.5, 0.3, p = 3, b_j = 0.4, width = 0.10)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_src() returns a tidy data.frame with necessary_N", {
  res <- ss_aipe_src(rho2_Y_X = 0.5, Rho2_j_X_without_j = 0.3,
                     p = 3, beta_j = 0.4, width = 0.20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_src() requires more N for tighter widths", {
  wide  <- ss_aipe_src(0.5, 0.3, p = 3, beta_j = 0.4, width = 0.30)$value[1]
  tight <- ss_aipe_src(0.5, 0.3, p = 3, beta_j = 0.4, width = 0.10)$value[1]
  expect_gt(tight, wide)
})
