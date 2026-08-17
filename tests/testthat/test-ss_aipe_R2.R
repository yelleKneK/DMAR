## ss_aipe_R2() -- sample size planning for the squared multiple correlation.

test_that("ss_aipe_R2() returns a tidy data.frame with a positive necessary_N", {
  res <- suppressWarnings(ss_aipe_R2(population_R2 = 0.30, width = 0.20, p = 5))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_R2() requires more N for tighter target widths", {
  wide  <- suppressWarnings(ss_aipe_R2(0.30, width = 0.30, p = 5)$value[1])
  tight <- suppressWarnings(ss_aipe_R2(0.30, width = 0.10, p = 5)$value[1])
  expect_gt(tight, wide)
})

test_that("ss_aipe_R2() requires more N for more predictors at the same R^2", {
  few  <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 2)$value[1])
  many <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 10)$value[1])
  expect_gt(many, few)
})

test_that("ss_aipe_R2() accepts assurance > 0.5 and bumps N", {
  no_assurance <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 5)$value[1])
  with_85      <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 5, assurance = 0.85)$value[1])
  expect_gte(with_85, no_assurance)
})
