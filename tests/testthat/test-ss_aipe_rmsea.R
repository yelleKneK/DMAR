## ss_aipe_rmsea() -- sample size planning for the RMSEA confidence interval.
##
## suppressMessages throughout: during the iterative search the intermediate
## sample sizes can put the lower confidence limit of the noncentrality
## parameter at its bound, and ci_rmsea() notes the clamp with a message
## that is not what these tests check.

test_that("ss_aipe_rmsea() returns a tidy data.frame with a positive necessary_N", {
  res <- suppressMessages(ss_aipe_rmsea(RMSEA = 0.05, df = 20, width = 0.05))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_rmsea() requires more N for tighter widths", {
  wide  <- suppressMessages(
    ss_aipe_rmsea(RMSEA = 0.05, df = 20, width = 0.08))$value[1]
  tight <- suppressMessages(
    ss_aipe_rmsea(RMSEA = 0.05, df = 20, width = 0.02))$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_rmsea() requires fewer N at higher df (more model information)", {
  df_low  <- suppressMessages(
    ss_aipe_rmsea(RMSEA = 0.05, df = 10, width = 0.05))$value[1]
  df_high <- suppressMessages(
    ss_aipe_rmsea(RMSEA = 0.05, df = 60, width = 0.05))$value[1]
  expect_gt(df_low, df_high)
})
