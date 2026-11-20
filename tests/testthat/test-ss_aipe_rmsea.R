## ss_aipe_rmsea() -- sample size planning for the RMSEA confidence interval.
##
## suppressMessages throughout: during the iterative search the intermediate
## sample sizes can put the lower confidence limit of the noncentrality
## parameter at its bound, and ci_rmsea() notes the clamp with a message
## that is not what these tests check.

test_that("ss_aipe_rmsea() returns a tidy data.frame with a positive necessary_N", {
  res <- suppressMessages(ss_aipe_rmsea(true_rmsea = 0.05, df = 20, width = 0.05))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_rmsea() requires more N for tighter widths", {
  skip_on_cran()
  wide  <- suppressMessages(
    ss_aipe_rmsea(true_rmsea = 0.05, df = 20, width = 0.08))$value[1]
  tight <- suppressMessages(
    ss_aipe_rmsea(true_rmsea = 0.05, df = 20, width = 0.02))$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_rmsea() requires fewer N at higher df (more model information)", {
  skip_on_cran()
  df_low  <- suppressMessages(
    ss_aipe_rmsea(true_rmsea = 0.05, df = 10, width = 0.05))$value[1]
  df_high <- suppressMessages(
    ss_aipe_rmsea(true_rmsea = 0.05, df = 60, width = 0.05))$value[1]
  expect_gt(df_low, df_high)
})

test_that("the default plans a 90 percent interval, and 0.95 restores DMAR 1.0.x", {
  skip_on_cran()
  # Since DMAR 1.1.0 the planned interval is the conventional 90 percent
  # one, the level ci_rmsea() forms by default. For an RMSEA of .045 on 5
  # df and a width of .060 the plan is 802 at 90 percent, against the 1068
  # DMAR 1.0.x gave at its default of 95 percent.
  expect_equal(eval(formals(ss_aipe_rmsea)$conf_level), 0.90)
  at_90 <- suppressMessages(ss_aipe_rmsea(true_rmsea = .045, df = 5, width = .060))
  at_95 <- suppressMessages(ss_aipe_rmsea(true_rmsea = .045, df = 5, width = .060,
                                          conf_level = 0.95))
  expect_equal(at_90$value, 802)
  expect_equal(at_95$value, 1068)
  expect_equal(attr(at_90, "conf_level"), 0.90)
  expect_true("Confidence level: 90%" %in%
                utils::capture.output(print(at_90)))
  # The plan meets the width at its level: the 90 percent interval at the
  # planned N is no wider than .060, and one case fewer is.
  w <- function(N) {
    ci <- suppressMessages(ci_rmsea(.045, 5, N, conf_level = 0.90))
    ci$value[3] - ci$value[1]
  }
  expect_lte(w(802), .060)
  expect_gt(w(801), .060)
})

test_that("a plan at the default level gives no warning", {
  skip_on_cran()  # two more plans; the default itself is checked above
  expect_no_warning(suppressMessages(
    ss_aipe_rmsea(true_rmsea = .045, df = 5, width = .060, conf_level = 0.90)))
  expect_no_warning(suppressMessages(
    ss_aipe_rmsea(true_rmsea = .045, df = 5, width = .060)))
})
