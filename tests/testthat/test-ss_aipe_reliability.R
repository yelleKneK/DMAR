## ss_aipe_reliability() -- sample size for the alpha / omega reliability CI.

test_that("ss_aipe_reliability(parallel/Normal Theory) returns a tidy necessary_N", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                        width = 0.10, i = 4, lambda = 0.7, psi_square = 0.51,
                        initial_iter = 50, final_iter = 200)
  ))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_reliability(Congeneric / Factor Analytic) plans a sample size", {
  # Regression test: the Factor Analytic path read the legacy cfa_1() list
  # layout ($factor_loadings / $parameter_cov), which no longer exists, so the
  # path errored. It now routes through the maintained omega internals.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  res <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = 0.10, i = 5,
                        lambda = c(.4, .4, .3, .3, .5),
                        psi_square = c(.2, .4, .3, .3, .2),
                        conf_level = 0.95)
  ))
  expect_s3_class(res, "data.frame")
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_reliability() requires more N for a tighter width target", {
  skip_on_cran()
  set.seed(113)
  wide <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                        width = 0.20, i = 4, lambda = 0.7, psi_square = 0.51,
                        initial_iter = 50, final_iter = 200)
  ))$value[1]
  set.seed(113)
  tight <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                        width = 0.05, i = 4, lambda = 0.7, psi_square = 0.51,
                        initial_iter = 50, final_iter = 200)
  ))$value[1]
  expect_gt(tight, wide)
})
