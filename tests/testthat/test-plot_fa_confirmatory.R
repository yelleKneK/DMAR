plot_cfa_k_result <- function(...) {
  data(holzinger_swineford, package = "DMAR", envir = environment())
  fa_confirmatory(holzinger_swineford,
        list(verbal = c("t6_paragraph_comprehension",
                        "t7_sentence", "t9_word_meaning"),
             deduction = c("t20_deduction", "t22_problem_reasoning",
                        "t23_series_completion")),
        ...)
}

test_that("plot_fa_confirmatory() returns a ggplot for each parameter type", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("ggplot2")

  res <- plot_cfa_k_result()
  expect_s3_class(plot_fa_confirmatory(res), "ggplot")
  expect_s3_class(plot_fa_confirmatory(res, what = "errors"), "ggplot")

  res_means <- plot_cfa_k_result(meanstructure = TRUE)
  expect_s3_class(plot_fa_confirmatory(res_means, what = "intercepts"), "ggplot")
})

test_that("plot_fa_confirmatory() displays equated estimates", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("ggplot2")

  res <- plot_cfa_k_result(equal_loading = TRUE)
  p <- plot_fa_confirmatory(res)
  expect_s3_class(p, "ggplot")
  # Every item of a factor sits at the single equated loading (the values
  # can differ in the last floating point bits across lavaan's rows).
  expect_lt(diff(range(p$data$estimate[p$data$factor == "verbal"])), 1e-8)
})

test_that("plot_fa_confirmatory() validates its input", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("ggplot2")

  # A listwise fit has no mean structure; the FIML default has one.
  res <- plot_cfa_k_result(missing = "listwise")
  expect_error(plot_fa_confirmatory(res, what = "intercepts"), "mean structure")
  expect_error(plot_fa_confirmatory(data.frame(term = "lambda_f_1")),
               "fa_confirmatory\\(\\) result")

  meas <- plot_cfa_k_result(output = "measurement")
  expect_error(plot_fa_confirmatory(meas), "fa_confirmatory\\(\\) result")
})

test_that("plot_cfa_k() is plot_fa_confirmatory() under its 1.0.0 name", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("ggplot2")
  expect_identical(formals(plot_cfa_k), formals(plot_fa_confirmatory))
  expect_s3_class(plot_cfa_k(plot_cfa_k_result()), "ggplot")
})
