test_that("plot_R2() returns a ggplot with CI and n", {
  skip_if_not_installed("ggplot2")
  p <- plot_R2(R2 = 0.25, N = 100, p = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_R2() works without CI or n", {
  skip_if_not_installed("ggplot2")
  p <- plot_R2(R2 = 0.10, show_ci = FALSE, show_n = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_R2() handles small R2 (label placement)", {
  skip_if_not_installed("ggplot2")
  p <- plot_R2(R2 = 0.05, N = 200, p = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_R2() handles large R2", {
  skip_if_not_installed("ggplot2")
  p <- plot_R2(R2 = 0.85, N = 50, p = 8)
  expect_s3_class(p, "ggplot")
})

test_that("plot_R2() errors on invalid R2", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_R2(R2 = -0.1), "between 0 and 1")
  expect_error(plot_R2(R2 = 1.5),  "between 0 and 1")
})

test_that("plot_R2() warns when CI requested without N or K", {
  skip_if_not_installed("ggplot2")
  expect_warning(plot_R2(R2 = 0.25, show_ci = TRUE), "suppressing CI")
})
