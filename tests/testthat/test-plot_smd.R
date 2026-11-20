test_that("plot_smd() returns a ggplot object from summary values", {
  skip_if_not_installed("ggplot2")
  p <- plot_smd(smd = 0.50, n_1 = 50, n_2 = 50)
  expect_s3_class(p, "ggplot")
})

test_that("plot_smd() returns a ggplot from raw data", {
  skip_if_not_installed("ggplot2")
  set.seed(113)
  g1 <- rnorm(30, mean = 0.5)
  g2 <- rnorm(30, mean = 0.0)
  p <- plot_smd(group_1 = g1, group_2 = g2)
  expect_s3_class(p, "ggplot")
})

test_that("plot_smd() works without CI and without n", {
  skip_if_not_installed("ggplot2")
  p <- plot_smd(smd = 0.80, show_ci = FALSE, show_n = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_smd() warns when smd and raw data both provided", {
  skip_if_not_installed("ggplot2")
  g1 <- rnorm(20); g2 <- rnorm(20)
  expect_warning(plot_smd(smd = 0.5, group_1 = g1, group_2 = g2),
                 "computing from the data")
})

test_that("plot_smd() errors without smd or raw data", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_smd(), "Provide either raw data")
})

test_that("plot_smd() handles negative d", {
  skip_if_not_installed("ggplot2")
  p <- plot_smd(smd = -0.60, n_1 = 40, n_2 = 40)
  expect_s3_class(p, "ggplot")
})

test_that("plot_smd() accepts custom group labels and title", {
  skip_if_not_installed("ggplot2")
  p <- plot_smd(smd = 0.50, n_1 = 30, n_2 = 30,
                group_labels = c("Treatment", "Control"),
                title = "My Custom Title")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "My Custom Title")
})

test_that("plot_smd() writes the confidence level as it is, not rounded to a whole percent", {
  skip_if_not_installed("ggplot2")
  # round(100 * conf_level) labeled a 97.5% interval "98% CI" and a 99.9%
  # interval "100% CI".
  labels_of <- function(p) {
    unlist(lapply(ggplot2::ggplot_build(p)$data,
                  function(d) as.character(d$label)))
  }
  lab <- labels_of(plot_smd(smd = 0.5, n_1 = 50, n_2 = 50,
                            conf_level = 0.975))
  expect_true(any(startsWith(lab, "97.5% CI [")))
  expect_false(any(startsWith(lab, "98% CI")))
  lab <- labels_of(plot_smd(smd = 0.5, n_1 = 50, n_2 = 50,
                            conf_level = 0.999))
  expect_true(any(startsWith(lab, "99.9% CI [")))
  expect_false(any(startsWith(lab, "100% CI")))
})
