test_that("plot_forest() returns a ggplot with study and pooled rows", {
  skip_if_not_installed("ggplot2")
  d  <- teacher_expectancy$d
  ne <- teacher_expectancy$n_experimental
  nc <- teacher_expectancy$n_control
  vi <- (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc))
  p <- plot_forest(d, vi, labels = teacher_expectancy$author)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 19 + 1)               # studies + pooled row
  expect_true("Random effects" %in% p$data$label)
  # Fixed effect variant labels its pooled row accordingly.
  pf <- plot_forest(d, vi, method = "fe")
  expect_true("Common effect" %in% pf$data$label)
})

test_that("plot_forest() validates labels", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_forest(c(.2, .3), c(.01, .02), labels = "one"),
               "one label per study")
})

test_that("the @examples data are genuinely heterogeneous", {
  # The example promises a prediction interval visibly wider than the
  # confidence interval; hold the example's simulated studies to that.
  set.seed(113)
  k <- 12
  n <- sample(20:100, k)
  theta <- rnorm(k, mean = 0.4, sd = 0.35)
  d <- rnorm(k, mean = theta, sd = sqrt(2 / n))
  v <- 2 / n + d^2 / (4 * n)
  fit <- meta_es(d, v)
  tau2 <- fit$value[fit$term == "tau2"]
  ci_width <- fit$value[fit$term == "upper_limit"] -
    fit$value[fit$term == "lower_limit"]
  pi_width <- fit$value[fit$term == "prediction_upper"] -
    fit$value[fit$term == "prediction_lower"]
  expect_gt(tau2, 0)
  expect_gt(pi_width / ci_width, 1.3)
})

test_that("plot_forest() draws intervals without deprecated ggplot2 geoms", {
  skip_if_not_installed("ggplot2")
  set.seed(113)
  k <- 12
  n <- sample(20:100, k)
  theta <- rnorm(k, mean = 0.4, sd = 0.35)
  d <- rnorm(k, mean = theta, sd = sqrt(2 / n))
  v <- 2 / n + d^2 / (4 * n)
  expect_silent({
    p <- plot_forest(d, v)
    built <- ggplot2::ggplot_build(p)
  })
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_false("GeomErrorbarh" %in% geoms)
  expect_true("GeomErrorbar" %in% geoms)
})
