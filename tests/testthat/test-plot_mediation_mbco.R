# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for plot_mediation_mbco(): the conditional-effect display of a
# moderated mediation_mbco() analysis. The central check is agreement
# between the drawn curve and the table's probe rows: both must
# evaluate the same stored moderator polynomial.

test_that("plot_mediation_mbco draws the stored conditional effects", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("ggplot2")

  set.seed(113)
  n <- 250
  x <- rnorm(n)
  w <- rnorm(n)
  m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
  y <- 0.5 * m + 0.2 * x + 0.1 * w + rnorm(n)
  d <- data.frame(x = x, w = w, m = m, y = y)
  res <- mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w", data = d,
                        x = "x", y = "y", moderator = "w",
                        ci_method = "wald")

  p <- plot_mediation_mbco(res, seed = 113)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("effect_label", "w_value", "estimate",
                    "band_lower", "band_upper") %in% names(p$data)))
  # Both moderated effects are drawn; the unmoderated direct effect is
  # not.
  expect_identical(nlevels(p$data$effect_label), 2L)
  expect_false(any(grepl("^x -> y$", levels(p$data$effect_label))))
  # The band brackets the curve.
  expect_true(all(p$data$band_lower <= p$data$estimate + 1e-8))
  expect_true(all(p$data$band_upper >= p$data$estimate - 1e-8))

  # The curve is the same quantity as the table's probe rows: evaluate
  # the plot at the probed moderator value and compare.
  v_mean <- unname(attr(res, "moderation")$values["mean"])
  p2 <- plot_mediation_mbco(res, effects = "indirect_via_m",
                            from = v_mean, to = v_mean + 1,
                            n_grid = 2, B = 200, seed = 113)
  expect_equal(p2$data$estimate[1L],
               res$estimate[res$term == "indirect_via_m_at_mean"],
               tolerance = 1e-8)

  # Effect selection is validated with the available names.
  expect_error(plot_mediation_mbco(res, effects = "nope"),
               "Unknown effect")
  expect_error(plot_mediation_mbco(res, from = 2, to = -2),
               "smaller than")

  # The caption names the band's level as it is: round(100 * conf_level)
  # captioned a 97.5% band "98%".
  p3 <- plot_mediation_mbco(res, effects = "indirect_via_m",
                            conf_level = 0.975, n_grid = 2, B = 200,
                            seed = 113)
  expect_identical(p3$labels$caption,
                   "97.5% pointwise Monte Carlo confidence band")

  # An unmoderated result has nothing to draw and says so.
  res_plain <- mediation_mbco("m ~ x \n y ~ m + x",
                              data = d[, c("x", "m", "y")],
                              x = "x", y = "y", ci_method = "wald")
  expect_error(plot_mediation_mbco(res_plain), "moderator")
})

test_that("the band of a badly scaled moderator is that of the variables as given (review)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("ggplot2")
  # With the moderator times 10,000 mediation_mbco() fits a rescaled
  # metric (its "rescale" attribute). The conditional effects are the same
  # functions of the moderator in its own units, so the curve at the
  # matching points is the same, and the band agrees within Monte Carlo
  # error (the draws are taken in each fit's metric). Drawn from the
  # covariance matrix of the original metric, whose variances spanned 25
  # orders of magnitude, the band of the first draft was about 3000 times
  # too wide.
  set.seed(3812)
  n <- 300
  w <- stats::rnorm(n)
  x <- stats::rnorm(n)
  m <- .4 * x + .2 * x * w + stats::rnorm(n)
  y <- .3 * m + .1 * x + stats::rnorm(n)
  # lavaan's note that it rescaled the variables itself is a message.
  fit <- function(k) {
    suppressMessages(mediation_mbco(
      "m ~ x + w + x:w\ny ~ m + x", data = data.frame(x, w = w * k, m, y),
      x = "x", y = "y", moderator = "w", probe_values = c(-k, 0, k),
      ci_method = "wald"))
  }
  r0 <- fit(1)
  r1 <- fit(1e4)
  expect_false(is.null(attr(r1, "rescale")))
  b0 <- plot_mediation_mbco(r0, B = 2000, seed = 11, n_grid = 3)$data
  b1 <- plot_mediation_mbco(r1, B = 2000, seed = 11, n_grid = 3)$data
  expect_equal(b1$w_value, b0$w_value * 1e4, tolerance = 1e-10)
  expect_equal(b1$estimate, b0$estimate, tolerance = 1e-8)
  width0 <- b0$band_upper - b0$band_lower
  width1 <- b1$band_upper - b1$band_lower
  expect_true(all(abs(width1 / width0 - 1) < 0.2))
})
