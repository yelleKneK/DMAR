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

  # An unmoderated result has nothing to draw and says so.
  res_plain <- mediation_mbco("m ~ x \n y ~ m + x",
                              data = d[, c("x", "m", "y")],
                              x = "x", y = "y", ci_method = "wald")
  expect_error(plot_mediation_mbco(res_plain), "moderator")
})
