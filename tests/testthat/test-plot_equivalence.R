test_that("plot_equivalence() returns a ggplot from direct vectors", {
  skip_if_not_installed("ggplot2")
  p <- plot_equivalence(estimate = c(-1.0, 3.5, -3.6, -1.5, -8.0),
                        lower    = c(-3.2, -1.4, -4.8, -6.6, -10.5),
                        upper    = c( 1.2,  8.4, -2.4,  3.6,  -5.5),
                        names    = c("A", "B", "C", "D", "E"),
                        delta_upper = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_equivalence() accepts equivalence_c() results and reads their
           bounds", {
  skip_if_not_installed("ggplot2")
  res <- list(
    focal  = equivalence_c(psi_hat = -5.28, se = 2.49, df_error = 399,
                    delta_upper = 5),
    within = equivalence_c(psi_hat = -0.53, se = 2.66, df_error = 399,
                    delta_upper = 5)
  )
  p <- plot_equivalence(res)
  expect_s3_class(p, "ggplot")
  # A single result works too.
  p1 <- plot_equivalence(res$focal)
  expect_s3_class(p1, "ggplot")
})

test_that("plot_equivalence() verdicts match equivalence_c()", {
  skip_if_not_installed("ggplot2")
  res <- equivalence_c(psi_hat = -0.53, se = 2.66, df_error = 399, delta_upper = 5)
  p <- plot_equivalence(res)
  drawn <- as.character(p$data$verdict)
  expect_identical(drawn, attr(res, "verdict"))
})

test_that("plot_equivalence() maps all five verdicts and honors the palette", {
  skip_if_not_installed("ggplot2")
  # One interval per verdict against symmetric bounds of 5.
  args <- list(estimate = c(-1.0, 3.5, 7.0, -1.5, -8.0),
               lower    = c(-3.2, -1.4, 5.5, -6.6, -10.5),
               upper    = c( 1.2,  8.4, 8.5,  3.6,  -5.5),
               names    = c("A", "B", "C", "D", "E"),
               delta_upper = 5)
  p <- do.call(plot_equivalence, args)
  expect_setequal(as.character(p$data$verdict),
                  c("Equivalent", "Superior", "Noninferior only",
                    "Inconclusive", "Inferior"))
  # The default palette is the neutral colorblind-safe okabe_ito set: every
  # point color is drawn from it, and it is not the Notre Dame palette.
  b   <- ggplot2::ggplot_build(p)
  col <- b$data[[length(b$data)]]$colour           # the geom_pointrange layer
  expect_true(all(col %in% .dmar_palette(5, palette = "okabe_ito")))
  # A different palette produces different point colors.
  p2  <- do.call(plot_equivalence, c(args, palette = "tableau"))
  b2  <- ggplot2::ggplot_build(p2)
  col2 <- b2$data[[length(b2$data)]]$colour
  expect_false(isTRUE(all.equal(sort(col), sort(col2))))
})

test_that("plot_equivalence() refuses mixed bounds and missing inputs", {
  skip_if_not_installed("ggplot2")
  res <- list(
    a = equivalence_c(psi_hat = 0, se = 1, df_error = 50, delta_upper = 5),
    b = equivalence_c(psi_hat = 0, se = 1, df_error = 50, delta_upper = 3)
  )
  expect_error(plot_equivalence(res), "different equivalence bounds")
  expect_error(plot_equivalence(estimate = 1, lower = 0, upper = 2),
               "delta_upper")
  expect_error(plot_equivalence(estimate = c(1, 2), lower = 0, upper = 2,
                                delta_upper = 5),
               "same length")
})

test_that("plot_equivalence() honors an explicit bound override (MEDIUM-06)", {
  skip_if_not_installed("ggplot2")
  t1 <- equivalence_c(psi_hat = 1, se = 1, df_error = 30, delta_lower = 5, delta_upper = 5)
  t2 <- equivalence_c(psi_hat = 1, se = 1, df_error = 30, delta_lower = 3, delta_upper = 8)
  # Differing stored bounds, but an explicit common region is supplied: proceed.
  p <- plot_equivalence(list(A = t1, B = t2), delta_lower = 5, delta_upper = 5)
  expect_s3_class(p, "ggplot")
  # Without the override, differing stored bounds still error.
  expect_error(plot_equivalence(list(A = t1, B = t2)), "different equivalence bounds")
})
