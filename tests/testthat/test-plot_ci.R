test_that("plot_ci() returns a ggplot from explicit values", {
  skip_if_not_installed("ggplot2")
  p <- plot_ci(estimate = 0.45, lower = 0.15, upper = 0.75,
               names = "d", n = 60)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ci() parses ci_smd() output", {
  skip_if_not_installed("ggplot2")
  ci_result <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  p <- plot_ci(ci_result, n = 100, reference_line = 0)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ci() parses ci_omega_squared() output", {
  skip_if_not_installed("ggplot2")
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(omega_result <- ci_omega_squared(fit),
                 "below the alpha_lower critical value")
  p <- plot_ci(omega_result, reference_line = 0)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ci() handles multiple explicit effects", {
  skip_if_not_installed("ggplot2")
  p <- plot_ci(
    estimate = c(0.45, 0.20, -0.10),
    lower    = c(0.10, -0.15, -0.45),
    upper    = c(0.80,  0.55,  0.25),
    names    = c("Outcome A", "Outcome B", "Outcome C"),
    n        = c(100, 80, 120),
    reference_line = 0
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_ci() errors without lower/upper and no ci data.frame", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_ci(estimate = 0.5), "Confidence limits are required")
})

test_that("plot_ci() works without n (show_n has no effect)", {
  skip_if_not_installed("ggplot2")
  p <- plot_ci(estimate = 0.30, lower = 0.05, upper = 0.55, show_n = TRUE)
  expect_s3_class(p, "ggplot")
})
