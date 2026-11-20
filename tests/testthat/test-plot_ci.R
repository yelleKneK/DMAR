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

test_that("plot_ci() writes the confidence level as it is, not rounded to a whole percent", {
  skip_if_not_installed("ggplot2")
  # round(100 * conf_level) labeled a 97.5% interval "98%" and a 99.9%
  # interval "100%", a level only an unbounded interval has.
  p <- plot_ci(estimate = 0.45, lower = 0.15, upper = 0.75,
               conf_level = 0.975)
  expect_identical(p$labels$title, "97.5% Confidence Interval")
  p <- plot_ci(estimate = 0.45, lower = 0.15, upper = 0.75,
               conf_level = 0.999)
  expect_identical(p$labels$title, "99.9% Confidence Interval")
})

test_that("plot_ci() titles a DMAR table with the level the table records", {
  skip_if_not_installed("ggplot2")
  # The title came from the conf_level argument, whose default of .95 titled
  # a 97.5% or a 90% table "95% Confidence Interval".
  p <- plot_ci(ci = ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.975))
  expect_identical(p$labels$title, "97.5% Confidence Interval")
  expect_null(p$labels$subtitle)
  ci90 <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.90)
  expect_identical(plot_ci(ci = ci90)$labels$title, "90% Confidence Interval")
  expect_identical(plot_ci(ci = ci90, conf_level = 0.90)$labels$title,
                   "90% Confidence Interval")
  # A level that contradicts the table's is refused: the limits drawn are
  # the table's, at its level.
  expect_error(plot_ci(ci = ci90, conf_level = 0.95),
               "not the level the 'ci' table records for its limits, 90%")
  # A wide table records its level too.
  pyg <- pygmalion
  pyg$grade <- factor(pyg$grade)
  om <- ci_omega_squared(stats::aov(iq_8 ~ treatment + grade, data = pyg),
                         conf_level = 0.90)
  expect_identical(plot_ci(om)$labels$title, "90% Confidence Interval")
  # A frame that records no level keeps the argument's level.
  hand <- data.frame(term = c("ci_lower", "d", "ci_upper"),
                     value = c(0.1, 0.4, 0.7))
  expect_identical(plot_ci(ci = hand)$labels$title, "95% Confidence Interval")
  expect_identical(plot_ci(ci = hand, conf_level = 0.99)$labels$title,
                   "99% Confidence Interval")
  # Unequal tail areas are named in a subtitle, and a one-sided interval
  # says which limit is the confidence bound.
  p <- plot_ci(ci = ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0.01,
                           alpha_upper = 0.04))
  expect_identical(p$labels$title, "95% Confidence Interval")
  expect_identical(p$labels$subtitle,
                   "Unequal tail areas: 0.01 below, 0.04 above")
  p <- plot_ci(ci = ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0,
                           alpha_upper = 0.05))
  expect_identical(p$labels$subtitle,
                   "One-sided: the upper limit is the 95% upper confidence bound")
  p <- plot_ci(ci = ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0.10,
                           alpha_upper = 0))
  expect_identical(p$labels$title, "90% Confidence Interval")
  expect_identical(p$labels$subtitle,
                   "One-sided: the lower limit is the 90% lower confidence bound")
  # Rows that record different levels cannot share one title.
  mixed <- data.frame(effect = c("a", "b"), estimate = c(0.2, 0.3),
                      ci_lower = c(0.1, 0.1), ci_upper = c(0.3, 0.5),
                      conf_level = c(0.95, 0.90))
  expect_error(plot_ci(mixed), "record different confidence levels")
})

test_that("plot_ci() draws each estimate of an ANOVA table that has limits", {
  skip_if_not_installed("ggplot2")
  pyg <- pygmalion
  pyg$grade <- factor(pyg$grade)
  # A factorial fit: the totals have no limits, so the partial values
  # alone are drawn, one row per effect.
  fac <- ci_omega_squared(stats::aov(iq_8 ~ treatment + grade, data = pyg))
  parsed <- DMAR:::.parse_dmar_ci(fac)
  expect_identical(parsed$names, c("treatment", "grade"))
  expect_identical(parsed$lower, fac$omega_squared_partial_ci_lower)
  expect_s3_class(plot_ci(fac), "ggplot")
  # A one-way fit: both estimands have limits, and each row names its
  # estimate.
  one <- ci_eta_squared(stats::aov(iq_gain ~ grade, data = pyg))
  parsed <- DMAR:::.parse_dmar_ci(one)
  expect_identical(parsed$names, c("grade: eta_squared_total",
                                   "grade: eta_squared_partial"))
  expect_identical(parsed$n, c(310L, 310L))
})

test_that("plot_ci() draws the estimates of ci_eta_squared_generalized(), not the midpoints", {
  # Through the first builds of 1.1.0 the parser looked for an
  # omega_squared column and, finding none in this table, plotted the
  # midpoint of each interval.
  g <- suppressWarnings(ci_eta_squared_generalized(
    aov(iq_8 ~ treatment * factor(grade), data = pygmalion),
    method = "parametric"))
  parsed <- DMAR:::.parse_dmar_ci(g)
  expect_identical(parsed$estimate, g$eta_squared_generalized)
  expect_identical(parsed$lower, g$ci_lower)
})
