# Oracles used here, none of which is DMAR's own output:
#   (a) the quadratic roots computed by hand from coef() and vcov() of the
#       fitted lm, with base R's qf() / qt() for the critical value;
#   (b) the definitional check that at a boundary the group difference sits
#       exactly on the critical value, |D(x)| / SE[D(x)] = t_crit, evaluated
#       through a linear combination of the coefficients rather than through
#       the quadratic;
#   (c) the published boundaries for the Pygmalion teacher-expectancy
#       example in Maxwell, Delaney, and Kelley (2027, Chapter 9 extension),
#       X = 103.9 and X = 128.7.

# Roots of D(x)^2 = t_crit^2 * Var[D(x)], from the fitted model's
# coefficients and covariance matrix. Two-group case, treatment contrasts,
# so the group difference and the slope difference are single coefficients,
# named by 'i' and 'j'.
hand_roots <- function(fit, conf_level = 0.95, method = "simultaneous",
                       i = "gtreated", j = "x:gtreated") {
  b  <- coef(fit)
  V  <- vcov(fit)
  df <- fit$df.residual
  d0 <- b[[i]]; d1 <- b[[j]]
  v00 <- V[i, i]; v01 <- V[i, j]; v11 <- V[j, j]
  tcrit <- if (method == "simultaneous") {
    sqrt(2 * qf(conf_level, 2, df))
  } else {
    qt(1 - (1 - conf_level) / 2, df)
  }
  A <- d1^2 - tcrit^2 * v11
  B <- 2 * (d0 * d1 - tcrit^2 * v01)
  C <- d0^2 - tcrit^2 * v00
  disc <- B^2 - 4 * A * C
  if (disc <= 0) return(c(NA_real_, NA_real_))
  sort(c((-B - sqrt(disc)) / (2 * A), (-B + sqrt(disc)) / (2 * A)))
}

bound <- function(res, which_bound, pair = NULL) {
  keep <- res$term == which_bound
  if (!is.null(pair)) keep <- keep & res$pair == pair
  res$value[keep]
}

sim_two_group <- function(n = 120) {
  set.seed(113)
  g <- factor(rep(c("control", "treated"), each = n / 2))
  x <- rnorm(n, 50, 10)
  y <- 10 + 0.5 * x + (g == "treated") * (0.35 * x - 14) + rnorm(n, 0, 6)
  data.frame(y, x, g)
}


test_that("regions_of_significance() reproduces the hand-computed roots", {
  d <- sim_two_group()
  fit <- lm(y ~ x * g, data = d)

  for (m in c("simultaneous", "pointwise")) {
    res <- regions_of_significance(fit, method = m)
    target <- hand_roots(fit, method = m)
    expect_equal(bound(res, "lower_bound"), target[1], tolerance = 1e-10)
    expect_equal(bound(res, "upper_bound"), target[2], tolerance = 1e-10)
  }
})


test_that("the critical value is Potthoff's simultaneous value by default", {
  d <- sim_two_group()
  fit <- lm(y ~ x * g, data = d)
  df <- fit$df.residual

  res_sim <- regions_of_significance(fit)
  res_pt  <- regions_of_significance(fit, method = "pointwise")
  expect_equal(bound(res_sim, "critical_value"),
               sqrt(2 * qf(0.95, 2, df)), tolerance = 1e-12)
  expect_equal(bound(res_pt, "critical_value"),
               qt(0.975, df), tolerance = 1e-12)
  # Scanning the whole covariate costs something, always.
  expect_gt(bound(res_sim, "critical_value"), bound(res_pt, "critical_value"))

  res_90 <- regions_of_significance(fit, conf_level = 0.90)
  expect_equal(bound(res_90, "critical_value"),
               sqrt(2 * qf(0.90, 2, df)), tolerance = 1e-12)
  expect_equal(bound(res_90, "conf_level"), 0.90)
  expect_equal(bound(res_sim, "df_error"), df)
})


test_that("at a boundary the difference sits exactly on the critical value", {
  # An independent route to the same claim: form D(x) and its standard error
  # as a linear combination of the coefficients and check the ratio, rather
  # than solving the quadratic.
  d <- sim_two_group()
  fit <- lm(y ~ x * g, data = d)
  b <- coef(fit); V <- vcov(fit)

  t_ratio <- function(x) {
    cvec <- c(0, 0, 1, x)
    abs(sum(cvec * b)) / sqrt(drop(t(cvec) %*% V %*% cvec))
  }

  for (m in c("simultaneous", "pointwise")) {
    res <- regions_of_significance(fit, method = m)
    tcrit <- bound(res, "critical_value")
    for (nm in c("lower_bound", "upper_bound")) {
      expect_equal(t_ratio(bound(res, nm)), tcrit, tolerance = 1e-9)
    }
    # Just inside the nonsignificant side, the ratio drops below the
    # critical value; just outside, it exceeds it. Here the leading
    # coefficient is positive, so significance lies outside the bounds.
    expect_equal(bound(res, "region_code"), 1)
    mid <- mean(c(bound(res, "lower_bound"), bound(res, "upper_bound")))
    expect_lt(t_ratio(mid), tcrit)
    expect_gt(t_ratio(bound(res, "upper_bound") + 1), tcrit)
  }
})


test_that("the Pygmalion example reproduces the published boundaries", {
  # Maxwell, Delaney, and Kelley (2027), Chapter 9 extension: the
  # simultaneous region of significance for the teacher-expectancy effect
  # runs from a pretest IQ of about 103.9 to about 128.7, with post-test IQ
  # averaged over the two follow-up assessments as the outcome.
  data(pygmalion)
  pygmalion$iq_post <- (pygmalion$iq_4 + pygmalion$iq_8) / 2
  fit <- lm(iq_post ~ iq_pre * treatment, data = pygmalion)

  res <- regions_of_significance(fit)
  expect_equal(bound(res, "lower_bound"), 103.8918578933, tolerance = 1e-8)
  expect_equal(bound(res, "upper_bound"), 128.6798173561, tolerance = 1e-8)
  expect_equal(round(bound(res, "lower_bound"), 1), 103.9)
  expect_equal(round(bound(res, "upper_bound"), 1), 128.7)

  # Hand computation from the fitted model, independent of the book.
  target <- hand_roots(fit, i = "treatmentBloomer",
                       j = "iq_pre:treatmentBloomer")
  expect_equal(bound(res, "lower_bound"), target[1], tolerance = 1e-10)
  expect_equal(bound(res, "upper_bound"), target[2], tolerance = 1e-10)

  # The slope difference does not itself clear the simultaneous critical
  # value, so the parabola opens downward and the groups differ *between*
  # the boundaries, both of which are inside the observed pretest range.
  expect_equal(bound(res, "region_code"), 2)
  expect_equal(bound(res, "n_boundaries"), 2)
  expect_equal(bound(res, "lower_bound_in_range"), 1)
  expect_equal(bound(res, "upper_bound_in_range"), 1)
  expect_match(attr(res, "pairs")$region, "103.9 < iq_pre < 128.7")

  # Difference line: Bloomer minus Control, as the pair label says.
  b <- coef(fit)
  expect_equal(bound(res, "difference_intercept"),
               unname(b[["treatmentBloomer"]]), tolerance = 1e-12)
  expect_equal(bound(res, "difference_slope"),
               unname(b[["iq_pre:treatmentBloomer"]]), tolerance = 1e-12)
  expect_identical(unique(res$pair), "Bloomer - Control")
})


test_that("a boundary outside the observed covariate range is flagged", {
  data(pygmalion)
  pygmalion$iq_post <- (pygmalion$iq_4 + pygmalion$iq_8) / 2
  fit <- lm(iq_post ~ iq_pre * treatment, data = pygmalion)

  res <- regions_of_significance(fit, method = "pointwise")
  # The pointwise upper boundary is about 174.7, past the largest observed
  # pretest IQ of 158, so it is an extrapolation and is marked as one.
  expect_gt(bound(res, "upper_bound"), max(pygmalion$iq_pre))
  expect_equal(bound(res, "upper_bound_in_range"), 0)
  expect_equal(bound(res, "lower_bound_in_range"), 1)
  expect_match(attr(res, "pairs")$region, "extrapolation")
})


test_that("the formula interface and the fitted model agree", {
  d <- sim_two_group()
  from_formula <- regions_of_significance(y ~ x * g, data = d)
  from_fit     <- regions_of_significance(lm(y ~ x * g, data = d))
  expect_equal(from_formula$value, from_fit$value)
  expect_identical(from_formula$pair, from_fit$pair)
  expect_warning(regions_of_significance(lm(y ~ x * g, data = d), data = d),
                 "ignored")
})


test_that("every pair of groups is reported, and matches a hand computation", {
  set.seed(113)
  n <- 150
  g <- factor(rep(c("control", "low", "high"), each = n / 3))
  x <- rnorm(n, 50, 10)
  y <- 2 + 0.5 * x + (g == "high") * (0.4 * x - 15) + rnorm(n, 0, 5)
  d <- data.frame(y, x, g)
  fit <- lm(y ~ x * g, data = d)

  res <- regions_of_significance(fit)
  # Levels are control, high, low; pairs are later level minus earlier one.
  expect_identical(unique(res$pair),
                   c("high - control", "low - control", "low - high"))

  # Hand computation for "low - high", the pair that involves no reference
  # level and so is not read off a single coefficient. With treatment
  # contrasts the coefficients are (Intercept), x, ghigh, glow, x:ghigh,
  # x:glow, so the contrast picks out glow - ghigh and x:glow - x:ghigh.
  L <- rbind(c(0, 0, -1, 1, 0, 0), c(0, 0, 0, 0, -1, 1))
  dvec <- as.numeric(L %*% coef(fit))
  Vd   <- L %*% vcov(fit) %*% t(L)
  tcrit <- sqrt(2 * qf(0.95, 2, fit$df.residual))
  A <- dvec[2]^2 - tcrit^2 * Vd[2, 2]
  B <- 2 * (dvec[1] * dvec[2] - tcrit^2 * Vd[1, 2])
  C <- dvec[1]^2 - tcrit^2 * Vd[1, 1]
  roots <- sort(c((-B - sqrt(B^2 - 4 * A * C)) / (2 * A),
                  (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)))

  expect_equal(bound(res, "lower_bound", "low - high"), roots[1],
               tolerance = 1e-10)
  expect_equal(bound(res, "upper_bound", "low - high"), roots[2],
               tolerance = 1e-10)
  expect_equal(bound(res, "difference_intercept", "low - high"), dvec[1],
               tolerance = 1e-12)
  expect_equal(bound(res, "difference_slope", "low - high"), dvec[2],
               tolerance = 1e-12)

  # The two groups generated from the same line never differ significantly.
  expect_equal(bound(res, "region_code", "low - control"), 4)
  expect_equal(bound(res, "n_boundaries", "low - control"), 0)
  expect_true(is.na(bound(res, "lower_bound", "low - control")))
  expect_true(is.na(bound(res, "lower_bound_in_range", "low - control")))
  expect_match(attr(res, "pairs")$region[2], "do not differ at any value")
})


test_that("predictors that do not interact with the group are allowed", {
  set.seed(113)
  n <- 120
  g <- factor(rep(c("a", "b"), each = n / 2))
  block <- factor(rep(c("I", "II", "III"), length.out = n))
  x <- rnorm(n, 50, 10)
  y <- 0.5 * x + (g == "b") * (0.3 * x - 12) +
    as.numeric(block) * 2 + rnorm(n, 0, 5)
  d <- data.frame(y, x, g, block)

  fit <- lm(y ~ x * g + block, data = d)
  res <- regions_of_significance(fit)
  # The block term cancels out of the group difference, so the boundaries
  # solve the same quadratic built from the group and interaction
  # coefficients (positions 3 and 4 of the coefficient vector here).
  target <- hand_roots(fit, i = "gb", j = "x:gb")
  expect_equal(bound(res, "lower_bound"), target[1], tolerance = 1e-10)
  expect_equal(bound(res, "upper_bound"), target[2], tolerance = 1e-10)
})


test_that("the degenerate linear geometry is solved correctly", {
  # When the slope difference sits exactly on the critical value the
  # quadratic collapses to a line and there is a single boundary. With
  # t_crit = 2, d0 = 1, d1 = 2, Var(d0) = Var(d1) = 1, Cov = 0:
  #   f(x) = (1 + 2x)^2 - 4(1 + x^2) = 4x - 3,
  # so the single boundary is at x = 0.75 and the groups differ above it.
  sol <- DMAR:::.regions_solve(d0 = 1, d1 = 2, v00 = 1, v01 = 0, v11 = 1,
                               tcrit = 2)
  expect_equal(sol$lower, 0.75)
  expect_true(is.na(sol$upper))
  expect_identical(sol$n, 1L)
  expect_identical(sol$code, 5L)

  # Mirror image: the difference line points the other way, so the groups
  # differ below the single boundary.
  sol2 <- DMAR:::.regions_solve(d0 = -1, d1 = 2, v00 = 1, v01 = 0, v11 = 1,
                                tcrit = 2)
  expect_equal(sol2$lower, -0.75)
  expect_identical(sol2$code, 6L)
})


test_that("the result is a tidy dmar_tbl with metadata on attributes", {
  d <- sim_two_group()
  res <- regions_of_significance(y ~ x * g, data = d)

  expect_s3_class(res, "dmar_tbl")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("pair", "term", "value") %in% names(res)))
  expect_true(is.numeric(res$value))
  expect_true(is.character(res$term))

  expect_identical(attr(res, "method"), "simultaneous")
  expect_identical(attr(res, "covariate"), "x")
  expect_identical(attr(res, "group"), "g")
  expect_identical(attr(res, "outcome"), "y")
  expect_equal(attr(res, "covariate_range"), range(d$x))
  expect_identical(attr(res, "conf_level"), 0.95)
  expect_true(is.data.frame(attr(res, "pairs")))
  expect_true(is.character(attr(res, "pairs")$region))
  expect_true(is.data.frame(attr(res, "geometry")))

  # Printing works and does not error.
  expect_output(print(res), "lower_bound")
})


test_that("regions_of_significance() rejects models it cannot interpret", {
  set.seed(113)
  n <- 80
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  g = factor(rep(c("a", "b"), each = n / 2)),
                  h = factor(rep(c("u", "v"), length.out = n)),
                  g01 = rep(c(0, 1), each = n / 2))

  # No interaction at all.
  expect_error(regions_of_significance(y ~ x + g, data = d),
               "No covariate-by-group interaction")
  # A grouping variable stored as a number is a numeric predictor to R.
  expect_error(regions_of_significance(y ~ x * g01, data = d),
               "No covariate-by-group interaction")
  # More than one numeric-by-factor interaction.
  expect_error(regions_of_significance(y ~ x * g + x * h, data = d),
               "more than one")
  # The grouping factor interacting with something else as well.
  expect_error(regions_of_significance(y ~ x * g + g:h, data = d),
               "also appears in the term")
  expect_error(regions_of_significance(y ~ x * g * h, data = d),
               "more than one")
  # Bad inputs.
  expect_error(regions_of_significance(y ~ x * g, data = d,
                                       conf_level = 1.2),
               "conf_level")
  expect_error(regions_of_significance(y ~ x * g), "data.frame")
  expect_error(regions_of_significance("not a model"), "fitted lm")
})


test_that("plot_regions_of_significance() returns a ggplot", {
  skip_if_not_installed("ggplot2")
  d <- sim_two_group()
  res <- regions_of_significance(y ~ x * g, data = d)

  p <- plot_regions_of_significance(res)
  expect_s3_class(p, "ggplot")
  # It accepts the model or the formula directly as well.
  expect_s3_class(plot_regions_of_significance(lm(y ~ x * g, data = d)),
                  "ggplot")
  expect_s3_class(plot_regions_of_significance(y ~ x * g, data = d), "ggplot")

  # The drawn band is D(x) +/- t_crit * SE[D(x)], so it touches zero
  # exactly at an in-range boundary.
  built <- ggplot2::ggplot_build(p)
  ribbon <- built$data[[1]]
  expect_true(all(c("ymin", "ymax") %in% names(ribbon)))
  expect_equal(range(ribbon$x), range(d$x), tolerance = 1e-8)

  # Vertical lines mark the in-range boundaries.
  lower <- bound(res, "lower_bound")
  upper <- bound(res, "upper_bound")
  in_range <- c(lower, upper)
  in_range <- in_range[in_range >= min(d$x) & in_range <= max(d$x)]
  vlines <- built$data[[4]]
  expect_equal(sort(vlines$xintercept), sort(in_range), tolerance = 1e-8)

  expect_error(plot_regions_of_significance(res, n_points = 2), "n_points")
  expect_error(plot_regions_of_significance(res, facet = "yes"), "facet")
})
