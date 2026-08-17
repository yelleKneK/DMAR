# Tests for the heterogeneous-slope factorial ANCOVA composite planner. The
# one-factor two-level case must reproduce the two-group composite ANCOVA, the
# population-values and effect-size interfaces must agree, and a three-effect
# composite must match a direct simulation of the full heterogeneous-slope model.

test_that("the one-factor two-level case reproduces ss_power_composite_ancova_2group", {
  het <- suppressWarnings(ss_power_composite_factorial_ancova_het(
    factor_levels = 2, means = c(-0.25, 0.25), correlations = c(0.1, 0.5),
    sigma = 1, n_per_cell = 100,
    effects = list(list(type = "mean", factors = 1),
                   list(type = "covariate"),
                   list(type = "slope", factors = 1))))
  two <- suppressWarnings(ss_power_composite_ancova_2group(
    smd = 0.5, rho = c(0.1, 0.5), sigma = 1, n = 100,
    composite_terms = c("group", "covariate", "group_by_covariate")))
  # These correlations differ in absolute value, so both planners report the
  # composite under its approximate name; naming the row keeps the comparison
  # from succeeding on two zero-length vectors.
  got  <- het$value[het$term == "approximate_composite_power"]
  want <- two$value[two$term == "approximate_composite_power"]
  expect_length(got, 1L)
  expect_length(want, 1L)
  expect_equal(got, want, tolerance = 1e-10)
})

test_that("the population-values and effect-size interfaces agree", {
  M <- matrix(c(10, 12, 11, 13), nrow = 2, byrow = TRUE)
  R <- matrix(c(0.55, 0.55, 0.15, 0.20), nrow = 2, byrow = TRUE)
  pop <- suppressWarnings(ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), means = M, correlations = R, sigma = 4, sd_cov = 2,
    n_per_cell = 40,
    effects = list(list(type = "mean", factors = 1),
                   list(type = "covariate"),
                   list(type = "slope", factors = 1))))
  fs <- ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), n_per_cell = 40,
    effects = list(list(type = "mean", factors = 1,
                        f = pop$value[pop$term == "f_1"]),
                   list(type = "covariate",
                        f = pop$value[pop$term == "f_covariate"]),
                   list(type = "slope", factors = 1,
                        f = pop$value[pop$term == "f_cov_x_1"])))
  # The population values have cell correlations of differing absolute value,
  # so that table is relabeled; the effect size interface says nothing about
  # them, so its table keeps the exact names. The powers themselves agree.
  got  <- pop$value[pop$term == "approximate_composite_power"]
  want <- fs$value[fs$term == "composite_power"]
  expect_length(got, 1L)
  expect_length(want, 1L)
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("sd_cov does not change any power (a correlation is scale free)", {
  M <- matrix(c(10, 12, 11, 13), nrow = 2, byrow = TRUE)
  R <- matrix(c(0.55, 0.55, 0.15, 0.20), nrow = 2, byrow = TRUE)
  base <- function(sd) suppressWarnings(ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), means = M, correlations = R, sigma = 4, sd_cov = sd,
    n_per_cell = 40,
    effects = list(list(type = "covariate"),
                   list(type = "slope", factors = 1))))$value
  expect_equal(base(1), base(7))
})

test_that("prod(marginals) <= composite <= min(marginals)", {
  p <- ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 3), n_per_cell = 25,
    effects = list(list(type = "mean", factors = 2, f = 0.25),
                   list(type = "covariate", f = 0.30),
                   list(type = "slope", factors = 2, f = 0.20)))
  comp <- p$value[p$term == "composite_power"]
  marg <- p$value[p$term %in% c("power_2", "power_covariate", "power_cov_x_2")]
  expect_gte(comp, prod(marg) - 1e-9)
  expect_lte(comp, min(marg) + 1e-9)
})

test_that("tidy() summarizes the per-cell size and composite power", {
  p <- ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), desired_power = 0.80,
    effects = list(list(type = "mean", factors = 1, f = 0.25),
                   list(type = "covariate", f = 0.40)))
  expect_s3_class(p, "dmar_ss_power")
  expect_s3_class(p, "dmar_composite_power_factorial_het")
  expect_equal(generics::tidy(p)$estimate,
               p$value[p$term == "necessary_n_per_cell"])
  expect_equal(generics::tidy(p)$power, p$value[p$term == "composite_power"])
})

test_that("plot() returns a ggplot on both interfaces", {
  skip_if_not_installed("ggplot2")
  M <- matrix(c(10, 12, 11, 13), nrow = 2, byrow = TRUE)
  R <- matrix(c(0.55, 0.55, 0.15, 0.20), nrow = 2, byrow = TRUE)
  pop <- suppressWarnings(ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), means = M, correlations = R, sigma = 4,
    n_per_cell = 40, effects = list(list(type = "mean", factors = 1),
                                    list(type = "slope", factors = 1))))
  gg <- plot(pop)
  expect_s3_class(gg, "ggplot")
  expect_true(any(vapply(gg$layers,
    function(L) inherits(L$geom, "GeomLine"), logical(1))))
  fs <- ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), n_per_cell = 40,
    effects = list(list(type = "mean", factors = 1, f = 0.25),
                   list(type = "covariate", f = 0.40)))
  expect_s3_class(plot(fs), "ggplot")
})

test_that("the heterogeneous-slope taxonomy is validated", {
  expect_error(ss_power_composite_factorial_ancova_het(factor_levels = c(2, 2),
    effects = list(list(type = "wrong", factors = 1, f = 0.2))), "must be one of")
  expect_error(ss_power_composite_factorial_ancova_het(factor_levels = c(2, 2),
    effects = list(list(type = "covariate", factors = 1, f = 0.2))),
    "spans no factors")
  expect_error(ss_power_composite_factorial_ancova_het(factor_levels = c(2, 2),
    means = matrix(0, 2, 2), sigma = 1,
    effects = list(list(type = "mean", factors = 1))), "needs 'correlations'")
  expect_error(ss_power_composite_factorial_ancova_het(factor_levels = c(2, 2),
    correlations = matrix(0.3, 2, 2),
    effects = list(list(type = "covariate"))), "positive 'sigma'")
  expect_error(ss_power_composite_factorial_ancova_het(factor_levels = c(2, 2),
    correlations = matrix(1.2, 2, 2), sigma = 1,
    effects = list(list(type = "covariate"))), "in \\(-1, 1\\)")
})

test_that("a three-effect heterogeneous-slope composite matches simulation", {
  skip_on_cran()
  skip_if_not_installed("car")
  set.seed(113)
  sd_cov <- 1.5; sigma <- 2; ncell <- 90; cells <- 4; N <- ncell * cells
  muM  <- matrix(c(10, 10.4, 10.2, 10.6), nrow = 2, byrow = TRUE)
  rhoM <- matrix(c(0.55, 0.55, 0.15, 0.20), nrow = 2, byrow = TRUE)
  analytic <- suppressWarnings(ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), means = muM, correlations = rhoM, sigma = sigma,
    sd_cov = sd_cov, n_per_cell = ncell,
    effects = list(list(type = "mean", factors = 1, label = "A"),
                   list(type = "covariate"),
                   list(type = "slope", factors = 1, label = "cov_x_1"))))
  # The cell correlations differ in absolute value, so the composite is
  # reported under its approximate name.
  comp <- analytic$value[analytic$term == "approximate_composite_power"]
  expect_length(comp, 1L)

  grid <- expand.grid(A = 1:2, B = 1:2)         # factor 1 fastest, matches as.vector
  mu_v <- as.vector(muM); beta_v <- as.vector(rhoM) * sigma / sd_cov
  rho_v <- as.vector(rhoM)
  G <- 2500
  rej <- t(replicate(G, {
    A <- factor(rep(grid$A, each = ncell)); B <- factor(rep(grid$B, each = ncell))
    cell <- rep(seq_len(cells), each = ncell)
    X <- rnorm(N, 0, sd_cov)
    Y <- mu_v[cell] + beta_v[cell] * X +
      rnorm(N, 0, sigma * sqrt(1 - rho_v[cell]^2))
    contrasts(A) <- stats::contr.sum(2); contrasts(B) <- stats::contr.sum(2)
    an <- car::Anova(lm(Y ~ A * B * X), type = 3)
    p  <- an[, "Pr(>F)"]; rn <- rownames(an)
    c(p[rn == "A"], p[rn == "X"], p[rn == "A:X"]) < 0.05
  }))
  mc <- mean(rej[, 1] & rej[, 2] & rej[, 3])
  se <- sqrt(mc * (1 - mc) / G)
  expect_lt(abs(comp - mc), 4 * se + 0.02)   # slack for the O(1/N) conditioning
})

test_that("cell correlations differing in absolute value warn (HIGH-01)", {
  cm  <- c(-0.25, 0.25)
  expect_warning(
    ss_power_composite_factorial_ancova_het(
      factor_levels = 2,
      effects = list(list(type = "mean", factors = 1)),
      means = cm, correlations = c(0.1, 0.5), sigma = 1, n_per_cell = 30),
    "approximation")
  # Cell correlations of equal absolute value leave the residual variances
  # equal, which is the exact case, and stay silent.
  expect_no_warning(
    ss_power_composite_factorial_ancova_het(
      factor_levels = 2,
      effects = list(list(type = "mean", factors = 1)),
      means = cm, correlations = c(-0.4, 0.4), sigma = 1, n_per_cell = 30))
})

test_that("cell correlations differing in absolute value relabel the powers", {
  cm <- c(-0.25, 0.25)
  eff <- list(list(type = "mean", factors = 1), list(type = "covariate"))
  approx <- suppressWarnings(ss_power_composite_factorial_ancova_het(
    factor_levels = 2, effects = eff, means = cm, correlations = c(0.1, 0.5),
    sigma = 1, desired_power = 0.80))
  expect_true(all(c("approximate_n_per_cell", "approximate_N",
                    "approximate_composite_power", "approximate_power_1",
                    "approximate_power_covariate") %in% approx$term))
  expect_false(any(c("necessary_n_per_cell", "necessary_N", "composite_power")
                   %in% approx$term))
  expect_true(attr(approx, "approximate"))
  # The effect sizes, degrees of freedom, and noncentralities are not powers and
  # keep their names.
  expect_true(all(c("f_1", "df_1", "noncentral_parm_1", "cells", "alpha_level",
                    "desired_power") %in% approx$term))
  # tidy() and glance() read fixed row names, so they have to be taught the
  # relabeled ones or they silently report nothing.
  td <- generics::tidy(approx)
  expect_equal(td$estimate, approx$value[approx$term == "approximate_n_per_cell"])
  expect_equal(td$power,
               approx$value[approx$term == "approximate_composite_power"])
  gl <- generics::glance(approx)
  expect_equal(gl$estimate, td$estimate)
  expect_equal(gl$power, td$power)

  # Equal absolute cell correlations are exact and unchanged. Their average is
  # zero here, so the composite names the mean effect and the slope
  # heterogeneity, both of which vary.
  exact <- ss_power_composite_factorial_ancova_het(
    factor_levels = 2,
    effects = list(list(type = "mean", factors = 1),
                   list(type = "slope", factors = 1)),
    means = cm, correlations = c(-0.4, 0.4),
    sigma = 1, desired_power = 0.80)
  expect_true(all(c("necessary_n_per_cell", "necessary_N", "composite_power",
                    "power_1", "power_cov_x_1") %in% exact$term))
  expect_false(any(grepl("^approximate_", exact$term)))
  expect_false(attr(exact, "approximate"))
  expect_equal(generics::tidy(exact)$estimate,
               exact$value[exact$term == "necessary_n_per_cell"])
})

test_that("the effect size interface cannot see the cell correlations", {
  # Effect sizes state the effects directly and say nothing about the cells'
  # residual variances, so the function has nothing to detect the approximation
  # from and labels the table exact. The help page says so; this pins it down.
  fs <- ss_power_composite_factorial_ancova_het(
    factor_levels = c(2, 2), n_per_cell = 40,
    effects = list(list(type = "mean", factors = 1, f = 0.25),
                   list(type = "covariate", f = 0.40)))
  expect_false(attr(fs, "approximate"))
  expect_true("composite_power" %in% fs$term)
})
