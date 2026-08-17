# Tests for the factorial composite power planners. The composite of a set of
# balanced factorial F tests that share one error estimate is evaluated by
# integrating the joint conditional rejection probability over the chi square
# error; a single effect must reproduce the ordinary noncentral F power, and a
# two-effect composite must match a direct simulation of the joint rejection.

test_that("a single effect reproduces ss_power_factorial_anova exactly", {
  got <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), effects = list(list(factors = 2, f = 0.25)),
    n_per_cell = 20)
  ref <- ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 2,
                                  f = 0.25, n_per_cell = 20)
  expect_equal(got$value[got$term == "composite_power"],
               ref$value[ref$term == "actual_power"], tolerance = 1e-10)
})

test_that("a single effect with a covariate reproduces ss_power_factorial_ancova", {
  got <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 4, 3), effects = list(list(factors = 1, f = 0.10)),
    covariate_R2 = 0.30, n_covariates = 2, n_per_cell = 25)
  ref <- ss_power_factorial_ancova(factor_levels = c(2, 4, 3),
                                   effect_indices = 1, f = 0.10,
                                   covariate_R2 = 0.30, n_covariates = 2,
                                   n_per_cell = 25)
  expect_equal(got$value[got$term == "composite_power"],
               ref$value[ref$term == "actual_power"], tolerance = 1e-10)
})

test_that("the ANOVA wrapper is the ANCOVA with no covariate", {
  a <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 2), n_per_cell = 30,
    effects = list(list(factors = 1, f = 0.3), list(factors = c(1, 2), f = 0.25)))
  b <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 2), n_per_cell = 30, covariate_R2 = 0, n_covariates = 0,
    effects = list(list(factors = 1, f = 0.3), list(factors = c(1, 2), f = 0.25)))
  expect_equal(a$value, b$value)
  expect_true("factor_levels" %in% names(formals(ss_power_composite_factorial_anova)))
  expect_false("covariate_R2" %in% names(formals(ss_power_composite_factorial_anova)))
})

test_that("partial_eta_squared and the equivalent f agree", {
  pes <- 0.06
  f   <- sqrt(pes / (1 - pes))
  a <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = 25,
    effects = list(list(factors = 1, partial_eta_squared = pes),
                   list(factors = 2, f = 0.2)))
  b <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = 25,
    effects = list(list(factors = 1, f = f),
                   list(factors = 2, f = 0.2)))
  expect_equal(a$value[a$term == "composite_power"],
               b$value[b$term == "composite_power"], tolerance = 1e-12)
})

test_that("prod(marginals) <= composite <= min(marginals)", {
  p <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 2, 3), n_per_cell = 15,
    effects = list(list(factors = 1, f = 0.30, label = "A"),
                   list(factors = c(1, 3), f = 0.25, label = "AxC")))
  comp <- p$value[p$term == "composite_power"]
  marg <- p$value[p$term %in% c("power_A", "power_AxC")]
  expect_gte(comp, prod(marg) - 1e-9)
  expect_lte(comp, min(marg) + 1e-9)
})

test_that("the necessary-n recommendation is minimal", {
  p <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), desired_power = 0.80,
    effects = list(list(factors = 1, f = 0.25), list(factors = 2, f = 0.20)))
  n  <- p$value[p$term == "necessary_n_per_cell"]
  at <- function(k) ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = k,
    effects = list(list(factors = 1, f = 0.25),
                   list(factors = 2, f = 0.20)))$value[3]  # composite_power row
  expect_gte(at(n), 0.80)
  expect_lt(at(n - 1), 0.80)
})

test_that("tidy() and glance() summarize the per-cell size and composite power", {
  p <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 3), desired_power = 0.80, covariate_R2 = 0.25,
    n_covariates = 1,
    effects = list(list(factors = 1, f = 0.25), list(factors = 2, f = 0.20)))
  expect_s3_class(p, "dmar_ss_power")
  expect_s3_class(p, "dmar_composite_power_factorial")
  expect_equal(generics::tidy(p)$estimate,
               p$value[p$term == "necessary_n_per_cell"])
  expect_equal(generics::tidy(p)$power, p$value[p$term == "composite_power"])
})

test_that("plot() returns a ggplot", {
  skip_if_not_installed("ggplot2")
  p <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = 30,
    effects = list(list(factors = 1, f = 0.25), list(factors = 2, f = 0.20)))
  expect_s3_class(plot(p), "ggplot")
})

test_that("inputs are validated", {
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(1, 3),
    effects = list(list(factors = 1, f = 0.2))), "at least 2")
  expect_error(ss_power_composite_factorial_ancova(factor_levels = c(2, 3),
    effects = list(list(factors = 1, f = 0.2)), covariate_R2 = 0.3,
    n_covariates = 0), "how many covariates")
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    effects = list(list(factors = 5, f = 0.2))), "indices into")
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    effects = list(list(factors = 1))), "exactly one of")
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    effects = list(list(factors = 1, f = 0.2), list(factors = 1, f = 0.3))),
    "named more than once")
})

test_that("the two-effect composite matches a direct simulation", {
  skip_on_cran()
  # 2 x 2 ANOVA. Cell means m + a A_i + b B_j + g A_i B_j with A_i, B_j in
  # {-1, +1} give Cohen's f = |a| for main A and |g| for the A:B interaction at
  # sigma = 1, so calling the planner with those f reproduces the design.
  set.seed(113)
  m <- 10; a <- 0.35; b <- 0.20; g <- 0.30; ncell <- 30
  analytic <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 2), n_per_cell = ncell,
    effects = list(list(factors = 1,       f = a, label = "A"),
                   list(factors = c(1, 2), f = g, label = "AxB")))
  comp <- analytic$value[analytic$term == "composite_power"]

  Ai <- c(-1, -1, 1, 1); Bj <- c(-1, 1, -1, 1)
  mu <- m + a * Ai + b * Bj + g * Ai * Bj
  N  <- 4 * ncell
  fac_A <- factor(rep(c(1, 1, 2, 2), each = ncell))
  fac_B <- factor(rep(c(1, 2, 1, 2), each = ncell))
  G <- 4000
  reject_both <- replicate(G, {
    y  <- rnorm(N, mean = rep(mu, each = ncell), sd = 1)
    pv <- anova(lm(y ~ fac_A * fac_B))[["Pr(>F)"]]
    as.numeric(pv[1] < 0.05 & pv[3] < 0.05)   # A and A:B both reject
  })
  mc <- mean(reject_both); se <- sqrt(mc * (1 - mc) / G)
  expect_lt(abs(comp - mc), 4 * se)
})

test_that("effects from cell means match the equivalent effect sizes", {
  # 2 x 2 means m + a A + b B + g AB give f_A = |a|/sigma, f_AB = |g|/sigma.
  m <- 10; a <- 0.35; b <- 0.20; g <- 0.30; sigma <- 2
  M <- outer(c(-1, 1), c(-1, 1), function(ai, bj) m + a * ai + b * bj + g * ai * bj)
  from_means <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 2), n_per_cell = 50, means = M, sigma = sigma,
    effects = list(list(factors = 1, label = "A"),
                   list(factors = c(1, 2), label = "AxB")))
  from_es <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 2), n_per_cell = 50,
    effects = list(list(factors = 1, f = a / sigma, label = "A"),
                   list(factors = c(1, 2), f = g / sigma, label = "AxB")))
  expect_equal(from_means$value[from_means$term == "composite_power"],
               from_es$value[from_es$term == "composite_power"], tolerance = 1e-10)
  expect_equal(from_means$value[from_means$term == "f_A"], a / sigma, tolerance = 1e-12)
  expect_equal(from_means$value[from_means$term == "f_AxB"], g / sigma, tolerance = 1e-12)
})

test_that("a vector of means equals the array of means", {
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)  # 2 x 3
  a <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = 30, means = M, sigma = 4,
    effects = list(list(factors = 1), list(factors = 2)))
  b <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = 30, means = as.vector(M), sigma = 4,
    effects = list(list(factors = 1), list(factors = 2)))
  expect_equal(a$value, b$value)
})

test_that("plot() draws the mean pattern (with error bars) when means are given", {
  skip_if_not_installed("ggplot2")
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  p <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), n_per_cell = 30, means = M, sigma = 4,
    effects = list(list(factors = 1), list(factors = 2)))
  gg <- plot(p)
  expect_s3_class(gg, "ggplot")
  expect_true(any(vapply(gg$layers,
    function(L) inherits(L$geom, "GeomErrorbar"), logical(1))))
})

test_that("the means interface is validated", {
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    means = M, effects = list(list(factors = 1))), "positive 'sigma'")
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    means = M, sigma = 4, effects = list(list(factors = 1, f = 0.2))),
    "do not also give")
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    sigma = 4, effects = list(list(factors = 1, f = 0.2))), "only used with")
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 2),
    means = M, sigma = 4, effects = list(list(factors = 1))), "must match")
  # An effect with no variation in the means (factor 1 flat) cannot be planned.
  flat <- matrix(c(1, 2, 3, 1, 2, 3), nrow = 2, byrow = TRUE)   # varies in factor 2 only
  expect_error(ss_power_composite_factorial_anova(factor_levels = c(2, 3),
    means = flat, sigma = 1, effects = list(list(factors = 1, label = "A")),
    n_per_cell = 20), "do not vary")
})
