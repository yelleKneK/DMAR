# Tests for the general one-way / factorial ANCOVA composite dispatcher. It must
# forward each slopes mode to the matching factorial planner unchanged, reproduce
# the two-group planner in the one-factor two-level heterogeneous case, and reject
# arguments that belong to the other slopes mode.

test_that("slopes = \"homogeneous\" forwards to the factorial ANCOVA unchanged", {
  eff <- list(list(factors = 1, f = 0.25), list(factors = 2, f = 0.20))
  gen <- ss_power_composite_ancova(
    factor_levels = c(2, 3), effects = eff,
    covariate_R2 = 0.25, n_covariates = 1, n_per_cell = 30)
  dir <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 3), effects = eff,
    covariate_R2 = 0.25, n_covariates = 1, n_per_cell = 30)
  expect_equal(gen$term, dir$term)
  expect_equal(gen$value, dir$value)
  expect_s3_class(gen, "dmar_composite_power_factorial")
})

test_that("slopes = \"heterogeneous\" forwards to the het planner unchanged", {
  eff <- list(list(type = "mean", factors = 1, f = 0.30),
              list(type = "covariate", f = 0.40),
              list(type = "slope", factors = 1, f = 0.20))
  gen <- ss_power_composite_ancova(
    factor_levels = 4, slopes = "heterogeneous", effects = eff, n_per_cell = 25)
  dir <- ss_power_composite_factorial_ancova_het(
    factor_levels = 4, effects = eff, n_per_cell = 25)
  expect_equal(gen$term, dir$term)
  expect_equal(gen$value, dir$value)
  expect_s3_class(gen, "dmar_composite_power_factorial_het")
})

test_that("the one-factor two-level heterogeneous case reproduces the two-group
           planner", {
  gen <- suppressWarnings(ss_power_composite_ancova(
    factor_levels = 2, slopes = "heterogeneous",
    means = c(-0.25, 0.25), correlations = c(0.1, 0.5), sigma = 1,
    n_per_cell = 100,
    effects = list(list(type = "mean", factors = 1),
                   list(type = "covariate"),
                   list(type = "slope", factors = 1))))
  two <- suppressWarnings(ss_power_composite_ancova_2group(
    smd = 0.5, rho = c(0.1, 0.5), sigma = 1, n = 100,
    composite_terms = c("group", "covariate", "group_by_covariate")))
  # The correlations differ in absolute value, so both tables report the
  # composite under its approximate name. Reading the exact name here would
  # compare two zero-length vectors and pass without checking anything.
  got  <- gen$value[gen$term == "approximate_composite_power"]
  want <- two$value[two$term == "approximate_composite_power"]
  expect_length(got, 1L)
  expect_length(want, 1L)
  expect_equal(got, want, tolerance = 1e-10)
})

test_that("the general dispatcher supports an a-group one-way ANCOVA", {
  # More than two groups, which the two-group planner cannot express.
  res <- ss_power_composite_ancova(
    factor_levels = 5, slopes = "heterogeneous",
    effects = list(list(type = "mean", factors = 1, f = 0.30),
                   list(type = "covariate", f = 0.40),
                   list(type = "slope", factors = 1, f = 0.20)),
    desired_power = 0.80)
  n <- res$value[res$term == "necessary_n_per_cell"]
  expect_true(is.finite(n) && n >= 2)
  expect_gte(res$value[res$term == "composite_power"], 0.80)
  expect_equal(res$value[res$term == "cells"], 5)
})

test_that("the general dispatcher rejects arguments from the other slopes mode", {
  eff_h <- list(list(factors = 1, f = 0.25))
  eff_t <- list(list(type = "mean", factors = 1, f = 0.25),
                list(type = "covariate", f = 0.30))
  # Heterogeneous-only argument in a homogeneous call.
  expect_error(
    ss_power_composite_ancova(factor_levels = 3, effects = eff_h,
                              correlations = c(0.2, 0.3, 0.4)),
    "correlations")
  expect_error(
    ss_power_composite_ancova(factor_levels = 3, effects = eff_h, sd_cov = 2),
    "sd_cov")
  # Homogeneous-only argument in a heterogeneous call.
  expect_error(
    ss_power_composite_ancova(factor_levels = 3, slopes = "heterogeneous",
                              effects = eff_t, covariate_R2 = 0.2),
    "covariate_R2")
  expect_error(
    ss_power_composite_ancova(factor_levels = 3, slopes = "heterogeneous",
                              effects = eff_t, n_covariates = 1),
    "n_covariates")
})

test_that("the general dispatcher accepts the population-means interface", {
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  gen <- ss_power_composite_ancova(
    factor_levels = c(2, 3), means = M, sigma = 4,
    effects = list(list(factors = 1, label = "A"),
                   list(factors = 2, label = "B")),
    n_per_cell = 30)
  dir <- ss_power_composite_factorial_ancova(
    factor_levels = c(2, 3), means = M, sigma = 4,
    effects = list(list(factors = 1, label = "A"),
                   list(factors = 2, label = "B")),
    n_per_cell = 30)
  expect_equal(gen$value, dir$value)
})
