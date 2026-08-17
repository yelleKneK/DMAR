# Tests for the no-covariate ANOVA composite entry point. It forwards to the
# factorial ANOVA planner for both one-way and factorial designs and by both the
# effect-size and population-means interfaces, but drops the covariate rows the
# ANCOVA engine echoes so that its output surfaces no covariate.

test_that("ss_power_composite_anova() forwards to the factorial ANOVA planner", {
  eff <- list(list(factors = 1, f = 0.25), list(factors = 2, f = 0.20))
  gen <- ss_power_composite_anova(factor_levels = c(2, 3), effects = eff,
                                  n_per_cell = 30)
  dir <- ss_power_composite_factorial_anova(factor_levels = c(2, 3),
                                            effects = eff, n_per_cell = 30)
  # Identical to the factorial ANOVA planner save for the covariate rows, which
  # the ANOVA entry point drops.
  expect_false(any(gen$term %in% c("covariate_R2", "n_covariates")))
  shared <- !dir$term %in% c("covariate_R2", "n_covariates")
  expect_equal(gen$term, dir$term[shared])
  expect_equal(gen$value, dir$value[shared])
  expect_s3_class(gen, "dmar_composite_power_factorial")
  expect_s3_class(gen, "dmar_ss_power")
})

test_that("ss_power_composite_anova() exposes no covariate arguments", {
  # The no-covariate entry point must not surface a covariate in its interface;
  # a design with a covariate is directed to ss_power_composite_ancova instead.
  fmls <- names(formals(ss_power_composite_anova))
  expect_false(any(grepl("cov", fmls, ignore.case = TRUE)))
  expect_setequal(fmls, c("factor_levels", "effects", "means", "sigma",
                          "desired_power", "alpha_level", "n_per_cell"))
})

test_that("ss_power_composite_anova() plans an a-group one-way design", {
  res <- ss_power_composite_anova(
    factor_levels = 4,
    effects = list(list(factors = 1, f = 0.30)),
    desired_power = 0.80)
  # A one-way composite of a single effect is the ordinary one-way ANOVA power,
  # so the per-group size matches the dedicated one-way planner.
  ref <- ss_power_one_way_anova(a = 4, f = 0.30, desired_power = 0.80)
  expect_equal(res$value[res$term == "necessary_n_per_cell"],
               ref$value[ref$term == "n_per_group"], tolerance = 1)
})

test_that("ss_power_composite_anova() accepts the population-means interface", {
  M <- matrix(c(10, 12, 11, 13, 12, 16), nrow = 2, byrow = TRUE)
  gen <- ss_power_composite_anova(
    factor_levels = c(2, 3), means = M, sigma = 4,
    effects = list(list(factors = 1, label = "A"),
                   list(factors = 2, label = "B")),
    n_per_cell = 30)
  dir <- ss_power_composite_factorial_anova(
    factor_levels = c(2, 3), means = M, sigma = 4,
    effects = list(list(factors = 1, label = "A"),
                   list(factors = 2, label = "B")),
    n_per_cell = 30)
  shared <- !dir$term %in% c("covariate_R2", "n_covariates")
  expect_equal(gen$value, dir$value[shared])
  # The means attribute survives the covariate-row strip, so plot() still draws
  # the cell-mean pattern.
  expect_equal(attr(gen, "means"), M)
  expect_s3_class(plot(gen), "ggplot")
})
