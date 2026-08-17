# The closed-form width the planner inverts, recomputed independently:
# Wald width from the delta method standard error, with the b-hat
# variance carrying the 1/(1 - a^2) variance inflation factor for the
# correlation between M and X in the Y equation.
.wald_width <- function(n, a, b, conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  var_a <- (1 - a^2) / (n - 2)
  var_b <- (1 - b^2) / ((n - 3) * (1 - a^2))
  2 * z * sqrt(a^2 * var_b + b^2 * var_a)
}

test_that("ss_aipe_indirect_effect() closed form returns documented rows and the minimal N", {
  res <- ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20)
  expect_identical(res$term,
                   c("necessary_N", "expected_width", "a", "b", "ab",
                     "width_target", "conf_level"))
  expect_identical(attr(res, "ci_method"), "closed_form")
  N <- res$value[res$term == "necessary_N"]
  expect_equal(N, 116)
  expect_equal(res$value[res$term == "expected_width"],
               .wald_width(116, 0.40, 0.40), tolerance = 1e-12)
  # Minimality: the returned N meets the target and N - 1 does not.
  expect_lte(.wald_width(N, 0.40, 0.40), 0.20)
  expect_gt(.wald_width(N - 1, 0.40, 0.40), 0.20)
})

test_that("ss_aipe_indirect_effect() closed form carries the variance inflation factor for a", {
  # The b-hat variance grows with |a| through 1/(1 - a^2), so a larger
  # |a| with everything else fixed must not be planned as if b-hat were
  # estimated from an orthogonal design. These anchors also pin the
  # correction of the earlier planning variance, which omitted the
  # factor and returned N = 159 and N = 339 here.
  res_6 <- ss_aipe_indirect_effect(a = 0.60, b = 0.40, width = 0.20)
  expect_equal(res_6$value[res_6$term == "necessary_N"], 224)
  res_7 <- ss_aipe_indirect_effect(a = 0.70, b = 0.30, width = 0.15)
  expect_equal(res_7$value[res_7$term == "necessary_N"], 632)
})

test_that("ss_aipe_indirect_effect() closed form delivers the target width in simulation", {
  skip_on_cran()
  res <- ss_aipe_indirect_effect(a = 0.60, b = 0.40, width = 0.20)
  N <- res$value[res$term == "necessary_N"]
  set.seed(42)
  sens <- ss_aipe_indirect_effect_sensitivity(
    true_a = 0.60, true_b = 0.40, specified_N = N, width = 0.20,
    method = "closed_form", G = 400, print_iter = FALSE)
  med <- sens$value[sens$term == "median_ci_width"]
  expect_lt(abs(med - 0.20), 0.20 * 0.04)
})

test_that("ss_aipe_indirect_effect() Monte Carlo path is reproducible and coherent with the closed form", {
  skip_on_cran()
  res <- ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.30,
                                 method = "monte_carlo",
                                 G = 200, B = 1000, seed = 113)
  expect_identical(attr(res, "ci_method"), "monte_carlo")
  expect_identical(res$term,
                   c("necessary_N", "expected_width", "a", "b", "ab",
                     "width_target", "conf_level"))
  expect_lte(res$value[res$term == "expected_width"], 0.30)
  # Same seed, same plan.
  res_2 <- ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.30,
                                   method = "monte_carlo",
                                   G = 200, B = 1000, seed = 113)
  expect_identical(res$value, res_2$value)
  # The Monte Carlo interval is a little wider than the Wald interval,
  # so its plan sits at or modestly above the closed-form answer.
  N_cf <- ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.30)
  N_cf <- N_cf$value[N_cf$term == "necessary_N"]
  N_mc <- res$value[res$term == "necessary_N"]
  expect_gte(N_mc, N_cf)
  expect_lte(N_mc, N_cf + 20)
})

test_that("ss_aipe_indirect_effect() Monte Carlo plan delivers the target median width", {
  skip_on_cran()
  res <- ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.30,
                                 method = "monte_carlo",
                                 G = 200, B = 1000, seed = 113)
  N <- res$value[res$term == "necessary_N"]
  set.seed(42)
  sens <- ss_aipe_indirect_effect_sensitivity(
    true_a = 0.40, true_b = 0.40, specified_N = N, width = 0.30,
    method = "monte_carlo", B = 1000, G = 400, print_iter = FALSE)
  med <- sens$value[sens$term == "median_ci_width"]
  expect_lt(abs(med - 0.30), 0.30 * 0.05)
})

test_that("ss_aipe_indirect_effect() restores the caller's random number generator state", {
  skip_on_cran()
  set.seed(1)
  reference <- stats::runif(1)
  set.seed(1)
  invisible(ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.35,
                                    method = "monte_carlo",
                                    G = 20, B = 200, seed = 5))
  expect_identical(stats::runif(1), reference)
})

test_that("ss_aipe_indirect_effect() rejects bad inputs", {
  expect_error(ss_aipe_indirect_effect(a = 1.2, b = 0.4, width = 0.2),
               "in \\(-1, 1\\)")
  expect_error(ss_aipe_indirect_effect(a = 0.4, b = 0.4, width = 0),
               "positive")
  expect_error(ss_aipe_indirect_effect(a = 0.4, b = 0.4, width = 0.2,
                                       method = "sobel"),
               "closed_form")
  expect_error(ss_aipe_indirect_effect(a = 0.4, b = 0.4, width = 0.2,
                                       G = 5),
               "'G'")
  expect_error(ss_aipe_indirect_effect(a = 0.4, b = 0.4, width = 0.2,
                                       B = 50),
               "'B'")
  expect_error(ss_aipe_indirect_effect(a = 0.4, b = 0.4, width = 0.001,
                                       n_max = 500),
               "n_max")
})
