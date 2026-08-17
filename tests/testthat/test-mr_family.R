## mr_smd() / mr_cv() -- multi-stage (Robbins-Monro-style) sample size estimators.

test_that("mr_smd() returns a tidy data.frame with numeric value column", {
  res <- mr_smd(structural_cost = 0.0001, epsilon = 0.05, d = 0.5,
                n = 30, sampling_cost = 0.001)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("risk", "n1", "n2", "d") %in% res$term))
  expect_true(is.numeric(res$value))
  # is_satisfied is non-numeric metadata and lives on an attribute,
  # not as a row with a stringly-typed value.
  expect_false("is_satisfied" %in% res$term)
  expect_true(is.logical(attr(res, "is_satisfied")))
  expect_true(is.numeric(attr(res, "criterion")))
})

test_that("mr_smd() risk equals its closed form at a pinned design", {
  # A = structural_cost / epsilon^2 = 0.0001 / 0.05^2 = 0.04.
  # Risk = A * (2/n + d^2 / (4n)) + sampling_cost * 2n
  #      = 0.04 * (2/30 + 0.25/120) + 0.001 * 60 = 0.06275.
  res <- mr_smd(structural_cost = 0.0001, epsilon = 0.05, d = 0.5,
                n = 30, sampling_cost = 0.001)
  A <- 0.0001 / 0.05^2
  risk_hand <- A * (2 / 30 + 0.5^2 / (4 * 30)) + 0.001 * 60
  expect_equal(res$value[res$term == "risk"], risk_hand,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "risk"], 0.06275, tolerance = 1e-10)
})

test_that("mr_smd() accepts A directly and skips epsilon/structural_cost", {
  res <- mr_smd(A = 1, d = 0.5, n = 30, sampling_cost = 0.001)
  expect_s3_class(res, "data.frame")
  expect_true("risk" %in% res$term)
})

test_that("mr_smd() raises an error when A and epsilon are both supplied", {
  expect_error(
    mr_smd(A = 1, epsilon = 0.05, d = 0.5, n = 30, sampling_cost = 0.001),
    "specified 'A'"
  )
})

test_that("mr_cv() pilot mode returns a small pilot sample size", {
  set.seed(113)
  x <- rnorm(40, mean = 10, sd = 2)
  res <- mr_cv(data = x, A = 1, sampling_cost = 0.001, pilot = TRUE, m0 = 4)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("pilot_ss" %in% res$term)
  expect_gt(res$value[res$term == "pilot_ss"], 0)
})


test_that("mr_cv() non-pilot mode returns numeric value column with attributes", {
  set.seed(113)
  x <- rlnorm(60, meanlog = 1, sdlog = 0.3)
  res <- mr_cv(data = x, A = 1, sampling_cost = 0.001, gamma = 0.5,
               pilot = FALSE, m0 = 60, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true(is.numeric(res$value))
  expect_false("is_satisfied" %in% res$term)
  expect_true(is.logical(attr(res, "is_satisfied")))
  expect_true(is.numeric(attr(res, "criterion")))
})


test_that("mr_cv() verbose mode keeps numeric value column with extra term rows", {
  set.seed(113)
  x <- rlnorm(60)
  res <- mr_cv(data = x, A = 1, sampling_cost = 0.001, gamma = 0.5,
               pilot = FALSE, m0 = 60, verbose = TRUE)
  expect_true(is.numeric(res$value))
  expect_true(all(c("risk", "n", "cv", "v2_n", "mu_3_n", "mu_4_n",
                    "criterion") %in% res$term))
})
