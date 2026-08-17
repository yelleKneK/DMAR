test_that("output has the documented long-format schema", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 10, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1
  )
  expect_s3_class(d, "data.frame")
  expect_identical(names(d),
                   c("id", "population", "occasion", "target_time", "time",
                     "true_score", "y"))
  expect_equal(nrow(d), 10 * 5)
  expect_s3_class(d$id, "factor")
  expect_s3_class(d$population, "factor")
  expect_equal(nlevels(d$id), 10L)
  expect_equal(nlevels(d$population), 1L)
  expect_type(d$y, "double")
  expect_equal(sort(unique(d$occasion)), 1:5)
})

test_that("attributes record the error variance, per-occasion reliability, T, and order", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 10, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 2
  )
  expect_equal(attr(d, "error_variance"), 2)
  expect_equal(attr(d, "polynomial_order"), 1L)
  expect_equal(dim(attr(d, "random_covariance")), c(2L, 2L))
  rbo <- attr(d, "reliability_by_occasion")
  expect_length(rbo, 5L)
  expect_true(all(rbo > 0 & rbo < 1))
})

test_that("a list of coefficient vectors yields one group per element", {
  set.seed(113)
  two <- simulate_longitudinal_polynomial(
    n = c(8, 12), target_times = 0:4,
    fixed_coefficients = list(control = c(10, 1.0), treatment = c(10, 1.8)),
    random_variances = c(4, 0.25), error_variance = 1
  )
  expect_equal(nlevels(two$population), 2L)
  expect_equal(as.integer(table(two$population)), c(8L * 5L, 12L * 5L))
  expect_equal(nlevels(two$id), 20L)
})

test_that("order 0 is a flat line: the true score is constant within subject", {
  set.seed(113)
  flat <- simulate_longitudinal_polynomial(
    n = 15, target_times = 0:4, fixed_coefficients = 5,
    random_variances = 2, error_variance = 1
  )
  expect_equal(attr(flat, "polynomial_order"), 0L)
  within_sd <- tapply(flat$true_score, flat$id, stats::sd)
  expect_true(all(within_sd < 1e-10))
})

test_that("with no random effects, no jitter, and no error the data are the exact polynomial", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 3, target_times = c(0, 1, 2, 3),
    fixed_coefficients = c(2, -1, 0.5),   # 2 - t + 0.5 t^2
    random_variances = 0, error_variance = 0
  )
  expected <- 2 - d$time + 0.5 * d$time^2
  expect_equal(d$true_score, expected)
  expect_equal(d$y, expected)                       # error_variance = 0 => y == true_score
  # Every subject is identical (fixed trajectory, no randomness).
  by_subj <- split(d$y, d$id)
  expect_true(all(vapply(by_subj, function(z) max(abs(z - by_subj[[1]])), 0) == 0))
})

test_that("error_variance = 0 makes y identical to the true score even with random effects", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 10, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 0
  )
  expect_equal(d$y, d$true_score)
})

test_that("timing jitter perturbs the actual time but leaves the nominal target", {
  set.seed(113)
  no_jit <- simulate_longitudinal_polynomial(
    n = 10, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1, timing_sd = 0
  )
  expect_equal(no_jit$time, no_jit$target_time)

  set.seed(113)
  jit <- simulate_longitudinal_polynomial(
    n = 10, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1, timing_sd = 0.2
  )
  expect_false(isTRUE(all.equal(jit$time, jit$target_time)))
  # The nominal schedule is preserved exactly regardless of jitter.
  expect_equal(unique(jit$target_time), as.numeric(0:4))
})

test_that("the reliability route hits the requested AVERAGE per-occasion reliability", {
  set.seed(113)
  target <- 0.75
  d <- simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), reliability = target
  )
  rbo <- attr(d, "reliability_by_occasion")
  expect_equal(mean(rbo), target, tolerance = 1e-6)
  # Reconstruct from the solved error variance and the implied true-score
  # variance c(t)'T c(t): the per-occasion reliabilities must match.
  Tcov  <- attr(d, "random_covariance")
  s2e   <- attr(d, "error_variance")
  basis <- outer(0:4, 0:1, `^`)
  tv    <- rowSums((basis %*% Tcov) * basis)
  expect_equal(unname(rbo), unname(tv / (tv + s2e)), tolerance = 1e-8)
  # With a random slope the reliability is genuinely not constant across waves.
  expect_gt(stats::sd(rbo), 0)
})

test_that("results are reproducible under a fixed seed and vary without one", {
  set.seed(113)
  a <- simulate_longitudinal_polynomial(
    n = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1, timing_sd = 0.1
  )
  set.seed(113)
  b <- simulate_longitudinal_polynomial(
    n = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1, timing_sd = 0.1
  )
  expect_identical(a$y, b$y)
  expect_identical(a$time, b$time)
})

test_that("invalid arguments error informatively", {
  # M < P + 1
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = c(0, 1), fixed_coefficients = c(1, 1, 1),
      error_variance = 1),
    "at least 3 measurement occasions"
  )
  # neither error_variance nor reliability
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(4, 0.25)),
    "exactly one of"
  )
  # both error_variance and reliability
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(4, 0.25), error_variance = 1, reliability = 0.8),
    "exactly one of"
  )
  # reliability with no between-subject variance
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = 0, reliability = 0.8),
    "positive between-subject variance"
  )
  # unequal-length coefficient vectors across groups
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4,
      fixed_coefficients = list(c(10, 1), c(10, 1, 0.5)), error_variance = 1),
    "same length"
  )
  # random_variances of the wrong length
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(4, 0.25, 0.1), error_variance = 1),
    "length 1 or 2"
  )
  # bad n
  expect_error(
    simulate_longitudinal_polynomial(
      n = 2.5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      error_variance = 1),
    "positive integer"
  )
})

test_that("an unachievable target reliability is rejected with a clear message", {
  # Only the slope is random, and one target occasion sits at t = 0, where the
  # true-score variance is 0; the average reliability then cannot reach 1.
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(0, 0.25), reliability = 0.95),
    "not achievable"
  )
})

test_that("random_correlation is validated and feeds the level-two covariance", {
  # A valid correlation produces a non-diagonal T with the right off-diagonal.
  set.seed(113)
  R <- matrix(c(1, -0.5, -0.5, 1), 2, 2)
  d <- simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), random_correlation = R, error_variance = 1
  )
  Tcov <- attr(d, "random_covariance")
  expect_equal(Tcov[1, 2], -0.5 * sqrt(4) * sqrt(0.25))

  # Non-symmetric / wrong-size / non-PSD correlation matrices are rejected.
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(4, 0.25),
      random_correlation = matrix(c(1, 0.2, 0.3, 1), 2, 2), error_variance = 1),
    "symmetric"
  )
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
      random_variances = c(4, 0.25),
      random_correlation = diag(3), error_variance = 1),
    "2-by-2 matrix"
  )
})

test_that("the simulated data compose with plot_trajectories column arguments", {
  skip_if_not_installed("ggplot2")
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 6, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1
  )
  p <- plot_trajectories(d, id = "id", time = "time", outcome = "y",
                         group = "population")
  expect_s3_class(p, "ggplot")
})

test_that("the independent default still reports a diagonal error covariance", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 10, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 2
  )
  Sig <- attr(d, "error_covariance")
  expect_equal(dim(Sig), c(5L, 5L))
  expect_equal(Sig, diag(2, 5))                       # sigma^2 = 2 on the diagonal
  expect_equal(attr(d, "error_variance"), 2)          # scalar, backward compatible
})

test_that("an AR(1) error structure builds the right covariance and leaves reliability alone", {
  set.seed(113)
  indep <- simulate_longitudinal_polynomial(
    n = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1
  )
  set.seed(113)
  ar <- simulate_longitudinal_polynomial(
    n = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1,
    error_structure = "ar1", error_correlation = 0.5
  )
  Sig <- attr(ar, "error_covariance")
  expect_equal(Sig, 0.5 ^ abs(outer(0:4, 0:4, "-")))  # sigma^2 = 1, rho^|j-k|
  expect_equal(diag(Sig), rep(1, 5))                  # marginal variances unchanged
  # Classical reliability depends only on the diagonal, so the AR(1) correlation
  # leaves the per-occasion reliabilities identical to the independent design.
  expect_equal(attr(ar, "reliability_by_occasion"),
               attr(indep, "reliability_by_occasion"))
  # The per-occasion error variance is still 1, but the table is now a vector.
  expect_equal(unname(attr(ar, "error_variance")), rep(1, 5))
})

test_that("compound symmetry and Toeplitz structures build the documented matrices", {
  cs <- simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 2,
    error_structure = "compound_symmetry", error_correlation = 0.3
  )
  Scs <- attr(cs, "error_covariance")
  expect_equal(diag(Scs), rep(2, 5))
  expect_true(all(abs(Scs[upper.tri(Scs)] - 0.3 * 2) < 1e-10))   # rho * sigma^2

  toe <- simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1,
    error_structure = "toeplitz", error_correlation = c(0.5, 0.25, 0.125, 0.0625)
  )
  expect_equal(attr(toe, "error_covariance")[1, ],
               c(1, 0.5, 0.25, 0.125, 0.0625))
})

test_that("a full covariance matrix is used verbatim as the error covariance", {
  Sig <- 0.6 ^ abs(outer(0:4, 0:4, "-"))             # an arbitrary PSD covariance
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 6, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = Sig
  )
  expect_equal(attr(d, "error_covariance"), Sig)
  expect_equal(unname(attr(d, "error_variance")), diag(Sig))
})

test_that("a heteroscedastic error-variance vector sets the per-occasion variances", {
  v <- c(0.5, 1, 1.5, 2, 2.5)
  d <- simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = v
  )
  expect_equal(unname(attr(d, "error_variance")), v)
  expect_equal(diag(attr(d, "error_covariance")), v)
})

test_that("the reliability route still hits its target when an AR(1) structure is added", {
  target <- 0.75
  d <- simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), reliability = target,
    error_structure = "ar1", error_correlation = 0.4
  )
  expect_equal(mean(attr(d, "reliability_by_occasion")), target, tolerance = 1e-6)
  Sig <- attr(d, "error_covariance")
  # The solved error variance is homoscedastic; the AR(1) correlation rides on
  # top (off-diagonal / diagonal = rho).
  expect_true(max(abs(diag(Sig) - diag(Sig)[1])) < 1e-9)
  expect_equal(Sig[1, 2] / Sig[1, 1], 0.4)
})

test_that("results are reproducible under a fixed seed with a correlated error structure", {
  set.seed(113)
  a <- simulate_longitudinal_polynomial(
    n = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1,
    error_structure = "ar1", error_correlation = 0.5
  )
  set.seed(113)
  b <- simulate_longitudinal_polynomial(
    n = 8, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), error_variance = 1,
    error_structure = "ar1", error_correlation = 0.5
  )
  expect_identical(a$y, b$y)
})

test_that("invalid error-structure arguments are rejected with clear messages", {
  base <- function(...) simulate_longitudinal_polynomial(
    n = 5, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = c(4, 0.25), ...)
  # A named structure needs its correlation parameter.
  expect_error(base(error_variance = 1, error_structure = "ar1"),
               "requires 'error_correlation'")
  # AR(1) correlation out of range.
  expect_error(base(error_variance = 1, error_structure = "ar1",
                    error_correlation = 1.2), "\\(-1, 1\\)")
  # Toeplitz correlation of the wrong length.
  expect_error(base(error_variance = 1, error_structure = "toeplitz",
                    error_correlation = c(0.4, 0.2)), "length 4")
  # A correlation supplied for the independent default.
  expect_error(base(error_variance = 1, error_correlation = 0.3),
               "only used when")
  # error_variance vector of the wrong length.
  expect_error(base(error_variance = c(1, 2, 3)), "per-occasion variances")
  # A full matrix together with a structure double-specifies the errors.
  expect_error(base(error_variance = diag(5), error_structure = "ar1",
                    error_correlation = 0.5), "already specifies")
  # A non-symmetric / non-PSD error covariance matrix.
  bad <- diag(5); bad[1, 2] <- 0.4
  expect_error(base(error_variance = bad), "symmetric")
  npsd <- matrix(0.9, 5, 5); diag(npsd) <- 1; npsd[1, 5] <- npsd[5, 1] <- -0.9
  expect_error(base(error_variance = npsd), "positive semidefinite")
})

test_that("an AR(1) error structure induces the intended autocorrelation (Monte Carlo)", {
  skip_on_cran()
  set.seed(113)
  rho <- 0.6
  # No random effects and no jitter, so the true score is a fixed polynomial and
  # y minus true_score is exactly the level-one error; its lag-k correlation
  # should track rho^k.
  d <- simulate_longitudinal_polynomial(
    n = 8000, target_times = 0:4, fixed_coefficients = c(10, 1.5),
    random_variances = 0, error_variance = 1,
    error_structure = "ar1", error_correlation = rho
  )
  e <- matrix(d$y - d$true_score, ncol = 5, byrow = TRUE)
  Rhat <- cor(e)
  # Absolute deviations (the sampling SE of each correlation is about 0.007).
  expect_lt(abs(Rhat[1, 2] - rho),   0.04)   # lag 1
  expect_lt(abs(Rhat[1, 3] - rho^2), 0.04)   # lag 2
  expect_lt(abs(Rhat[1, 4] - rho^3), 0.04)   # lag 3
})

test_that("unit-specific schedules draw per-unit times inside the range", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 15, time_range = c(40, 90), occasions = c(5, 9),
    fixed_coefficients = c(10, 0.5), random_variances = c(4, 0.01),
    error_variance = 2
  )
  expect_identical(attr(d, "schedule"), "unit_specific")
  expect_true(all(d$time >= 40 & d$time <= 90))
  counts <- as.integer(table(d$id)[table(d$id) > 0])
  expect_true(all(counts >= 5 & counts <= 9))
  expect_true(all(tapply(d$time, d$id, function(x) !is.unsorted(x)),
                  na.rm = TRUE))
  expect_identical(d$target_time, d$time)
  expect_true(is.na(attr(d, "reliability_by_occasion")))
  expect_equal(attr(d, "error_variance"), 2)
})

test_that("a scalar 'occasions' gives every unit the same count", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 8, time_range = c(0, 10), occasions = 4,
    fixed_coefficients = c(5, 1), error_variance = 1
  )
  expect_true(all(table(d$id)[table(d$id) > 0] == 4))
})

test_that("the unit-specific branch is deterministic under a seed", {
  set.seed(7)
  a <- simulate_longitudinal_polynomial(
    n = 6, time_range = c(0, 5), occasions = c(2, 4),
    fixed_coefficients = c(2, 1), error_variance = 1)
  set.seed(7)
  b <- simulate_longitudinal_polynomial(
    n = 6, time_range = c(0, 5), occasions = c(2, 4),
    fixed_coefficients = c(2, 1), error_variance = 1)
  expect_identical(a, b)
})

test_that("the schedule branch refuses what needs a shared grid", {
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, time_range = c(0, 10), occasions = 4,
      fixed_coefficients = c(1, 1), reliability = 0.8,
      random_variances = 1),
    "needs a shared measurement schedule")
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, time_range = c(0, 10), occasions = 4,
      fixed_coefficients = c(1, 1), error_variance = 1,
      error_structure = "ar1", error_correlation = 0.5),
    "must be \"independent\"")
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, time_range = c(0, 10), occasions = 4,
      fixed_coefficients = c(1, 1), error_variance = 1,
      timing_sd = 0.2),
    "does not combine with 'time_range'")
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, time_range = c(0, 10), occasions = 4,
      fixed_coefficients = c(1, 1),
      error_variance = matrix(1, 4, 4)),
    "single")
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, fixed_coefficients = c(1, 1), error_variance = 1),
    "exactly one of 'target_times'")
})

test_that("every unit keeps enough occasions to identify the polynomial", {
  # A degree-2 model with occasions allowing as few as 2 must stop.
  expect_error(
    simulate_longitudinal_polynomial(
      n = 5, time_range = c(0, 10), occasions = c(2, 5),
      fixed_coefficients = c(1, 1, 1), error_variance = 1),
    "at least 3 measurement occasions per unit")
})

test_that("analysis_of_change() recovers a polynomial from unit-specific times", {
  set.seed(113)
  d <- simulate_longitudinal_polynomial(
    n = 60, time_range = c(0, 10), occasions = c(4, 7),
    fixed_coefficients = c(b0 = 10, b1 = 1.5),
    random_variances = c(1, 0.04), error_variance = 0.25
  )
  fit <- analysis_of_change(d, id = "id", time = "time", outcome = "y",
                            model = "polynomial", order = 1)
  expect_lt(abs(fit$estimate[fit$term == "b0"] - 10), 0.6)
  expect_lt(abs(fit$estimate[fit$term == "b1"] - 1.5), 0.15)
})
