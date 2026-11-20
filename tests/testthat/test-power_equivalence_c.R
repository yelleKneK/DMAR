test_that("power_equivalence_c() reduces to power_equivalence_md() for two
           equal groups", {
  # Phillips (1990) Table 1, 5th row, 5th column: expected 0.8029678.
  ref <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -.2, upper_bound = .2, true_effect = .05,
                              true_sd = .20, n_per_treatment = 24, df_error = 22)
  res <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 24, true_sd = .20,
                             contrast_lower = .2, contrast_upper = .2,
                             true_contrast = .05, df_error = 22)
  expect_equal(res$value, ref$value, tolerance = 1e-10)
  expect_equal(res$value, 0.8029678, tolerance = 1e-6)
})

test_that("power_equivalence_c() handles unequal group sizes with a shared
           error term", {
  # Two groups of 61 and 113 sharing a five-group pooled error term. The
  # result must agree exactly with power_equivalence_md() evaluated at the
  # effective sample size solving 2 / n_eff = sum(c_j^2 / n_j), and sits
  # near .281 for this design.
  n_eff <- 2 / (1 / 61 + 1 / 113)
  ref <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -5, upper_bound = 5, true_effect = 0,
                              true_sd = 15.673, n_per_treatment = n_eff, df_error = 399)
  res <- power_equivalence_c(c_weights = c(1, -1), n_per_group = c(61, 113),
                             true_sd = 15.673, contrast_upper = 5,
                             true_contrast = 0, df_error = 399)
  expect_equal(res$value, ref$value, tolerance = 1e-12)
  expect_equal(res$value, 0.281, tolerance = 1e-3)
})

test_that("power_equivalence_c() noninferiority power matches the closed
           form", {
  k  <- sqrt(sum(c(1, -1)^2 / c(61, 113)))
  nu <- 399
  ncp <- (0 + 5) / (15.673 * k)
  expected <- pt(qt(.95, nu), nu, ncp = ncp, lower.tail = FALSE)
  res <- power_equivalence_c(c_weights = c(1, -1), n_per_group = c(61, 113),
                             true_sd = 15.673, contrast_upper = 5,
                             true_contrast = 0, df_error = 399,
                             side = "noninferiority")
  expect_equal(res$value, expected, tolerance = 1e-12)
})

test_that("power_equivalence_c() equivalence power is zero for an infeasible
           design", {
  # The CI can never fit inside the bounds when sigma dwarfs them.
  res <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 5, true_sd = 100,
                             contrast_upper = 1)
  expect_equal(res$value, 0)
})

test_that("power_equivalence_c() power increases with n and decreases as
           true_contrast approaches a bound", {
  p_small <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 50,
                                 true_sd = 10, contrast_upper = 5)$value
  p_large <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 200,
                                 true_sd = 10, contrast_upper = 5)$value
  expect_gt(p_large, p_small)

  p_center <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 200,
                                  true_sd = 10, contrast_upper = 5,
                                  true_contrast = 0)$value
  p_edge   <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 200,
                                  true_sd = 10, contrast_upper = 5,
                                  true_contrast = 3)$value
  expect_gt(p_center, p_edge)
})

test_that("power_equivalence_c() handles a genuine three-group fractional
           contrast (the summation power_equivalence_md cannot reach)", {
  cw <- c(0.5, 0.5, -1); n <- c(30, 40, 50); sigma <- 8
  nu <- sum(n) - length(cw)               # 117
  k  <- sqrt(sum(cw^2 / n))               # SE factor summed over three weights

  # Noninferiority reduces to a single noncentral t probability, so the
  # analytic value can be hand-built independently of the function.
  ncp      <- (0 + 5) / (sigma * k)
  expected <- pt(qt(.95, nu), nu, ncp = ncp, lower.tail = FALSE)
  res <- power_equivalence_c(c_weights = cw, n_per_group = n, true_sd = sigma,
                             contrast_upper = 5, side = "noninferiority")
  expect_equal(res$value, expected, tolerance = 1e-12)

  # Equivalence power on the same contrast agrees with a Monte Carlo
  # declaration rate from equivalence_c().
  skip_on_cran()
  set.seed(113)
  analytic <- power_equivalence_c(c_weights = cw, n_per_group = n, true_sd = sigma,
                                  contrast_upper = 5)$value
  G <- 4000
  decl <- replicate(G, {
    y     <- lapply(n, function(nj) rnorm(nj, mean = 0, sd = sigma))
    means <- vapply(y, mean, numeric(1))
    mse   <- sum(mapply(function(yj, mj) sum((yj - mj)^2), y, means)) / nu
    tt <- equivalence_c(means = means, s_anova = sqrt(mse), c_weights = cw, n_per_group = n,
                 contrast_upper = 5)
    tt$value[tt$term == "equivalent"]
  })
  mc <- mean(decl); se <- sqrt(mc * (1 - mc) / G)
  expect_lt(abs(analytic - mc), 4 * se)
})

test_that("power_equivalence_c() validates its inputs", {
  expect_error(power_equivalence_c(c_weights = c(2, -2), n_per_group = 20, true_sd = 1,
                                   contrast_upper = 1),
               "positive weights must sum to 1")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = 20, true_sd = 1),
               "contrast_upper")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = c(20, 20, 20),
                                   true_sd = 1, contrast_upper = 1),
               "lengths")
})

test_that("power_equivalence_c() is invariant to the measurement scale", {
  # CRITICAL-01 regression. The equivalence power is a probability and must
  # not change when every response-scale input (sigma, bounds, true_contrast) is
  # multiplied by a common factor: dollars versus millions of dollars, or
  # meters versus micrometers, describe the same study. A legacy absolute
  # '- 1e-6' nudge on the integration endpoint broke this, driving the power
  # to 0 on a small raw scale.
  power_at_scale <- function(s) {
    power_equivalence_c(c_weights = c(1, -1), n_per_group = 24, true_sd = 0.20 * s,
                        contrast_lower = 0.2 * s, contrast_upper = 0.2 * s,
                        true_contrast = 0.05 * s, df_error = 22)$value
  }
  base <- power_at_scale(1)
  for (s in c(1e-6, 1e-3, 1e3, 1e6, 1e9)) {
    expect_equal(power_at_scale(s), base, tolerance = 1e-10)
  }
  # And the invariant value is still Phillips's (1990) published power.
  expect_equal(base, 0.8029678, tolerance = 1e-6)
})

test_that("power_equivalence_c() computes noninferiority power without contrast_upper", {
  # DMAR 1.0.x stopped for want of the upper bound even when only -contrast_lower
  # is tested. The noninferiority power never reads contrast_upper, so the
  # margin alone now suffices, and a call that gives contrast_upper alone, as
  # 1.0.x required, still takes it as the margin.
  args <- list(c_weights = c(1, -1), n_per_group = c(61, 113), true_sd = 15.673,
               side = "noninferiority", df_error = 399)
  lower_only <- do.call(power_equivalence_c, c(args, contrast_lower = 5))
  upper_only <- do.call(power_equivalence_c, c(args, contrast_upper = 5))
  both       <- do.call(power_equivalence_c,
                        c(args, contrast_lower = 5, contrast_upper = 50))
  k <- sqrt(1 / 61 + 1 / 113)
  expected <- pt(qt(.95, 399), 399, ncp = 5 / (15.673 * k),
                 lower.tail = FALSE)
  expect_equal(lower_only$value, expected, tolerance = 1e-12)
  expect_identical(upper_only$value, lower_only$value)
  expect_identical(both$value, lower_only$value)

  # With neither bound the message names the margin; equivalence still
  # needs contrast_upper, and a bound that is given is checked even where the
  # test does not use it.
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = 20, true_sd = 1,
                                   side = "noninferiority"),
               "'contrast_lower' must be specified")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = 20, true_sd = 1,
                                   contrast_lower = 1),
               "'contrast_upper' must be specified")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = 20, true_sd = 1,
                                   contrast_lower = 1, contrast_upper = -1,
                                   side = "noninferiority"),
               "'contrast_upper' must be a single positive number")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n_per_group = 20, true_sd = 1,
                                   contrast_lower = NA_real_,
                                   side = "noninferiority"),
               "'contrast_lower' must be a single positive number")
})

test_that("power_equivalence_c() paired noninferiority power reproduces its anchor two ways", {
  # 50 units, difference scores with SD 1, a margin of .5, and a true
  # contrast of 0: the standard error is 1 / sqrt(50) on 49 df, and the
  # power is 1 - pt(qt(.95, 49), 49, ncp = .5 * sqrt(50)) = .9672067. The
  # same number comes from a route that never calls the noncentral t: the
  # test rejects when Z + lambda > t_crit * sqrt(V / 49), with Z standard
  # normal and V chi square on 49 df, so the power is the integral of
  # pnorm(lambda - t_crit * sqrt(v / 49)) against the chi square density.
  res <- power_equivalence_c(N = 50, true_sd = 1, contrast_lower = .5,
                             side = "noninferiority", design = "paired")
  expect_identical(res$term, "power")
  lambda <- .5 * sqrt(50)
  t_crit <- qt(.95, 49)
  by_integral <- integrate(function(v)
                             pnorm(lambda - t_crit * sqrt(v / 49)) *
                               dchisq(v, 49),
                           0, Inf, rel.tol = 1e-12)$value
  expect_equal(res$value, by_integral, tolerance = 1e-9)
  expect_equal(res$value, 1 - pt(qt(.95, 49), 49, ncp = .5 * sqrt(50)),
               tolerance = 1e-12)
  expect_equal(res$value, 0.9672067, tolerance = 1e-7)

  # A supplied df_error replaces n - 1, as it replaces N - J for groups.
  r40 <- power_equivalence_c(N = 50, true_sd = 1, contrast_lower = .5,
                             side = "noninferiority", design = "paired",
                             df_error = 40)
  expect_equal(r40$value,
               pt(qt(.95, 40), 40, ncp = .5 * sqrt(50), lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("power_equivalence_c() paired equivalence power matches the two-group form and a simulation", {
  # The mean of n difference scores with SD sigma has standard error
  # sigma / sqrt(n), which is power_equivalence_md()'s sigma * sqrt(2 / m)
  # at m = 2n; with nu = n - 1 the two must agree.
  n <- 30
  res <- power_equivalence_c(N = n, true_sd = 1, contrast_upper = .5,
                             true_contrast = .1, design = "paired")
  ref <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -.5, upper_bound = .5, true_effect = .1,
                              true_sd = 1, n_per_treatment = 2 * n, df_error = n - 1)
  expect_equal(res$value, ref$value, tolerance = 1e-12)

  # And the declaration rate of the TOST on simulated difference scores.
  skip_on_cran()
  set.seed(113)
  G <- 20000
  d <- matrix(rnorm(G * n, mean = .1, sd = 1), nrow = G)
  m <- rowMeans(d)
  s <- sqrt(rowSums((d - m)^2) / (n - 1))
  h <- qt(.95, n - 1) * s / sqrt(n)
  mc <- mean(m - h > -.5 & m + h < .5)
  expect_lt(abs(res$value - mc), 4 * sqrt(mc * (1 - mc) / G))
})

test_that("power_equivalence_c() validates the design", {
  expect_error(power_equivalence_c(n_per_group = 20, true_sd = 1, contrast_upper = 1),
               "'c_weights' must be supplied")
  expect_error(power_equivalence_c(c_weights = c(1, -1), N = 20, true_sd = 1,
                                   contrast_upper = 1, design = "paired"),
               "Leave 'c_weights' NULL")
  expect_error(power_equivalence_c(N = c(20, 20), true_sd = 1, contrast_upper = 1,
                                   design = "paired"),
               "number of units")
  expect_error(power_equivalence_c(N = 1, true_sd = 1, contrast_upper = 1,
                                   design = "paired"),
               "at least 2")
  expect_error(power_equivalence_c(n_per_group = 20, true_sd = 1, contrast_upper = 1,
                                   design = "crossover"))
  # Naming the default design changes nothing.
  args <- list(c_weights = c(1, -1), n_per_group = 24, true_sd = .2, contrast_upper = .2,
               true_contrast = .05, df_error = 22)
  expect_identical(do.call(power_equivalence_c,
                           c(args, design = "independent"))$value,
                   do.call(power_equivalence_c, args)$value)
})

test_that("power_equivalence_c() noninferiority power is exact past a noncentrality of 37.62", {
  # Two per group, contrast_lower 40, true_sd 1: ncp = 40 on 2 df at alpha
  # .001. pt() (a normal approximation there) gave .9660579348 in DMAR
  # through 3338679.
  res <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 2, true_sd = 1,
                             contrast_lower = 40, true_contrast = 0,
                             side = "noninferiority", alpha_level = .001)
  expect_equal(res$value, 0.9591889337, tolerance = 1e-9)
  expect_equal(res$value, .oracle_nct(qt(.999, 2), 2, 40, lower_tail = FALSE),
               tolerance = 1e-12)
})
