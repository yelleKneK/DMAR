test_that("power_equivalence_c() reduces to power_equivalence_md() for two
           equal groups", {
  # Phillips (1990) Table 1, 5th row, 5th column: expected 0.8029678.
  ref <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              ltheta1 = -.2, ltheta2 = .2, ldiff = .05,
                              sigma = .20, n = 24, nu = 22)
  res <- power_equivalence_c(c_weights = c(1, -1), n = 24, sigma = .20,
                             delta_lower = .2, delta_upper = .2,
                             true_psi = .05, df_error = 22)
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
                              ltheta1 = -5, ltheta2 = 5, ldiff = 0,
                              sigma = 15.673, n = n_eff, nu = 399)
  res <- power_equivalence_c(c_weights = c(1, -1), n = c(61, 113),
                             sigma = 15.673, delta_upper = 5,
                             true_psi = 0, df_error = 399)
  expect_equal(res$value, ref$value, tolerance = 1e-12)
  expect_equal(res$value, 0.281, tolerance = 1e-3)
})

test_that("power_equivalence_c() noninferiority power matches the closed
           form", {
  k  <- sqrt(sum(c(1, -1)^2 / c(61, 113)))
  nu <- 399
  ncp <- (0 + 5) / (15.673 * k)
  expected <- pt(qt(.95, nu), nu, ncp = ncp, lower.tail = FALSE)
  res <- power_equivalence_c(c_weights = c(1, -1), n = c(61, 113),
                             sigma = 15.673, delta_upper = 5,
                             true_psi = 0, df_error = 399,
                             side = "noninferiority")
  expect_equal(res$value, expected, tolerance = 1e-12)
})

test_that("power_equivalence_c() equivalence power is zero for an infeasible
           design", {
  # The CI can never fit inside the bounds when sigma dwarfs them.
  res <- power_equivalence_c(c_weights = c(1, -1), n = 5, sigma = 100,
                             delta_upper = 1)
  expect_equal(res$value, 0)
})

test_that("power_equivalence_c() power increases with n and decreases as
           true_psi approaches a bound", {
  p_small <- power_equivalence_c(c_weights = c(1, -1), n = 50,
                                 sigma = 10, delta_upper = 5)$value
  p_large <- power_equivalence_c(c_weights = c(1, -1), n = 200,
                                 sigma = 10, delta_upper = 5)$value
  expect_gt(p_large, p_small)

  p_center <- power_equivalence_c(c_weights = c(1, -1), n = 200,
                                  sigma = 10, delta_upper = 5,
                                  true_psi = 0)$value
  p_edge   <- power_equivalence_c(c_weights = c(1, -1), n = 200,
                                  sigma = 10, delta_upper = 5,
                                  true_psi = 3)$value
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
  res <- power_equivalence_c(c_weights = cw, n = n, sigma = sigma,
                             delta_upper = 5, side = "noninferiority")
  expect_equal(res$value, expected, tolerance = 1e-12)

  # Equivalence power on the same contrast agrees with a Monte Carlo
  # declaration rate from equivalence_c().
  skip_on_cran()
  set.seed(113)
  analytic <- power_equivalence_c(c_weights = cw, n = n, sigma = sigma,
                                  delta_upper = 5)$value
  G <- 4000
  decl <- replicate(G, {
    y     <- lapply(n, function(nj) rnorm(nj, mean = 0, sd = sigma))
    means <- vapply(y, mean, numeric(1))
    mse   <- sum(mapply(function(yj, mj) sum((yj - mj)^2), y, means)) / nu
    tt <- equivalence_c(means = means, s_anova = sqrt(mse), c_weights = cw, n = n,
                 delta_upper = 5)
    tt$value[tt$term == "equivalent"]
  })
  mc <- mean(decl); se <- sqrt(mc * (1 - mc) / G)
  expect_lt(abs(analytic - mc), 4 * se)
})

test_that("power_equivalence_c() validates its inputs", {
  expect_error(power_equivalence_c(c_weights = c(2, -2), n = 20, sigma = 1,
                                   delta_upper = 1),
               "positive weights must sum to 1")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n = 20, sigma = 1),
               "delta_upper")
  expect_error(power_equivalence_c(c_weights = c(1, -1), n = c(20, 20, 20),
                                   sigma = 1, delta_upper = 1),
               "lengths")
})

test_that("power_equivalence_c() is invariant to the measurement scale", {
  # CRITICAL-01 regression. The equivalence power is a probability and must
  # not change when every response-scale input (sigma, bounds, true_psi) is
  # multiplied by a common factor: dollars versus millions of dollars, or
  # meters versus micrometers, describe the same study. A legacy absolute
  # '- 1e-6' nudge on the integration endpoint broke this, driving the power
  # to 0 on a small raw scale.
  power_at_scale <- function(s) {
    power_equivalence_c(c_weights = c(1, -1), n = 24, sigma = 0.20 * s,
                        delta_lower = 0.2 * s, delta_upper = 0.2 * s,
                        true_psi = 0.05 * s, df_error = 22)$value
  }
  base <- power_at_scale(1)
  for (s in c(1e-6, 1e-3, 1e3, 1e6, 1e9)) {
    expect_equal(power_at_scale(s), base, tolerance = 1e-10)
  }
  # And the invariant value is still Phillips's (1990) published power.
  expect_equal(base, 0.8029678, tolerance = 1e-6)
})
