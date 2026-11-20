# Regression tests for the quadrature behind power_equivalence_md() and
# power_equivalence_c().
#
# The defect these guard against: the power was obtained by integrating the
# chi density of the estimated error standard deviation over the raw interval
# (0, upper], where 'upper' is the standard deviation at which the confidence
# interval just fits inside the equivalence interval. For a precise design
# 'upper' sits far out in the right tail, so that interval is many times wider
# than the region where the error standard deviation has any mass. The initial
# Gauss-Kronrod nodes then all landed where the density had underflowed and
# QUADPACK returned message "OK", with a negligible error estimate, on an
# essentially zero integral. Designs whose true power is 1 were reported as
# 1.7e-17 and 9.6e-14, and the reported power was not even monotone in sigma.

# An oracle that is independent of the package's quadrature: the probability
# integral transform. Writing r for the error standard deviation divided by
# sigma, the power is E[phi(r) 1{r <= r_max}], which the substitution
# r = sqrt(qchisq(p, nu) / nu) turns into an average of phi over p uniform on
# (0, p_star), p_star = P(r <= r_max). A midpoint rule on that average uses no
# adaptive quadrature at all and agrees with the closed-form integration to
# about eight decimal places.
pit_midpoint_power <- function(alpha, theta1, theta2, diff, sigma, k, nu,
                               M = 50000L) {
  crit_t      <- stats::qt(1 - alpha, df = nu)
  one_over_se <- 1 / (sigma * k)
  z_lower     <- (theta1 - diff) * one_over_se
  z_upper     <- (theta2 - diff) * one_over_se
  r_max       <- (z_upper - z_lower) / (2 * crit_t)
  p_star      <- stats::pchisq(nu * r_max^2, df = nu)
  p <- (seq_len(M) - 0.5) / M * p_star
  r <- sqrt(stats::qchisq(p, df = nu) / nu)
  p_star * mean(pmax(stats::pnorm(z_upper - r * crit_t) -
                       stats::pnorm(z_lower + r * crit_t), 0))
}

# Monte Carlo oracle: draw the estimated contrast and the estimated error
# standard deviation, and count how often the interval lands inside the
# equivalence interval.
mc_equivalence_power <- function(alpha, theta1, theta2, diff, sigma, k, nu,
                                 G = 2e5) {
  se  <- sigma * k
  est <- stats::rnorm(G, mean = diff, sd = se)
  hw  <- stats::qt(1 - alpha, df = nu) * se *
    sqrt(stats::rchisq(G, df = nu) / nu)
  p <- mean((est - hw > theta1) & (est + hw < theta2))
  c(power = p, se = sqrt(p * (1 - p) / G))
}


test_that("power_equivalence_md() reports the true power for a precise design", {
  # Both designs returned essentially zero before the quadrature was fixed
  # (1.71e-17 and 9.59e-14). A Monte Carlo oracle with G = 1e6 declares
  # equivalence in every replication for each of them, that is, an estimated
  # power of 1.000000 with a standard error of 0.
  bio <- power_equivalence_md(alpha_level = .05, logscale = TRUE,
                              lower_bound = .8, upper_bound = 1.25, true_effect = 1.0,
                              true_sd = .02, n_per_treatment = 500, df_error = 400)
  expect_equal(bio$value, 1, tolerance = 1e-9)

  raw <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -.2, upper_bound = .2, true_effect = 0,
                              true_sd = 0.05, n_per_treatment = 500, df_error = 400)
  expect_equal(raw$value, 1, tolerance = 1e-9)
})


test_that("power_equivalence_md() matches a Monte Carlo oracle where the old
           quadrature collapsed", {
  skip_on_cran()
  # Three designs whose true power is far from both 0 and 1 and for which the
  # old integration returned 8.8e-06, 3.6e-14, and 3.1e-05. The expected
  # values are Monte Carlo estimates with G = 2e7 (standard errors given), so
  # they are independent of any quadrature.
  designs <- list(
    list(diff = .196, sigma = .02, n = 50,  nu = 400,  mc = 0.2588978, se = 9.79e-05),
    list(diff = .196, sigma = .02, n = 500, nu = 200,  mc = 0.9339972, se = 5.55e-05),
    list(diff = .19,  sigma = .03, n = 50,  nu = 1000, mc = 0.5081041, se = 1.12e-04)
  )
  for (d in designs) {
    got <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                lower_bound = -.2, upper_bound = .2, true_effect = d$diff,
                                true_sd = d$sigma, n_per_treatment = d$n, df_error = d$nu)$value
    expect_lt(abs(got - d$mc), 4 * d$se)
    # And against the transform-based oracle, which is exact to about 1e-8.
    expect_equal(got,
                 pit_midpoint_power(.05, -.2, .2, d$diff, d$sigma,
                                    sqrt(2 / d$n), d$nu),
                 tolerance = 1e-5)
  }
})


test_that("power_equivalence_md() agrees with a live Monte Carlo run", {
  skip_on_cran()
  set.seed(113)
  d <- mc_equivalence_power(.05, -.2, .2, diff = .196, sigma = .02,
                            k = sqrt(2 / 50), nu = 400, G = 2e5)
  got <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                              lower_bound = -.2, upper_bound = .2, true_effect = .196,
                              true_sd = .02, n_per_treatment = 50, df_error = 400)$value
  expect_lt(abs(got - d[["power"]]), 4 * d[["se"]])
})


test_that("power_equivalence_md() power is monotone in sigma", {
  # A more precise study cannot be less likely to declare equivalence, so the
  # power must be non-increasing in sigma. The old quadrature returned
  # 9.6e-14, 9.4e-07, and 3.4e-28 at sigma = 0.05, 0.01, and 0.02, a sequence
  # that is not monotone in either direction and so could not be explained by
  # underflow.
  sigmas <- c(1e-6, 1e-4, .001, .002, .005, .01, .02, .05, .1, .2, .3, .4, .6,
              1, 1.5, 2, 3, 5)
  pw <- vapply(sigmas, function(s)
    power_equivalence_md(alpha_level = .05, logscale = FALSE,
                         lower_bound = -.2, upper_bound = .2, true_effect = 0,
                         true_sd = s, n_per_treatment = 500, df_error = 400)$value, numeric(1))
  expect_true(all(diff(pw) <= 0))
  expect_equal(pw[1], 1, tolerance = 1e-12)
  expect_lt(pw[length(pw)], 1e-3)

  # Same check for the contrast form, on a three-group fractional contrast.
  pw_c <- vapply(sigmas, function(s)
    power_equivalence_c(c_weights = c(0.5, 0.5, -1), n_per_group = c(200, 400, 300),
                        true_sd = s, contrast_upper = .2, true_contrast = 0,
                        df_error = 897)$value, numeric(1))
  expect_true(all(diff(pw_c) <= 0))
})


test_that("power_equivalence_md() and power_equivalence_c() return a
           probability", {
  skip_on_cran()
  # The quadrature could overshoot: n = 100, nu = 98, sigma = 0.03 previously
  # returned 1.000000000034, and the return is documented as a power in
  # [0, 1].
  grid <- expand.grid(n = c(20, 60, 100, 200, 500), nu = c(10, 58, 98, 198, 400),
                      sigma = c(.001, .01, .02, .03, .05, .2, 1))
  pw <- vapply(seq_len(nrow(grid)), function(i)
    power_equivalence_md(alpha_level = .05, logscale = FALSE,
                         lower_bound = -.2, upper_bound = .2, true_effect = 0,
                         true_sd = grid$sigma[i], n_per_treatment = grid$n[i],
                         df_error = grid$nu[i])$value, numeric(1))
  expect_true(all(pw >= 0 & pw <= 1))
  expect_true(any(pw == 1))            # the clamp is actually exercised

  pw_c <- vapply(seq_len(nrow(grid)), function(i)
    power_equivalence_c(c_weights = c(1, -1), n_per_group = grid$n[i],
                        true_sd = grid$sigma[i], contrast_upper = .2,
                        df_error = grid$nu[i])$value, numeric(1))
  expect_true(all(pw_c >= 0 & pw_c <= 1))
})


test_that("the published equivalence powers are unchanged by the new
           quadrature", {
  # Phillips (1990, Table 1, row 5, column 5) and Diletti (1991, Table 1),
  # plus a low-power corner, all of which the old integration handled
  # correctly. The values below are the pre-fix returns; agreement to 1e-12
  # shows the new integration did not move any design that was already right.
  expect_equal(power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = -.2, upper_bound = .2, true_effect = .05,
                                    true_sd = .20, n_per_treatment = 24, df_error = 22)$value,
               0.80296783413873, tolerance = 1e-12)
  expect_equal(power_equivalence_md(alpha_level = .05, logscale = TRUE,
                                    lower_bound = .8, upper_bound = 1.25, true_effect = 1.05,
                                    true_sd = .20, n_per_treatment = 18, df_error = 16)$value,
               0.79227955226626, tolerance = 1e-12)
  expect_equal(power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                    true_sd = 1, n_per_treatment = 6, df_error = 4)$value,
               7.4953108424477e-05, tolerance = 1e-10)

  # The transform-based oracle reaches the same Phillips value.
  expect_equal(pit_midpoint_power(.05, -.2, .2, .05, .20, sqrt(2 / 24), 22),
               0.8029678, tolerance = 1e-6)
})


test_that("equivalence power is invariant to the measurement scale, including
           where the old quadrature collapsed", {
  skip_on_cran()
  # Dollars or millions of dollars describe the same study. The check is run
  # at a design whose power the old code reported as 9.6e-14.
  power_at_scale <- function(s) {
    power_equivalence_md(alpha_level = .05, logscale = FALSE,
                         lower_bound = -0.2 * s, upper_bound = 0.2 * s, true_effect = 0,
                         true_sd = 0.05 * s, n_per_treatment = 500, df_error = 400)$value
  }
  base <- power_at_scale(1)
  expect_equal(base, 1, tolerance = 1e-9)
  for (s in c(1e-9, 1e-6, 1e-3, 1e3, 1e6, 1e9)) {
    expect_equal(power_at_scale(s), base, tolerance = 1e-10)
  }

  # And at Phillips's design, on both interfaces.
  md_at_scale <- function(s)
    power_equivalence_md(alpha_level = .05, logscale = FALSE, lower_bound = -0.2 * s,
                         upper_bound = 0.2 * s, true_effect = 0.05 * s,
                         true_sd = 0.20 * s, n_per_treatment = 24, df_error = 22)$value
  c_at_scale <- function(s)
    power_equivalence_c(c_weights = c(1, -1), n_per_group = 24, true_sd = 0.20 * s,
                        contrast_lower = 0.2 * s, contrast_upper = 0.2 * s,
                        true_contrast = 0.05 * s, df_error = 22)$value
  for (s in c(1e-6, 1e-3, 1e3, 1e6, 1e9)) {
    expect_equal(md_at_scale(s), 0.8029678, tolerance = 1e-6)
    expect_equal(c_at_scale(s), 0.8029678, tolerance = 1e-6)
  }
})


test_that("power_equivalence_c() reports the true power for a precise
           contrast design", {
  skip_on_cran()
  # Two groups, the mean-difference design that used to return 9.59e-14.
  two_group <- power_equivalence_c(c_weights = c(1, -1), n_per_group = 500, true_sd = 0.05,
                                   contrast_upper = 0.2, true_contrast = 0,
                                   df_error = 400)$value
  expect_equal(two_group, 1, tolerance = 1e-9)
  expect_equal(two_group,
               power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                    true_sd = 0.05, n_per_treatment = 500, df_error = 400)$value,
               tolerance = 1e-14)

  # A three-group fractional contrast with unequal group sizes, evaluated
  # close enough to the bound that the power is genuinely intermediate. The
  # old integration returned 6.6e-06 here. The expected value is a Monte
  # Carlo estimate with G = 1e7 (standard error 1.38e-04).
  cw <- c(0.5, 0.5, -1); n <- c(200, 400, 300)
  got <- power_equivalence_c(c_weights = cw, n_per_group = n, true_sd = 0.03,
                             contrast_upper = 0.2, true_contrast = 0.195,
                             df_error = 897)$value
  expect_lt(abs(got - 0.7462059), 4 * 1.38e-04)
  expect_equal(got,
               pit_midpoint_power(.05, -.2, .2, .195, .03,
                                  sqrt(sum(cw^2 / n)), 897),
               tolerance = 1e-5)
})


test_that("equivalence power still behaves at the edges of its argument
           space", {
  # Small, fractional, and very large degrees of freedom, and a design so
  # imprecise that equivalence is unreachable. None of these may error, and
  # all must return a probability.
  for (nu in c(0.5, 1, 2, 3, 22, 1e4, 1e6)) {
    p_precise <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                      lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                      true_sd = .02, n_per_treatment = 500, df_error = nu)$value
    p_vague <- power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                    true_sd = 50, n_per_treatment = 500, df_error = nu)$value
    expect_true(is.finite(p_precise) && p_precise >= 0 && p_precise <= 1)
    expect_true(is.finite(p_vague) && p_vague >= 0 && p_vague <= 1)
    expect_gt(p_precise, p_vague)
  }
  expect_equal(power_equivalence_md(alpha_level = .05, logscale = FALSE,
                                    lower_bound = -.2, upper_bound = .2, true_effect = 0,
                                    true_sd = 1e9, n_per_treatment = 24, df_error = 22)$value, 0)
})

test_that("the shared quadrature takes the error degrees of freedom as df_error", {
  skip_on_cran()
  # .power_equivalence_tost_quad() named the degrees of freedom nu through
  # the first 1.1.0 builds, and .power_equivalence_c_fast() took nu and
  # sigma; they take df_error and true_sd, the names of the public
  # arguments they serve. The value agrees with the probability integral
  # transform oracle above.
  expect_identical(names(formals(.power_equivalence_tost_quad)),
                   c("power_t", "z_lower", "z_upper", "df_error"))
  expect_true(all(c("df_error", "true_sd") %in%
                    names(formals(.power_equivalence_c_fast))))
  expect_false(any(c("nu", "sigma") %in%
                     c(names(formals(.power_equivalence_tost_quad)),
                       names(formals(.power_equivalence_c_fast)))))
  k <- sqrt(2 / 24)
  t_crit <- stats::qt(.95, 22)
  got <- .power_equivalence_tost_quad(
    power_t = t_crit, z_lower = (-.2 - .05) / (.2 * k),
    z_upper = (.2 - .05) / (.2 * k), df_error = 22)
  # The midpoint rule of the oracle is good to about 1e-6 here, the
  # tolerance the same design is checked at above.
  expect_equal(got, pit_midpoint_power(.05, -.2, .2, .05, .2, k, 22),
               tolerance = 1e-6)
})
