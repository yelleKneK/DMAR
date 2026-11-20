# The quantile, power, and density routines of R/nc_t_internals.R (the
# distribution function and the switch rule are tested in test-ci_nc_t.R).
# The oracles are in helper-nc_t_oracle.R. The evidence for the values and
# rules is in work/round6/evidence/e2_nct_users.

test_that(".nc_t_quantile() is exact where qt() is a normal approximation", {
  q <- DMAR:::.nc_t_quantile
  # The review's cases: qt() gave 429.17833 and 83.93688 (and 29.07563 and
  # 31.74404 for the .025 quantiles).
  cases <- list(list(.975, 3, 50, 186.4938651548), list(.975, 10, 45, 79.0839051910),
                list(.025, 3, 50, 28.2828377277), list(.025, 10, 45, 31.3532061443))
  for (x in cases) {
    value <- expect_no_warning(q(x[[1]], x[[2]], x[[3]]))
    expect_equal(value, x[[4]], tolerance = 1e-10)
    tail <- if (x[[1]] > .5) .oracle_nct(value, x[[2]], x[[3]], lower_tail = FALSE) else
      .oracle_nct(value, x[[2]], x[[3]])
    expect_equal(tail, min(x[[1]], 1 - x[[1]]), tolerance = 1e-10)
  }
  expect_equal(q(.975, 3, 50),
               .oracle_nct_quantile(.975, 3, 50, 150, 250), tolerance = 1e-10)
  # With one degree of freedom and ncp = 38, qt() returns Inf.
  expect_identical(suppressWarnings(stats::qt(.975, 1, 38)), Inf)
  value <- q(.975, 1, 38)
  expect_equal(.oracle_nct(value, 1, 38, lower_tail = FALSE), .025, tolerance = 1e-10)
})

test_that(".nc_t_quantile() is qt() inside the range where qt() is accurate", {
  q <- DMAR:::.nc_t_quantile
  for (df in c(1, 10, 100, 1000)) for (ncp in c(-20, -1, 2, 30, 37.62))
    for (p in c(.001, .025, .5, .975, .999)) {
      expect_equal(q(p, df, ncp), suppressWarnings(stats::qt(p, df, ncp)),
                   tolerance = 1e-11, info = paste(p, df, ncp))
    }
  # qt() warns for some of these (its search meets a probability near 1);
  # .nc_t_quantile() does not.
  expect_match(tryCatch(stats::qt(.975, 2000, 37), warning = conditionMessage),
               "full precision")
  expect_no_warning(q(.975, 2000, 37))
  expect_identical(q(0, 5, 2), -Inf)
  expect_identical(q(1, 5, 2), Inf)
  expect_identical(q(.3, 7, 0), stats::qt(.3, 7))
  expect_identical(q(.3, 2e6, 0), stats::qt(.3, 2e6))
  expect_equal(q(.3, Inf, 2), stats::qnorm(.3, 2))
})

test_that(".nc_t_power() is the power past 37.62, by the oracle and by simulation", {
  pw <- DMAR:::.nc_t_power
  crit <- stats::qt(1 - .001 / 2, 2)
  # Two groups of 2 at d = 40 (ncp 40), alpha .001: pt() gave .7823613; 4e6
  # simulated pooled-variance t tests rejected .798202 of the time (MCSE
  # .000201), and .959419 (MCSE .000099) one-sided, where pt() gave .9660579
  # (evidence 08).
  two <- pw(crit, 2, 40)
  expect_equal(two, 0.7981439583, tolerance = 1e-9)
  expect_equal(two, .oracle_nct(crit, 2, 40, lower_tail = FALSE) +
                 .oracle_nct(-crit, 2, 40), tolerance = 1e-12)
  one <- pw(stats::qt(1 - .001, 2), 2, 40, two_sided = FALSE)
  expect_equal(one, 0.9591889337, tolerance = 1e-9)
  # The sign of the effect does not matter.
  expect_identical(pw(crit, 2, -40), two)
  # One degree of freedom, ncp 41.9 (the composite 2-group planner's
  # marginal power below): pt() gave .1666, the exact power is .0525.
  expect_equal(pw(stats::qt(1 - .001 / 2, 1), 1, 41.93139347),
               .oracle_nct(stats::qt(1 - .001 / 2, 1), 1, 41.93139347, lower_tail = FALSE) +
                 .oracle_nct(-stats::qt(1 - .001 / 2, 1), 1, 41.93139347),
               tolerance = 1e-12)
})

test_that(".power_noncentral_t() is the planners' earlier power inside the range of pt()", {
  # The formula ss_power_smd() and its siblings used through 3338679.
  old <- function(ncp, df, a, dir) {
    if (dir) {
      crit <- stats::qt(1 - a, df)
      if (ncp >= 0) 1 - stats::pt(crit, df, ncp) else stats::pt(-crit, df, ncp)
    } else {
      (1 - stats::pt(stats::qt(1 - a / 2, df), df, ncp)) +
        stats::pt(stats::qt(a / 2, df), df, ncp)
    }
  }
  for (df in c(2, 10, 100, 2000, 1e5)) for (ncp in c(-30, -2, 0, 0.5, 3, 33))
    for (a in c(.05, .001)) for (dir in c(FALSE, TRUE)) {
      expect_lt(abs(DMAR:::.power_noncentral_t(ncp, df, a, dir) -
                      suppressWarnings(old(ncp, df, a, dir))), 1e-13,
                label = paste(df, ncp, a, dir))
    }
})

test_that(".nc_t_density() is the noncentral t density, by the oracle", {
  dens <- DMAR:::.nc_t_density
  # Past 37.62, with many df near x = 0, and with more than 1,000 df, where
  # dt() was off by 5.5 percent, 2 percent, .26 percent, and 2.3e-10.
  for (x in list(c(50, 3, 50), c(45, 1000, 41), c(1e-3, 1e5, 2), c(5, 5000, 3),
                 c(0.05, 30, 1), c(-50, 3, -50), c(-0.5, 2e4, 1))) {
    expect_equal(dens(x[1], x[2], x[3]), .oracle_nct_density(x[1], x[2], x[3]),
                 tolerance = 1e-10, info = paste(x, collapse = " "))
  }
  expect_equal(dens(0, 10, 3), .oracle_nct_density(0, 10, 3), tolerance = 1e-14)
  expect_equal(dens(2, Inf, 1), stats::dnorm(2, 1))
  expect_identical(dens(2, 7, 0), stats::dt(2, 7))
  # Vectorized over x and over ncp.
  expect_equal(dens(c(1, 40, 60), 3, 45),
               vapply(c(1, 40, 60), dens, 0, df = 3, ncp = 45))
  expect_equal(dens(40, 49, c(0, 30, 50)),
               vapply(c(0, 30, 50), function(n) dens(40, 49, n), 0))
  expect_identical(dens(c(-Inf, Inf), 3, 45), c(0, 0))
})

test_that(".nc_t_density() is dt() where dt() is accurate, and does not warn", {
  dens <- DMAR:::.nc_t_density
  accurate <- DMAR:::.nc_t_dt_accurate
  expect_true(all(accurate(c(-5, 0.1, 30), 1000, 37.62)))
  expect_false(accurate(0.09, 10, 1))
  expect_false(accurate(1, 1001, 1))
  expect_false(accurate(1, 10, 37.63))
  expect_true(accurate(1e-5, 1e6, 0))
  # dt() warns here (its lower-tail pnt() sum is within 1e-10 of 1);
  # .nc_t_dt() asks for the upper tails and returns the same number.
  expect_match(tryCatch(stats::dt(30, 10, 0.5), warning = conditionMessage),
               "full precision")
  expect_identical(expect_no_warning(dens(30, 10, 0.5)),
                   suppressWarnings(stats::dt(30, 10, 0.5)))
  x <- c(-40, -3, -0.1, 0.1, 0.5, 2, 8, 30, 200, 1e4)
  for (df in c(1, 3, 30, 1000)) for (ncp in c(-37.62, -5, 0.5, 2, 20, 37.62)) {
    expect_identical(expect_no_warning(DMAR:::.nc_t_dt(x, df, ncp)),
                     suppressWarnings(stats::dt(x, df, ncp)),
                     info = paste(df, ncp))
    expect_identical(dens(x, df, ncp), suppressWarnings(stats::dt(x, df, ncp)))
  }
})

test_that("fewer than one degree of freedom is integrated, and agrees with two oracles", {
  # With 0 < df < 1 the density of S = sqrt(V / df) is infinite at 0, and
  # integrate() over s stopped with an unclassed "the integral is probably
  # divergent" (cv_t(.001, 0.5, ncp = 38) and cv_t(.05, 0.2, ncp = 45));
  # pt() with ncp was off inside 37.62 by as much as 3.8e-5 at df = 0.5 and
  # .14 at df = 0.05 (see .nc_t_pt_limit()). The integral is now over
  # u = s^df, and pt() is not used for a nonzero noncentrality there.
  expect_identical(DMAR:::.nc_t_pt_limit(0.5), -Inf)
  expect_identical(DMAR:::.nc_t_pt_limit(1), 37.62)
  expect_false(DMAR:::.nc_t_pt_accurate(0.99, 0.3))
  expect_true(DMAR:::.nc_t_pt_accurate(0.5, 0))
  expect_false(DMAR:::.nc_t_dt_accurate(1, 0.5, 3))
  prob <- DMAR:::.nc_t_prob
  for (x in list(c(-9.1565074e7, .5, -37.62), c(1e6, .5, 38), c(10, .5, 38),
                 c(-1e3, .2, 45), c(3, .8, -1), c(60, .1, 40), c(.5, .9, 2),
                 c(-1e5, .99, 5), c(1e9, .3, -.3))) {
    for (lower in c(TRUE, FALSE)) {
      value <- prob(x[1], x[2], x[3], lower_tail = lower)
      expect_lt(abs(value - .oracle_nct(x[1], x[2], x[3], lower)), 1e-12,
                label = paste(c(x, lower), collapse = " "))
      series <- .oracle_nct_series(x[1], x[2], x[3], lower)
      if (!is.na(series) && series > 1e-4)
        expect_lt(abs(value / series - 1), 1e-8, label = paste(c(x, lower), collapse = " "))
    }
  }
  # pt() inside 37.62 at the .0005 quantile of df = 0.5, ncp = -37.62:
  # 4.66e-4 where the probability is 5e-4.
  q <- -9.1565073644e7
  expect_equal(prob(q, .5, -37.62), 5e-4, tolerance = 1e-8)
  expect_equal(.oracle_nct_series(q, .5, -37.62), 5e-4, tolerance = 1e-8)
  expect_gt(abs(stats::pt(q, .5, -37.62) - 5e-4), 3e-5)
  # And far worse with fewer degrees of freedom, away from the far tails:
  # at the .2 quantile of df = 0.1, ncp = -20, pt() gives .130.
  q <- DMAR:::.nc_t_quantile(.2, .1, -20)
  expect_equal(.oracle_nct(q, .1, -20), .2, tolerance = 1e-9)
  expect_equal(stats::pt(q, .1, -20), .130, tolerance = 1e-2)
  # The density, against the conditioning oracle.
  for (x in list(c(-5, .5, -3), c(2, .2, 1), c(35, .9, 30), c(1e4, .5, 40), c(.05, .3, 1)))
    expect_lt(abs(DMAR:::.nc_t_density(x[1], x[2], x[3]) -
                    .oracle_nct_density(x[1], x[2], x[3])), 1e-12,
              label = paste(x, collapse = " "))
})

test_that("more than 1e7 degrees of freedom are integrated over the standardized chi square, by two oracles", {
  # Up to 1e12 degrees of freedom the oracle that conditions on Z; beyond,
  # the expansion of E{Phi(q S - ncp)} about S = 1 to second order, whose
  # error is of order q^4 / df^2 (below 1e-17 here). The integral over s
  # stopped with "roundoff error was detected" at 6.1e15 degrees of freedom.
  expansion <- function(q, df, ncp, lower_tail) {
    x <- q - ncp
    corr <- stats::dnorm(x) * (-q / (4 * df) - x * q^2 / (4 * df))
    if (lower_tail) stats::pnorm(x) + corr else stats::pnorm(-x) - corr
  }
  for (pt in list(c(2e7, 2, 0.5), c(1e9, 10, 12.5), c(1e11, -3, -1.5),
                  c(1e12, 37, 39))) {
    for (lt in c(TRUE, FALSE))
      expect_equal(DMAR:::.nc_t_cdf(pt[3], pt[1], pt[2], lower_tail = lt),
                   .oracle_nct(pt[3], pt[1], pt[2], lower_tail = lt),
                   tolerance = 1e-10, label = paste(pt, lt, collapse = " "))
  }
  for (pt in list(c(6.146334e15, 1.959964, 0), c(1.8e16, 0, 4.75),
                  c(1e14, 100, 95.25), c(1e13, -2, 2.9))) {
    for (lt in c(TRUE, FALSE)) {
      v <- DMAR:::.nc_t_cdf(pt[3], pt[1], pt[2], lower_tail = lt)
      o <- expansion(pt[3], pt[1], pt[2], lt)
      expect_lt(abs(v - o), 1e-13 + 1e-10 * min(o, 1 - o))
    }
  }
})

test_that("with more than 1e7 degrees of freedom the integral is accurate where q is large beside the spread of S", {
  # The integrand was given a rounded s = sqrt(1 + h), so q s - ncp moved
  # in steps of q times 1.1e-16 or 2.2e-16 (the spacing of doubles just
  # below and just above 1): where q sd(S) is .3 or more, integrate()
  # stopped with "roundoff error was detected" (dmar_computation_failed),
  # and where it returned the probability was off by as much as 2.1e-10
  # (1.2e-10 at the second point). With 4.1e31 degrees of freedom and
  # q = 2.7e15 those steps are .30 and .60 wide, and P(T <= 2.7e15) came
  # out .131 where it is .169. The range of the integral from qchisq() was
  # .31 of a standard deviation of the chi square wider than the 1e-14
  # quantiles at each end there, which does no harm: with s not rounded
  # and that range the probability is within 8e-16 of the pin
  # (work/round6/evidence/m4_nct, fix3/range_vs_s.log). Each
  # row is df, q (also k and x), ncp, and the 256-bit values of P(T <= q),
  # P(T > q), P(|T| <= q), and the density at q (work/round6/evidence/
  # m4_nct, test_pins.R, an integral over the standardized chi square with
  # its density from the log gamma function and q s - ncp formed in 256
  # bits). q sd(S) runs from .3 to 10. The integral leaves out the chi
  # square beyond its 1e-14 and 1 - 1e-14 quantiles and has an absolute
  # tolerance of 1e-13 a piece, so the values are compared to 1e-13.
  pins <- rbind(
    c(1e+14, 14142135.62, 14142138.62, 0.016947426312576442, 0.98305257368742338,
      0.016947426312576442, 0.029732574263189933),
    c(1e+14, 1.4e+07, 14000001, 0.23864466179665, 0.76135533820334977,
      0.23864466179665, 0.2202455546253026),
    c(1.8e+16, 56920997.88, 56920997.88, 0.49999999968959513, 0.5000000003104047,
      0.49999999968959513, 0.3821174024562391),
    c(1.8e+16, 569209978.8, 569209982, 0.15578603244706893, 0.84421396755293088,
      0.15578603244706893, 0.075605133555292983),
    c(1e+10, 1414213.562, 1414183, 0.99882118244980111, 0.0011788175501986583,
      0.99882118244980111, 0.00038956829245084563),
    c(1e+18, 1.4e+10, 13999999998, 0.57965366966618459, 0.42034633033381524,
      0.57965366966618459, 0.039293330634082128),
    c(4.1e+31, 2.7e+15, 2700000000000001, 0.16895351243256118, 0.8310464875674386,
      0.16895351243256118, 0.24154448527999231))
  for (i in seq_len(nrow(pins))) {
    p <- pins[i, ]
    label <- paste(format(p[1:3], digits = 12), collapse = " ")
    expect_lt(abs(DMAR:::.nc_t_cdf(p[2], p[1], p[3]) - p[4]), 1e-13, label = label)
    expect_lt(abs(DMAR:::.nc_t_cdf(p[2], p[1], p[3], lower_tail = FALSE) - p[5]),
              1e-13, label = label)
    expect_lt(abs(DMAR:::.nc_t_abs_probability(p[2], p[1], p[3]) - p[6]), 1e-13,
              label = label)
    expect_lt(abs(DMAR:::.nc_t_density_integral(p[2], p[1], p[3]) - p[7]), 1e-13,
              label = label)
  }
  # The quantile, which stopped the same way at 1e15 degrees of freedom;
  # its search stops within 1e-12 |ncp| (2.2e-5 here) of the root, where
  # the density of T is below .4.
  ncp <- 0.5 * sqrt(2e15)
  q <- DMAR:::.nc_t_quantile(.975, 1e15, ncp)
  expect_lt(abs(DMAR:::.nc_t_cdf(q, 1e15, ncp) - .975), .4e-12 * ncp)
})

test_that("the density at 0 is the ratio of gamma functions to 1e-14 at every df", {
  # exp(-ncp^2 / 2) / {B(df / 2, 1/2) sqrt(df)}, against 256-bit values
  # of exp(lgamma((df + 1) / 2) - lgamma(df / 2) - log(pi df) / 2 -
  # ncp^2 / 2) (work/round6/evidence/m4_nct, dens0_pins.R). Differencing
  # the two log gammas, the form used up to 1e7 degrees of freedom, was off
  # by 3.8e-11 at 1e5, 4.3e-9 at 5.6e6, and 9.7e-9 at 9.99e6.
  pins <- rbind(c(0.05, 0.035095627272231172), c(1, 0.10334008993412312),
                c(1e5, 0.1295172718723073), c(5.6e6, 0.12951758988385634),
                c(9.99e6, 0.1295175924247107), c(1e16, 0.12951759566589172))
  for (i in seq_len(nrow(pins)))
    expect_equal(DMAR:::.nc_t_density_integral(0, pins[i, 1], 1.5), pins[i, 2],
                 tolerance = 1e-14, label = format(pins[i, 1]))
})

test_that("a failure of the numerical integral is refused with a classed condition naming df", {
  e <- tryCatch(DMAR:::.nc_t_integral(function(s, d) rep(NaN, length(s)), 0.5),
                error = function(e) e)
  expect_s3_class(e, "dmar_computation_failed")
  expect_identical(e$df, 0.5)
  expect_match(conditionMessage(e), "df = 0.5", fixed = TRUE)
  e <- tryCatch(DMAR:::.nc_t_integral(function(s, d) rep(NaN, length(s)), 7),
                error = function(e) e)
  expect_s3_class(e, "dmar_computation_failed")
  expect_identical(e$df, 7)
  e <- tryCatch(DMAR:::.nc_t_integral(function(s, d) rep(NaN, length(s)), 1e12),
                error = function(e) e)
  expect_s3_class(e, "dmar_computation_failed")
  expect_identical(e$df, 1e12)
})

test_that("a quantile search that fails past the largest double is refused with a classed condition naming df", {
  far <- DMAR:::.nc_t_quantile_far
  stop_e <- simpleError("no sign change found in 1000 iterations")
  # A tail probability that cannot be computed.
  e <- tryCatch(far(function(q) NA_real_, .5, .001, 5, c(-1, 1), 1e-12, stop_e),
                error = function(e) e)
  expect_s3_class(e, "dmar_computation_failed")
  expect_identical(e$df, .001)
  expect_match(conditionMessage(e), "df = 0.001", fixed = TRUE)
  # A function that is not increasing as a distribution function is.
  e <- tryCatch(far(function(q) if (abs(q) >= 1e308) sign(q) else 0, .5, .01, 5,
                    c(-1, 1), 1e-12, stop_e), error = function(e) e)
  expect_s3_class(e, "dmar_computation_failed")
  expect_match(conditionMessage(e), "no sign change", fixed = TRUE)
  # -Inf and Inf when the tail at the largest double reaches the target.
  expect_identical(far(function(q) if (q < 0) .1 else 1, .2, .01, 5, c(-1, 1), 1e-12, stop_e), -Inf)
  expect_identical(far(function(q) if (q > 0) -.1 else -1, .8, .01, 5, c(-1, 1), 1e-12, stop_e), Inf)
  # A refusal of the integral inside the search is passed on as it is, not
  # taken for a failure of the bracket.
  local_mocked_bindings(
    .nc_t_prob = function(...) DMAR:::.dmar_refuse("the integral failed", "computation_failed", df = .5),
    .package = "DMAR")
  e <- tryCatch(DMAR:::.nc_t_quantile(.3, .5, 2), error = function(e) e)
  expect_s3_class(e, "dmar_computation_failed")
  expect_identical(conditionMessage(e), "the integral failed")
})
