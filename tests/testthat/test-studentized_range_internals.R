# Tests for the studentized range engine shared by cv_tukey_hsd(),
# ci_tukey_kramer(), ci_games_howell(), contrast_test(adjust = "tukey"), and
# the Bryant-Paulson functions (R/studentized_range_internals.R).
#
# Oracles, none of them package code:
#   - the closed form for two means, Q = sqrt(2) |T| on nu df;
#   - an independent quadrature of the studentized range (the 1.1.0 review's
#     F3 evidence, oracle_studentized_range.R: the upper tail integrated
#     directly, nested
#     adaptive quadrature in the chi variable S with the range integral done
#     afresh at every node, no spline), whose values are pinned below to the
#     digits printed there; two of them were reproduced by a second engine,
#     mpmath's arbitrary precision quadrature;
#   - a second independent quadrature (the 1.1.0 review's G2 evidence,
#     oracle_ts.R: the minimum of the k normals and the chi square variate
#     replaced by their probability integral transforms, both integrals by a
#     tanh-sinh rule, the tail and the distribution function each integrated
#     from its own integrand), which agrees with the first to 1e-11 where both
#     were run and supplies the far-tail and many-means values below;
#   - the studentized range table abridged from Table 29 of Pearson and
#     Hartley's Biometrika Tables for Statisticians, as printed in the earlier
#     appendix of Maxwell, Delaney, and Kelley (Table A.4);
#   - a raw-data simulation of the studentized range.

test_that("two means: the quantile and the distribution function are the closed forms", {
  for (nu in c(1, 1.5, 2, 3, 5, 8.5)) for (p in c(0.90, 0.95, 0.99, 0.999)) {
    expect_equal(.sr_quantile(p, 2, nu), sqrt(2) * qt(1 - (1 - p) / 2, nu),
                 tolerance = 1e-12, info = paste(nu, p))
  }
  for (nu in c(1, 2, 4)) for (q in c(0.5, 3, 20)) {
    expect_equal(.sr_cdf(q, 2, nu), 2 * pt(q / sqrt(2), nu) - 1,
                 tolerance = 1e-12, info = paste(nu, q))
  }
})

test_that("the exact quantiles at 2 and 3 error df match the independent quadrature", {
  # stats::qtukey() gives 44.526266 and 19.947761 here.
  expect_equal(.sr_quantile(0.99, 20, 2), 37.943462287637, tolerance = 1e-8)
  expect_equal(.sr_quantile(0.99, 20, 3), 19.764809074501, tolerance = 1e-8)
  expect_gt(qtukey(0.99, 20, 2) / 37.943462287637 - 1, 0.17)
  expect_gt(qtukey(0.99, 20, 3) / 19.764809074501 - 1, 0.009)
  # The upper tail at qtukey()'s value, by the same quadrature: the level
  # qtukey() delivers is .00727, not .01.
  expect_equal(.sr_cdf(qtukey(0.99, 20, 2), 20, 2, lower_tail = FALSE),
               0.0072732841435454, tolerance = 1e-7)
})

test_that("further small-df quantiles match the independent quadrature", {
  skip_on_cran()  # a spline of the range distribution per number of means
  # Columns: k, df, p, oracle quantile. qtukey() returns NaN at df = 1 and
  # 1.5, is 19.914 at (50, 2, .95), 22.150 at (10, 4, .999), 27.757 at
  # (6, 2, .99), and 11.440122 at (50, 12, .999).
  ref <- rbind(c( 5,   1, 0.95,  37.081501903677),
               c(10,   1, 0.99, 245.541642670005),
               c( 4, 1.5, 0.95,  14.40883847907),
               c(50,   2, 0.95,  20.046456412062),
               c( 6,   2, 0.99,  26.629041325501),
               c(10,   4, 0.999, 22.328431008143),
               c( 4,   6, 0.999, 10.965012923597),
               c(30,  12, 0.999, 10.675949325378),
               c(50,  12, 0.999, 11.440605620077))
  for (i in seq_len(nrow(ref))) {
    expect_equal(.sr_quantile(ref[i, 3], ref[i, 1], ref[i, 2]), ref[i, 4],
                 tolerance = 1e-8, info = paste(ref[i, 1:3], collapse = ", "))
  }
})

test_that("above 10000 df the exact quantile replaces qtukey()'s value", {
  skip_on_cran()  # a spline of the range distribution for ten means
  # From df = 25001 on, stats::qtukey(.95, 10, df) returns the infinite-df
  # value 4.474124; the quadrature gives the values below.
  expect_equal(.sr_quantile(0.95, 10, 30000), 4.474463176805, tolerance = 1e-8)
  expect_equal(.sr_quantile(0.99, 10, 1e5), 5.156803194664, tolerance = 1e-8)
  expect_equal(qtukey(0.95, 10, 30000), qtukey(0.95, 10, Inf))
  # Below 25001, ptukey()'s error grows from 15000 on: at 18000 df
  # qtukey(.999, 3, 18000) is 5.064422809 where the tanh-sinh quadrature gives
  # 5.0644103381282 (1.2e-5 off), and qtukey(.975, 100, 15000) is 6.334793825
  # where it gives 6.3346825564786 (1.1e-4 off).
  expect_true(all(.sr_needs_exact(c(3, 100), c(18000, 15000))))
  expect_equal(.sr_quantile(0.999, 3, 18000), 5.0644103381282, tolerance = 1e-9)
  expect_equal(.sr_quantile(0.975, 100, 15000), 6.3346825564786,
               tolerance = 1e-9)
  # At 10000 df and below, R's value is kept.
  expect_false(.sr_needs_exact(10, 10000))
  expect_identical(.sr_quantile(0.95, 10, 10000), qtukey(0.95, 10, 10000))
})

test_that("the large-df region ends at 1e8 df, where the exact path still holds", {
  skip_on_cran()  # a spline of the range distribution for 50 means
  expect_true(all(.sr_needs_exact(c(4, 50, 500), 1e8)))
  expect_false(any(.sr_needs_exact(c(4, 50, 500), 1e8 + 1)))
  # Two independent quadratures (the lane's and the reviewer's) agree to
  # 5e-13 at these cells; qtukey() returns the infinite-df values 3.633160 and
  # 6.940585 here.
  expect_equal(.sr_quantile(0.95, 4, 1e8), 3.63315963598318, tolerance = 1e-9)
  expect_equal(.sr_quantile(0.999, 50, 1e8), 6.94058540522843, tolerance = 1e-9)
})

test_that("above 1e8 df the values are qtukey()'s and ptukey()'s infinite-df values", {
  # From 25001 df on, R returns the infinite-df distribution, which at 1e8 df
  # is within 1e-6 of the exact one. The exact path is not used beyond 1e8:
  # uncapped, it broke down at huge df, where qchisq(1e-13, df) / df rounds
  # to 1 (at 1e40 df its .95 quantile for 4 means came out 11 rather than
  # 3.633, and cv_tukey_hsd(.05, df = 1e40, n_groups = 4) returned 7.778).
  p <- c(0.9, 0.95, 0.99, 0.999); q <- c(2, 3.633, 5)
  for (nu in c(1e8 + 1, 1e9, 1e40, 1e300)) {
    expect_false(.sr_needs_exact(4, nu))
    expect_identical(.sr_quantile(p, 4, nu), qtukey(p, 4, Inf))
    expect_identical(.sr_cdf(q, 4, nu), ptukey(q, 4, Inf))
    expect_identical(.sr_cdf(q, 4, nu, lower_tail = FALSE),
                     ptukey(q, 4, Inf, lower.tail = FALSE))
  }
  # The infinite-df .95 quantile for 4 means by both independent quadratures
  # is 3.63315957490263.
  expect_equal(.sr_quantile(0.95, 4, 1e40), 3.63315957490263, tolerance = 1e-8)
  for (nu in c(1e40, 1e300)) {
    out <- cv_tukey_hsd(.05, df = nu, n_groups = 4)
    expect_equal(out$value, 3.63315957490263 / sqrt(2), tolerance = 1e-8)
    expect_equal(out$area_less, 0.95, tolerance = 1e-8)
    expect_equal(out$area_greater, 0.05, tolerance = 1e-7)
    expect_equal(qbryant_paulson(0.95, n_covariates = 0, n_groups = 4, df = nu),
                 3.63315957490263, tolerance = 1e-8)
  }
  # The forced exact engine treats such df as infinite rather than failing.
  expect_identical(.sr_quantile_exact_one(0.95, 4, 1e300),
                   .sr_quantile_exact_one(0.95, 4, Inf))
  expect_identical(.sr_cdf_exact_one(3.6, 4, 1e40), .sr_cdf_exact_one(3.6, 4, Inf))
})

test_that("where qtukey() fails above 1e8 df, the infinite-df exact root is used", {
  skip_on_cran()  # a spline of the range distribution for 50 means
  # qtukey(.5, 50, df) stops with "convergence failed" at any df; the
  # infinite-df median for 50 means by both independent quadratures is
  # 4.45048133802192.
  expect_true(is.nan(suppressWarnings(qtukey(0.5, 50, 1e40))))
  for (nu in c(1e9, 1e40, Inf)) {
    expect_no_warning(v <- .sr_quantile(0.5, 50, nu))
    expect_equal(v, 4.45048133802192, tolerance = 1e-9)
  }
})

test_that("where qtukey() and ptukey() are accurate, DMAR returns them unchanged", {
  # Outside the regions .sr_needs_exact() and .sr_needs_exact_p() mark, the
  # values are R's own (the check on qtukey()'s value passes at these cells).
  # With more than 50 means the distribution function is R's only from its
  # .75 point up (.80 with more than 1000 means), so the probabilities are
  # taken above it.
  cells <- rbind(c(3, 27), c(3, 9), c(6, 9), c(7, 11), c(10, 11), c(11, 15),
                 c(20, 15), c(21, 21), c(50, 21), c(51, 31), c(101, 41),
                 c(251, 47), c(1001, 81), c(5001, 111),
                 c(5, 1000), c(5, 10000), c(5, Inf))
  for (i in seq_len(nrow(cells))) {
    k <- cells[i, 1]; nu <- cells[i, 2]
    expect_false(.sr_needs_exact(k, nu))
    p <- c(0.9, 0.95, 0.99)
    expect_identical(.sr_quantile(p, k, nu), qtukey(p, k, nu))
    q <- if (k > 50) qtukey(c(0.85, 0.9, 0.99), k, nu) else c(2, 4, 6)
    expect_identical(.sr_cdf(q, k, nu), ptukey(q, k, nu))
  }
  # One step inside the boundary of each row the exact computation is used.
  expect_true(all(.sr_needs_exact(
    c(2, 6, 7, 10, 11, 20, 21, 50, 51, 100, 101, 250, 251, 400, 401, 700, 701,
      1000, 1001, 1500, 1501, 2000, 2001, 3000, 3001, 5000, 5001),
    c(8, 8, 10, 10, 14, 14, 20, 20, 30, 30, 40, 40, 46, 46, 56, 56, 70, 70,
      80, 80, 90, 90, 100, 100, 110, 110, 110))))
  expect_true(.sr_needs_exact(5, 10001))
})

test_that("with 250 to 5000 means the region follows the measurements", {
  skip_on_cran()  # splines of the range distribution for 1000 and 200 means
  # 1000 means at 41 df (inside the region, which now reaches 70 df there):
  # qtukey(.999, 1000, 41) is 10.299349568 and qtukey(.995, 1000, 41)
  # 9.486225211, where both independent quadratures give 10.2994946166 and
  # 9.4861526011 (1.45e-4 and 7.3e-5 off).
  expect_true(.sr_needs_exact(1000, 41))
  expect_equal(.sr_quantile(c(0.999, 0.995), 1000, 41),
               c(10.2994946166086, 9.4861526010680), tolerance = 1e-9)
  # At 70 df, the last measured df of the 1000-means row, qtukey(.999) is
  # 9.509346417 and the quadrature 9.5093521652895 (5.7e-6 off); at 75 df,
  # outside the region, qtukey() is 1.8e-6 off and kept.
  expect_true(.sr_needs_exact(1000, 70))
  expect_equal(.sr_quantile(0.999, 1000, 70), 9.5093521652895, tolerance = 1e-9)
  expect_false(.sr_needs_exact(1000, 75))
  expect_identical(.sr_quantile(0.999, 1000, 75), qtukey(0.999, 1000, 75))
  # With more than 50 means ptukey() is off below its .75 point: at 200
  # means and 1000 df its value at the quadrature's .05 quantile,
  # 4.6306635297958, is .0499975665 (2.4e-6 off). The exact distribution
  # function is used there, for the probability and for the quantile.
  expect_lt(ptukey(4.6306635297958, 200, 1000), 0.05 - 2e-6)
  expect_equal(.sr_cdf(4.6306635297958, 200, 1000), 0.05, tolerance = 1e-8)
  expect_equal(.sr_cdf(4.6306635297958, 200, 1000, lower_tail = FALSE), 0.95,
               tolerance = 1e-8)
  expect_equal(.sr_quantile(0.05, 200, 1000), 4.6306635297958, tolerance = 1e-9)
})

test_that("with more than 1000 means the probability rule reaches the .80 point", {
  # Measured (G2 round 2): at infinite df, which R also uses for every df
  # above 25000, ptukey() is more than 1e-6 off just above its .75 point from
  # about 4500 means on, and just above 1e8 df, where DMAR returns R's
  # infinite-df value, from 2500 means on; at .80 and above it was within
  # 4.9e-7 with up to 10000 means.
  expect_true(.sr_needs_exact_p(0.7499, 1000))
  expect_false(.sr_needs_exact_p(0.79, 1000))
  expect_true(.sr_needs_exact_p(0.7999, 1001))
  expect_false(.sr_needs_exact_p(0.80, 1001))
  expect_equal(.sr_needs_exact_p(c(0.70, 0.78, 0.85, 0.78), c(50, 2000, 2000, 1000)),
               c(FALSE, TRUE, FALSE, FALSE))
})

test_that("the range integral finds the narrow peak of its integrand with many means", {
  # With 7000 means the integrand of the range's distribution function at
  # w = 6.218556 is a peak with a standard deviation of 0.08 just below
  # w / 2. Integrated over [-8.5, 8.5] in one piece it was missed (2.3e-17),
  # which put a wrong knot in the spline and moved the .00001 quantile at
  # infinite df by 6.7e-4. The tanh-sinh quadrature gives 9.546014059722e-6.
  expect_lt(abs(.sr_range_cdf(6.218556, 7000) / 9.546014059722e-06 - 1), 1e-9)
})

test_that("with thousands of means the distribution function is exact below .80", {
  skip_on_cran()  # splines of the range distribution for 300 to 7000 means
  # At infinite df ptukey() puts 7.6250574053 at .750001043 with 5000 means,
  # where the tanh-sinh quadrature gives .749999999975 (to 1e-15 at two
  # step sizes); R uses the infinite-df value at every df above 25000.
  expect_gt(ptukey(7.6250574053, 5000, Inf) - 0.75, 1e-6)
  expect_identical(ptukey(7.6250574053, 5000, 1e9), ptukey(7.6250574053, 5000, Inf))
  for (nu in c(Inf, 1e9)) {
    expect_lt(abs(.sr_cdf(7.6250574053, 5000, nu) - 0.75), 1e-8)
    expect_lt(abs(.sr_cdf(7.6250574053, 5000, nu) - 0.749999999975), 1e-9)
    expect_lt(abs(.sr_cdf(7.6250574053, 5000, nu, lower_tail = FALSE) - 0.25), 1e-8)
  }
  # Just above the 1e8 cap, 3000 means: R's infinite-df value is 1.07e-6 off
  # at the .75 point of the distribution at 1e8 + 1 df (the quadrature's
  # quantile there is 7.3678768376933); the exact infinite-df distribution,
  # which DMAR uses there, is 1.6e-7 off.
  q <- 7.3678768376933
  expect_gt(abs(ptukey(q, 3000, 1e8 + 1) - 0.75), 1e-6)
  expect_lt(abs(.sr_cdf(q, 3000, 1e8 + 1) - 0.75), 3e-7)
  # The far lower tail: before the range integral was split, the spline
  # knots missed by the one-piece integral put the .00001 quantile 6.7e-4
  # high with 7000 means at infinite df and 1.7e-6 high with 300 means at
  # 1000 df. The quadrature gives 6.2207266518362 and 3.9950894055082.
  expect_lt(abs(.sr_quantile(1e-5, 7000, Inf) - 6.2207266518362), 1e-6)
  expect_lt(abs(.sr_quantile(1e-5, 300, 1000) - 3.9950894055082), 1e-7)
})

test_that("beyond .999 and below .001 the quantile is the exact root at any df", {
  skip_on_cran()  # splines of the range distribution for 3, 10, and 2000 means
  expect_true(all(.sr_needs_exact_p(c(0.0009, 0.9991, 0.9999), 3)))
  expect_false(any(.sr_needs_exact_p(c(0.001, 0.5, 0.999), 50)))
  expect_true(all(.sr_needs_exact_p(c(0.001, 0.5, 0.7499), 51)))
  expect_false(.sr_needs_exact_p(0.75, 51))
  # Upper tail: qtukey(.9999, 10, 11) is 12.008462418 and
  # qtukey(.99999, 3, 95) 6.984008738 where the tanh-sinh quadrature gives
  # 12.008579776018 (review 2's quadrature: 12.008579776) and 6.9839926286659;
  # qtukey(.99999, 2000, 1000) is 12.929769658 where it gives 9.8717782115797,
  # and at infinite df qtukey(.99999, 2000, Inf) is about 13.75 against
  # 9.754723051.
  expect_equal(.sr_quantile(0.9999, 10, 11), 12.0085797760180, tolerance = 1e-8)
  expect_gt(abs(qtukey(0.9999, 10, 11) - 12.0085797760180), 1e-4)
  expect_equal(.sr_quantile(0.99999, 3, 95), 6.9839926286659, tolerance = 1e-8)
  expect_equal(.sr_quantile(0.99999, 2000, 1000), 9.8717782115797,
               tolerance = 1e-8)
  expect_gt(qtukey(0.99999, 2000, 1000) - 9.8717782115797, 3)
  expect_equal(.sr_quantile(0.99999, 2000, Inf), 9.754723051, tolerance = 1e-8)
  expect_gt(qtukey(0.99999, 2000, Inf) - 9.754723051, 3)
  expect_equal(.sr_quantile(0.99999, 2000, 1e9), 9.754723051, tolerance = 1e-8)
  # Lower tail: qtukey(.00001, 3, 150) is 9.2e-5 where the quadrature gives
  # 0.0060229704045, and at infinite df 0.0060229704 as well (to 1e-9).
  expect_lt(qtukey(1e-5, 3, 150), 1e-4)
  expect_equal(.sr_quantile(1e-5, 3, 150), 0.0060229704045, tolerance = 1e-6)
  expect_equal(.sr_quantile(1e-5, 3, Inf), 0.00602297, tolerance = 1e-6)
  # Two means: the closed form, where qtukey(.99999, 2, 5) is 0.48 too large.
  expect_equal(.sr_quantile(0.99999, 2, 5), sqrt(2) * qt(1 - 0.00001 / 2, 5),
               tolerance = 1e-12)
  expect_gt(qtukey(0.99999, 2, 5) - sqrt(2) * qt(1 - 0.00001 / 2, 5), 0.4)
})

test_that("a qtukey() value that stopped short of ptukey()'s root is replaced", {
  skip_on_cran()  # splines of the range distribution for 5, 10, and 150 means
  # qtukey()'s secant search can stop short of ptukey()'s root. Outside the
  # df region, at these cells, it did: qtukey(.975, 10, 12) is
  # 6.003766474, qtukey(.9975, 5, 23) 6.018665414, qtukey(.995, 5, 17)
  # 5.903134911, and qtukey(.975, 150, 300) 6.671163653, where both
  # independent quadratures give (the first) 6.003807292 and the tanh-sinh
  # quadrature 6.0185456918079, 5.9030792766994, and 6.6710808618315.
  ref <- rbind(c(10, 12, 0.975, 6.0038072922447),
               c( 5, 23, 0.9975, 6.0185456918079),
               c( 5, 17, 0.995, 5.9030792766994),
               c(150, 300, 0.975, 6.6710808618315))
  for (i in seq_len(nrow(ref))) {
    k <- ref[i, 1]; nu <- ref[i, 2]; p <- ref[i, 3]
    expect_false(.sr_needs_exact(k, nu))
    expect_false(.sr_needs_exact_p(p, k))
    expect_gt(abs(qtukey(p, k, nu) - ref[i, 4]), 4e-5)
    expect_gt(abs(.sr_qtukey_offset(qtukey(p, k, nu), p, k, nu)), 2e-6)
    expect_equal(.sr_quantile(p, k, nu), ref[i, 4], tolerance = 1e-9)
  }
  # Where qtukey() reached ptukey()'s root, the check passes and R's value is
  # returned unchanged.
  expect_lt(abs(.sr_qtukey_offset(qtukey(0.95, 10, 12), 0.95, 10, 12)), 2e-6)
  expect_identical(.sr_quantile(0.95, 10, 12), qtukey(0.95, 10, 12))
})

test_that("the exact computation agrees with qtukey() where qtukey() is accurate", {
  skip_on_cran()  # forced exact roots for a grid of means
  # At these cells qtukey() was within 1e-5 of the independent quadratures;
  # the exact engine, forced, must agree with it to 1e-5.
  for (k in c(3, 5, 10)) for (nu in c(12, 30, 120)) for (p in c(0.9, 0.95, 0.99)) {
    expect_lt(abs(.sr_quantile_exact_one(p, k, nu) - qtukey(p, k, nu)), 1e-5)
    q <- qtukey(p, k, nu)
    expect_lt(abs(.sr_cdf_exact_one(q, k, nu) - ptukey(q, k, nu)), 1e-6)
  }
})

test_that("a value qtukey() fails to compute comes from the exact root", {
  skip_on_cran()  # a spline of the range distribution for 51 means
  expect_true(is.nan(suppressWarnings(qtukey(0.25, 51, 32))))
  expect_no_warning(q <- .sr_quantile(0.25, 51, 32))
  expect_equal(q, 3.984255805519, tolerance = 1e-8)
  # In a vector, only the cell that failed is replaced; the others keep R's
  # value.
  both <- .sr_quantile(c(0.25, 0.95), 51, 32)
  expect_equal(both[1], 3.984255805519, tolerance = 1e-8)
  expect_identical(both[2], qtukey(0.95, 51, 32))
})

test_that("outside the domain the result is NaN with R's own warning", {
  expect_warning(out <- .sr_quantile(0.95, 3, 0), "NaN")
  expect_true(is.nan(out))
  expect_warning(out <- .sr_cdf(1, 3, -1), "NaN")
  expect_true(is.nan(out))
  expect_true(is.na(.sr_quantile(NA, 3, 10)))
})

test_that("the exact path reproduces the printed Pearson-Hartley table at 5 and 7 df", {
  skip_on_cran()  # a spline of the range distribution per number of means
  # Appendix Table A.4 as printed before the 2027 edition of Maxwell,
  # Delaney, and Kelley, abridged from Table 29 of Pearson and Hartley's
  # Biometrika Tables for Statisticians (3rd ed., Vol. 1). At 5 and 7 error df
  # every k up to 20 is in the exact region. One entry, k = 5 at 7 df and
  # alpha = .01, prints 7.01; its exact value is 7.004993, within 1e-5 of the
  # rounding point, so it is checked to 0.01.
  k <- c(2, 3, 4, 5, 6, 8, 10, 15, 20)
  printed <- rbind(
    "5 .05" = c(3.64, 4.60, 5.22, 5.67, 6.03, 6.58, 6.99, 7.72, 8.21),
    "5 .01" = c(5.70, 6.98, 7.80, 8.42, 8.91, 9.67, 10.24, 11.24, 11.93),
    "7 .05" = c(3.34, 4.16, 4.68, 5.06, 5.36, 5.82, 6.16, 6.76, 7.17),
    "7 .01" = c(4.95, 5.92, 6.54, 7.01, 7.37, 7.94, 8.37, 9.12, 9.65))
  df    <- c(5, 5, 7, 7)
  alpha <- c(.05, .01, .05, .01)
  for (r in seq_len(nrow(printed))) for (j in seq_along(k)) {
    expect_true(.sr_needs_exact(k[j], df[r]))
    v <- .sr_quantile(1 - alpha[r], k[j], df[r])
    if (df[r] == 7 && alpha[r] == .01 && k[j] == 5) {
      expect_equal(v, 7.004993, tolerance = 1e-6)
      expect_lt(abs(v - unname(printed[r, j])), 0.01)
    } else {
      expect_equal(round(v, 2), unname(printed[r, j]),
                   info = paste(rownames(printed)[r], k[j]))
    }
  }
})

test_that("a raw-data simulation confirms the exact quantile at 2 error df", {
  skip_on_cran()  # 200000 simulated studentized ranges
  # Q = (range of 20 iid N(0, 1)) / S, with S^2 the sample variance of 3 more
  # iid N(0, 1) draws (2 df). The upper tail beyond the exact .99 quantile
  # should be .01 within Monte Carlo error; beyond qtukey()'s value it is
  # about .0073.
  R <- 200000L; k <- 20L
  Q <- withr::with_seed(20261111, {
    Z <- matrix(rnorm(R * k), R, k)
    rng <- do.call(pmax, as.data.frame(Z)) - do.call(pmin, as.data.frame(Z))
    E <- matrix(rnorm(R * 3L), R, 3L)
    rng / sqrt(rowSums((E - rowMeans(E))^2) / 2)
  })
  expect_equal(sum(!is.finite(Q)), 0L)  # attempted R, usable R, failed 0
  mcse <- sqrt(0.01 * 0.99 / R)
  expect_lt(abs(mean(Q > .sr_quantile(0.99, 20, 2)) - 0.01), 4 * mcse)
  expect_lt(mean(Q > qtukey(0.99, 20, 2)), 0.01 - 4 * mcse)
})
