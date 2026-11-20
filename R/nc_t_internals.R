# The noncentral t distribution, computed exactly where R's pt(), qt(), and
# dt() are not.
#
# Shared by ci_nc_t() (and so every interval built on it: ci_smd(),
# ci_smd_c(), ci_sm(), ci_sc(), ci_sc_ancova(), ci_cv(), ci_reg_coef(),
# smd_trimmed(), ...), the assurance steps of the AIPE planners built on
# those intervals (ss_aipe_smd(), ss_aipe_sm(), ss_aipe_sc(),
# ss_aipe_sc_ancova(), ss_aipe_cv()), the equivalence planners
# (ss_aipe_equivalence_smd() and its sensitivity function), the critical
# values of cv_t(), the individual intervals of limits_of_agreement(), the
# power of design_consequences(), of the ss_power_* planners whose test is a
# t test (through .power_noncentral_t(), R/ss_power_internals.R), of the
# composite power functions (their marginal powers), and of
# power_equivalence_c() (noninferiority), the component powers of
# ss_power_indirect_effect(), and the likelihood of the Bayes t tests
# (bayes_one_sample_t(), bayes_paired_t(), bayes_independent_t()). No
# other code in R/ calls pt(), qt(), or dt() with a noncentrality.
#
# Fewer than one degree of freedom (0 < df < 1, which cv_t() and ci_nc_t()
# accept) is always integrated, by .nc_t_integral() over u = s^df: pt()
# with a noncentrality is not accurate there (see .nc_t_pt_limit()).
#
# A noncentral t variate on df degrees of freedom with noncentrality ncp is
# T = (Z + ncp) / S, S = sqrt(V / df), with Z standard normal and V chi
# square on df degrees of freedom, independent. Conditioning on V,
#   P(T <= q) = integral over v of Phi(q s - ncp) dchisq(v, df) dv,
#   P(T >  q) = integral over v of Phi(ncp - q s) dchisq(v, df) dv,
#   P(|T| <= k) = integral over v of
#                   {Phi(k s - ncp) - Phi(-k s - ncp)} dchisq(v, df) dv,
# with s = sqrt(v / df). Only the normal distribution function and the chi
# square density enter. Each tail is integrated directly, so a small tail
# probability keeps its relative accuracy. The integrals run between the
# 1e-14 and 1 - 1e-14 quantiles of the chi square, which leaves out at
# most 2e-14 of probability, with a relative tolerance of 1e-11 and an
# absolute one of 1e-13 (see .nc_t_integral()).
#
# R's pt() with ncp sums the series of Lenth (1989, AS 243) up to a
# noncentrality of 37.62 (del^2 > 2 log(2) 1021, R's pnt.c) and with at
# most 400,000 degrees of freedom, and otherwise returns the normal
# approximation of Abramowitz and Stegun (1964, 26.7.10). The approximation
# is poor: over the quantiles of the grid below, at noncentralities from
# 37.7 to 100, it was off by as much as .14 with 1 df, .0076 with 38, .0011
# with 1,000, 9e-5 with 1e4, and 1e-7 with 4e5; ci_smd() at n = 20 per
# group, d = 10, and 99 percent gave a width of 5.996 where the exact width
# is 6.096. The
# series itself stops short of convergence near the top of its range when
# the degrees of freedom are large. The range in which pt() is used, and
# the evidence for it, is stated at .nc_t_pt_accurate().

# .nc_t_integral(g, df, breaks = NULL, breaks_minus_1 = breaks - 1)
#
# The expectation of g(S, d) for S = sqrt(V / df), V chi square on df
# degrees of freedom, by stats::integrate() over s with the density of S,
# 2 df s dchisq(df s^2, df). With more than 1e7 degrees of freedom d is
# S - 1, computed without rounding S (see the integral over w below), and
# with fewer it is NULL; g forms each normal argument q s - ncp from s and
# d by .nc_t_affine(). Working in s rather than v removes the
# infinite density of the chi square at 0 when df < 2 (the density of S is
# finite at 0 for df >= 1). `breaks`, in increasing order, are points in s
# where g changes quickly, and `breaks_minus_1` the same points less 1,
# computed by the caller without rounding them near 1 (the integral over w
# places its breaks by them); the range is split there, so that a step in g
# that is narrow beside the spread of S (a large observed t with few
# degrees of freedom) lies at the end of a piece and is not missed. With
# the integral over v in one piece, the form ss_aipe_smd() used before,
# P(T > 1e6) at df = 1 and ncp = 3e4 came out 0 (it is .0239) and
# P(T <= -318310) at df = 1 and ncp = 0 came out 0 (it is 1e-6).
#
# With 0 < df < 1 the density of S is itself infinite at 0, like s^(df - 1),
# and a step of g near 0 (a quantile far in a tail, where q s - ncp crosses
# 0 at a small s) sat beside that singularity: integrate() stopped with
# "the integral is probably divergent" at df = 0.5, ncp = 38, and the
# .9995 quantile, and at df = 0.2, ncp = 45, and the .975 quantile. There
# the integral is taken over u = s^df instead. Since
# du = df s^(df - 1) ds, the density of S times ds is
#   c exp(-df s^2 / 2) du,  c = 2 df^(df / 2 - 1) / {2^(df / 2) Gamma(df / 2)},
# which is finite and smooth on [0, u_max]; the range runs from u = 0 (s =
# 0, so none of the lower tail of S is left out) to the 1 - 1e-14 quantile,
# and the breaks are mapped to u. Against the Poisson mixture of incomplete
# beta functions (Guenther, 1978) and, where its terms cancel, the
# Hermite-function integral of the e2 review (work/round6/evidence/
# g1_final_fixes/04): at 462 quantiles (df from .05 to .99, ncp from -45 to
# 45, probabilities from 1e-6 to 1 - 1e-6) the tail at the quantile was
# within 2.7e-11 of its target (relative) for probabilities from .0005 to
# .9995 and 1e-8 at 1e-6, every probability at 378 fixed points within
# 3.6e-13, and the density within 1.2e-12 of its largest value. Should
# integrate() still fail, the failure is refused with a condition of class
# dmar_computation_failed that names df, never passed on as integrate()'s
# own error.
#
# With more than 1e7 degrees of freedom, S is within about 8.5 sqrt(1 /
# (2 df)) of 1, and neither s nor v = df s^2 can be held finely enough in
# a double for the density to be smooth at the scale of that spread: with
# 6.1e15 degrees of freedom, df s^2 moves in steps of 1.35 where the chi
# square has a standard deviation of 1.1e8, integrate() stopped with
# "roundoff error was detected", and from 1e11 on the integral over s was
# off by as much as 1.6e-9. There the integral is taken over the
# standardized chi square, w = (v - df) / sqrt(2 df), so s = sqrt(1 + h)
# with h = w sqrt(2 / df). With k = df / 2, the density of w is exactly
#   exp{-stirlerr(k) + k (log1p(h) - h) - log1p(h)} / sqrt(2 pi),
# where stirlerr(k) = log Gamma(k + 1) - log{sqrt(2 pi k) (k / e)^k} =
# 1 / (12 k) - 1 / (360 k^3) + ... and k (log1p(h) - h) = w^2 (-1/2 + h / 3
# - h^2 / 4 + ...), since k h^2 = w^2. Neither v nor the difference
# log1p(h) - h is formed, so the density is smooth in w. With df > 1e7 and
# |w| < 9, |h| < .004, and the terms of the series past h^8 and of
# stirlerr() past k^-3 are below 1e-20.
#
# Nor is s formed for g: g is given d = s - 1 = h / (1 + sqrt(1 + h)),
# which is smooth in w, and forms q s - ncp as (q - ncp) + q d. A rounded
# s = sqrt(1 + h) is held only to about 1.1e-16, so q s - ncp formed from it
# moved in steps of q times that: where q is large beside the spread of S
# (q sd(S) of .3 or more at 1e14 degrees of freedom) integrate() stopped
# with "roundoff error was detected", and where it returned, P(T <= 1.4e7)
# at df = 1e14 and ncp = 1.4e7 + 1 was off by 1.2e-10. With fewer degrees
# of freedom s is the variable of integration (or a power of it) and is
# exact, and q s - ncp is formed directly: (q - ncp) + q (s - 1) cancels
# where s is far from 1, as it is with few degrees of freedom and a large
# q, and it made integrate() fail at df = 0.05 and 0.2 with q = 1e9.
#
# The range in w runs between the 1e-14 and 1 - 1e-14 quantiles of the chi
# square by the approximation of Wilson and Hilferty (1931), written in w
# without cancellation. It is within 1e-6 of the quantile of qchisq() from
# 1e7 to 1e20 degrees of freedom, so it leaves out the same probability,
# and it keeps its digits as df grows, where qchisq(p, df) - df loses them
# (at 4.1e31 degrees of freedom that put the lower end of the range .31 of
# a standard deviation of w from the quantile). The breaks are mapped to w
# from breaks_minus_1, as (s - 1) (s + 1) / sqrt(2 / df).
#
# The evidence (work/round6/evidence/m4_nct, wgrid_compare.log, against a
# 256-bit oracle that integrates over w with the chi square density from
# the log gamma function in Rmpfr and forms q s - ncp in 256 bits): at 13
# df from 1.0000001e7 to 1e32, q sd(S) from 1e-3 to 10, and q - ncp from
# -4 to 4 times sqrt(1 + (q sd(S))^2), the spread of T - ncp (1,820
# values), both tails of the distribution function, P(|T| <= q), and the
# density at q (in units of its largest value) were within 2.0e-14 of the
# oracle, the probability the range of the integral leaves out, and
# integrate() did not fail. With a rounded s and the range from qchisq()
# (the code before this change), 595 of the 1,820 failed, the error
# reached 2.1e-10 at 1e14 df and 3.1e-8 at 1e18, and from 4.1e31 on the
# probabilities were off by .17, all from the rounded s: with the range
# above and a rounded s they are as far off (fix3/range_vs_s.log). The range
# covers the planners: a plan of ss_aipe_sc() or ss_aipe_sc_ancova() up to
# 2^53 per group with J groups has at most J (2^53 - 1) degrees of
# freedom, and a vector of weights in R has at most 2^52 elements (R
# refuses only longer ones), so at most 2^52 (2^53 - 1), fewer than 4.1e31.
.nc_t_integral <- function(g, df, breaks = NULL, breaks_minus_1 = breaks - 1) {
  if (df > 1e7) {
    k <- df / 2
    r <- sqrt(2 / df)
    log_c <- -0.5 * log(2 * pi) - (1 / (12 * k) - 1 / (360 * k^3))
    integrand <- function(x) {
      h <- x * r
      series <- -1 / 2 + h * (1 / 3 + h * (-1 / 4 + h * (1 / 5 + h * (
        -1 / 6 + h * (1 / 7 + h * (-1 / 8 + h * (1 / 9 - h / 10)))))))
      d <- h / (1 + sqrt(1 + h))
      g(1 + d, d) * exp(log_c + x^2 * series - log1p(h))
    }
    a <- 2 / (9 * df)
    z <- stats::qnorm(1e-14)
    lower <- df * expm1(3 * log1p(z * sqrt(a) - a)) / sqrt(2 * df)
    upper <- df * expm1(3 * log1p(-z * sqrt(a) - a)) / sqrt(2 * df)
    d <- breaks_minus_1[breaks_minus_1 > -1]
    breaks <- d * (2 + d) / r
    points <- c(lower, sort(breaks[breaks > lower & breaks < upper]), upper)
    return(.nc_t_integrate_pieces(integrand, points, df))
  }
  upper <- sqrt(stats::qchisq(1e-14, df, lower.tail = FALSE) / df)
  if (df < 1) {
    log_c <- log(2) + (df / 2 - 1) * log(df) - (df / 2) * log(2) -
      lgamma(df / 2)
    integrand <- function(x) {
      s <- x^(1 / df)
      g(s, NULL) * exp(log_c - df * s^2 / 2)
    }
    breaks <- breaks[breaks > 0 & breaks < upper]
    points <- c(0, sort(breaks^df), upper^df)
  } else {
    lower <- sqrt(stats::qchisq(1e-14, df) / df)
    integrand <- function(x)
      g(x, NULL) * 2 * df * x * stats::dchisq(df * x^2, df)
    points <- c(lower, breaks[breaks > lower & breaks < upper], upper)
  }
  .nc_t_integrate_pieces(integrand, points, df)
}

# The integral of `integrand` over the pieces between successive `points`,
# with the tolerances above; a failure of integrate() is refused with a
# condition of class dmar_computation_failed.
.nc_t_integrate_pieces <- function(integrand, points, df) {
  total <- 0
  for (i in seq_len(length(points) - 1L)) {
    piece <- tryCatch(stats::integrate(
      integrand, points[i], points[i + 1L],
      rel.tol = 1e-11, abs.tol = 1e-13, subdivisions = 2000L)$value,
      error = function(e) .dmar_refuse(c(
        "The noncentral t distribution with df = ", format(df, digits = 7),
        " could not be computed reliably here: the numerical integral over ",
        "the sampling distribution of the standard deviation failed (",
        conditionMessage(e), ")."), "computation_failed", df = df))
    total <- total + piece
  }
  total
}

# .nc_t_affine(a, b, s, d)
#
# a s - b, for the integrand of .nc_t_integral(): (a - b) + a d where d =
# s - 1 is given (more than 1e7 degrees of freedom), a s - b where it is
# NULL.
.nc_t_affine <- function(a, b, s, d) if (is.null(d)) a * s - b else
  (a - b) + a * d

# .nc_t_cdf(q, df, ncp, lower_tail = TRUE)
#
# The exact noncentral t probability P(T <= q) (lower_tail = TRUE) or
# P(T > q), by the integral over S. Phi(q s - ncp) moves from 0 to 1 where
# q s - ncp runs from -8 to 8, so the range is split at
# s = (ncp - 8) / q and (ncp + 8) / q. With more than 1e7 degrees of
# freedom the argument is formed as (q - ncp) + q (s - 1) (see
# .nc_t_integral()). With df = Inf, T is normal with mean ncp and variance
# 1. An infinite q gives 0 or 1, as pt() does.
.nc_t_cdf <- function(q, df, ncp, lower_tail = TRUE) {
  if (is.infinite(df))
    return(stats::pnorm(q, mean = ncp, lower.tail = lower_tail))
  if (is.infinite(q)) return(as.numeric(xor(q < 0, lower_tail)))
  sign <- if (lower_tail) 1 else -1
  shift <- q - ncp
  edges <- if (q > 0) c(-8, 8) else if (q < 0) c(8, -8)
  value <- .nc_t_integral(
    function(s, d) stats::pnorm(sign * .nc_t_affine(q, ncp, s, d)), df,
    (ncp + edges) / q, (edges - shift) / q)
  min(max(value, 0), 1)
}

# .nc_t_abs_probability(k, df, ncp)
#
# P(|T| <= k), by the integral over S, split where either normal
# probability moves (see .nc_t_cdf()); with more than 1e7 degrees of
# freedom the arguments k s - ncp and -k s - ncp are formed as
# (k - ncp) + k (s - 1) and (-k - ncp) - k (s - 1). It is the same for ncp
# and -ncp, and |ncp| is used.
.nc_t_abs_probability <- function(k, df, ncp) {
  if (k <= 0) return(0)
  ncp <- abs(ncp)
  if (is.infinite(df))
    return(min(max(stats::pnorm(k - ncp) - stats::pnorm(-k - ncp), 0), 1))
  edges <- if (ncp <= 8) c(ncp - 8, 8 - ncp, ncp + 8) else
    c(8 - ncp, ncp - 8, ncp + 8)
  value <- .nc_t_integral(
    function(s, d) stats::pnorm(.nc_t_affine(k, ncp, s, d)) -
      stats::pnorm(.nc_t_affine(-k, ncp, s, d)), df,
    edges / k, (edges - k) / k)
  min(max(value, 0), 1)
}

# .nc_t_pt_limit(df) and .nc_t_pt_accurate(df, ncp)
#
# The largest |ncp| at which R's pt() gives the noncentral t distribution
# with df degrees of freedom to within 1e-9 of the integral (.nc_t_cdf()),
# and whether every noncentrality in `ncp` is within it: 37.62 with 1 to
# 2,000 degrees of freedom, 33 with more up to 400,000, and none beyond
# (there pt() is always the normal approximation; the limit is -Inf) or
# below 1. With 0 < df < 1 the series pt() sums is not accurate inside
# 37.62 either, and the error grows quickly as df falls: at the exact
# quantiles of probability .0005 to .9995 (|ncp| <= 37.5), the tail pt()
# gave was off by as much as .14 at df = 0.05, .070 at 0.1, .010 at 0.2,
# .0013 at 0.3, 2.4e-4 at 0.4, 3.8e-5 at 0.5, 2.8e-8 at 0.6, and 2.1e-11
# at 0.8, and the large errors are not confined to the far tails (.14 at
# the .8 quantile with df = 0.05; work/round6/evidence/g1_final_fixes/
# fix1/23, against an oracle that conditions on the normal variate, and
# review1/16d of that lane). So a nonzero noncentrality with fewer than one
# degree of freedom is always integrated (the integral over u = s^df of
# .nc_t_integral()). The grid below starts at df = 1.
#
# The evidence (work/round6/evidence/e1_nct_exact, run on R 4.6.1, arm64,
# where R's long double is a double; an x86 build sums the series in
# extended precision and is no less accurate). Both tails of pt() were
# compared with the integral at the exact quantiles of probability 1e-6,
# 1e-4, .001, .005, .025, .1, .5, .9, .975, .995, .999, .9999, and 1 - 1e-6
# of each distribution, over df = 1, 2, 3, 5, 10, 20, 50, 100, 200, 500,
# 1000, 1500, 2000, 2500, 3000, 4000, 5000, 7500, 1e4, 2e4, 5e4, 1e5, 2e5,
# and 4e5 and ncp = 0 to 28 by 2, 30 to 33 by .5, 33.2 to 37.6 by .2,
# 37.62, and -37.62, -36, -34, -30, and -20 (15,912 points), and then more
# finely near the edge (df from 2000 to 4e5, ncp by .05 or .06, quantiles
# out to 1e-7 and 1 - 1e-7). With at most 3,000 degrees of freedom the
# error was at most about 1e-12 up to 37.62, except with one degree of
# freedom, where it reached 2.3e-10 at the quantile 1 - 1e-6 (q near 2e7)
# and 2e-9 at 1 - 1e-7. With more, the series stops short of convergence
# in the upper tail near the top of its range: the largest noncentrality
# with every error below 1e-9 was 37.38 at df = 3250, 36.54 at 4000, 35.0
# at 5000, 34.55 at 1e4, 33.9 at 2e4, 33.55 at 5e4, 33.4 at 1e5, and 33.3
# at 4e5, and the error reached 1e-6 at 34, 1e-3 at 36, and .1 at 37.6
# with 1e5 degrees of freedom. Inside the rule the largest error found
# was 2e-9 (at df = 1 and the quantile 1 - 1e-7) and otherwise 1.6e-10.
# pt() warned (full precision may not have been achieved in 'pnt{final}')
# at 414 points, and never in the tail .nc_t_pt() asks for. The integral
# agreed with the oracle that conditions on the normal variate to 8.1e-13
# at every seventh point (2,274), and the two tails it gives summed to 1
# within 8.1e-13.
.nc_t_pt_limit <- function(df) {
  if (df < 1 || df > 4e5) -Inf else if (df <= 2000) 37.62 else 33
}

# With ncp = 0, pt() is the central t at any df (pnt() returns pt() there),
# which is accurate, so a zero noncentrality is always inside the range.
.nc_t_pt_accurate <- function(df, ncp)
  all(ncp == 0 | abs(ncp) <= .nc_t_pt_limit(df))

# .nc_t_pt(q, df, ncp, lower_tail = TRUE)
#
# pt(q, df, ncp, lower.tail = lower_tail), asked for in the tail that does
# not warn. pt() sums P(T <= q) directly when q >= 0 (and P(T > q), as the
# lower tail of -T, when q < 0) and warns that full precision may not have
# been achieved when that sum is within 1e-10 of 1 and is the tail asked
# for. Asking for the other tail and subtracting from 1 returns the same
# number to within rounding (pt() itself returns 1 minus the sum for it)
# and does not warn.
.nc_t_pt <- function(q, df, ncp, lower_tail = TRUE) {
  if (q >= 0) {
    upper <- stats::pt(q, df, ncp = ncp, lower.tail = FALSE)
    if (lower_tail) 1 - upper else upper
  } else {
    lower <- stats::pt(q, df, ncp = ncp)
    if (lower_tail) lower else 1 - lower
  }
}

# .nc_t_prob(q, df, ncp, lower_tail = TRUE)
#
# P(T <= q) or P(T > q): pt() where it is accurate (.nc_t_pt_accurate()),
# which is fast, and the integral (.nc_t_cdf()) elsewhere.
.nc_t_prob <- function(q, df, ncp, lower_tail = TRUE) {
  if (.nc_t_pt_accurate(df, ncp)) .nc_t_pt(q, df, ncp, lower_tail)
  else .nc_t_cdf(q, df, ncp, lower_tail)
}

# .nc_t_quantile(p, df, ncp)
#
# The q with P(T <= q) = p, by uniroot() on .nc_t_prob(), so it is exact
# where R's qt() with ncp, which inverts pt(), is not, and it does not warn
# as qt() does when its search meets a probability near 1. Each tail is
# matched directly, so a quantile far in either tail keeps the accuracy of
# the probability it inverts. p = 0 and 1 give -Inf and Inf, ncp = 0 the
# central quantile qt(p, df), and df = Inf the normal quantile.
#
# The tolerance, 1e-12 times max(1, |ncp|), puts the root where pt() puts
# it. The evidence is in work/round6/evidence/e2_nct_users. Inside the range
# of .nc_t_pt_accurate(), at the probabilities .001 to .999, the quantile
# agreed with qt() to 2.5e-12 of its size with 1 to 1,000 df, 1e-11 up to
# 5,000, and 1.6e-9 at 4e5 (05), and with the series oracle as closely as
# qt() did (2.7e-9, at df = 1e5, where pt() itself is no closer; 06).
# Outside the range it agreed with the series oracle to 1.4e-10 at those
# probabilities and 1e-8 at 1e-6 and 1 - 1e-6, and the conditioning oracle
# put the tail at the quantile within 5.1e-11 of its target (relative);
# qt() gave 429.18 there for the .975 quantile with ncp = 50 and df = 3,
# which is 186.49, and an infinite quantile at df = 1 with |ncp| >= 38.
# At 1e-6 and 1 - 1e-6 inside the range the quantile, like qt(), carries
# pt()'s absolute error, which there is a larger share of the tail: 1.8e-4
# of it at df = 1 (q near 3e7), and at df = 4e5 and ncp = 10,
# P(T <= 5.246472) is 9.99875e-7 where pt() gives 9.99995e-7.
#
# A quantile past the largest double is -Inf or Inf, as qt() returns it.
# With few degrees of freedom the tails are so heavy that this is common:
# with df = 0.001, P(T > 1.8e308) is about .49 whatever the noncentrality,
# and qt(5e-7, 0.01) is -Inf. uniroot()'s extension of the bracket doubles
# its step for at most 1,000 steps, which reaches about 1e301 times the
# first step, and then stops with the unclassed error "no sign change
# found in 1000 iterations" (cv_t(.05, df = .001, ncp = 5), and 64 of 468
# calls with df from .001 to .999999 in work/round6/evidence/
# g1_final_fixes/review1/06). That error is caught here and the tails at
# -/+ .Machine$double.xmax are computed: if P(T <= -xmax) is at least the
# lower target the quantile is -Inf, if P(T > xmax) is at least the upper
# one it is Inf, and otherwise the root is finite, beyond where the
# extension stopped, and is found by uniroot() between the first bracket
# and -/+ xmax. Any other failure of the search is refused with a
# condition of class dmar_computation_failed that names df; a refusal of
# .nc_t_integral() is passed on as it is.
.nc_t_quantile <- function(p, df, ncp) {
  if (p <= 0) return(-Inf)
  if (p >= 1) return(Inf)
  if (ncp == 0) return(stats::qt(p, df))
  if (is.infinite(df)) return(stats::qnorm(p, mean = ncp))
  spread <- sqrt(1 + ncp^2 / (2 * df))
  f <- if (p <= .5) function(q) .nc_t_prob(q, df, ncp) - p else
    function(q) (1 - p) - .nc_t_prob(q, df, ncp, lower_tail = FALSE)
  first <- ncp + c(-3, 3) * spread
  tol <- 1e-12 * max(1, abs(ncp))
  root <- tryCatch(
    stats::uniroot(f, interval = first, extendInt = "upX", tol = tol)$root,
    error = function(e) e)
  if (!inherits(root, "error")) return(root)
  if (inherits(root, "dmar_refusal")) stop(root)
  .nc_t_quantile_far(f, p, df, ncp, first, tol, root)
}

# .nc_t_quantile_far(f, p, df, ncp, first, tol, failure)
#
# The quantile when uniroot()'s extension of the bracket `first` found no
# sign change of the increasing function f (see .nc_t_quantile()): -Inf or
# Inf when the quantile is past the largest double, the root by uniroot()
# between `first` and -/+ xmax when it is not, and a dmar_computation_failed
# refusal naming df if that search fails too or f is not as it must be.
.nc_t_quantile_far <- function(f, p, df, ncp, first, tol, failure) {
  big <- .Machine$double.xmax
  refuse <- function(why) .dmar_refuse(c(
    "The ", format(p, digits = 7), " quantile of the noncentral t ",
    "distribution with df = ", format(df, digits = 7), " and ncp = ",
    format(ncp, digits = 7), " could not be computed reliably: ", why, "."),
    "computation_failed", df = df)
  f_low <- f(-big)
  f_high <- f(big)
  if (is.na(f_low) || is.na(f_high))
    refuse(c("its tail probability past the largest double could not be ",
             "computed"))
  if (f_low >= 0) return(-Inf)
  if (f_high <= 0) return(Inf)
  f_first <- c(f(first[1]), f(first[2]))
  if (anyNA(f_first) || (f_first[1] <= 0 && f_first[2] >= 0))
    refuse(c("the search for it stopped (", conditionMessage(failure), ")"))
  below <- f_first[1] > 0
  root <- tryCatch(
    stats::uniroot(f, interval = if (below) c(-big, first[1]) else
                     c(first[2], big),
                   f.lower = if (below) f_low else f_first[2],
                   f.upper = if (below) f_first[1] else f_high,
                   tol = tol, maxiter = 5000L)$root,
    error = function(e) e)
  if (inherits(root, "dmar_refusal")) stop(root)
  if (inherits(root, "error"))
    refuse(c("the search for it stopped (", conditionMessage(root), ")"))
  root
}

# .nc_t_power(crit, df, ncp, two_sided = TRUE)
#
# The power of a t test with critical value crit > 0 whose statistic is
# noncentral t: P(T > crit) for a test in the direction of the effect
# (two_sided = FALSE) and P(T > crit) + P(T < -crit) for a nondirectional
# one, with |ncp| as the noncentrality (the power of a test in the
# direction of the effect, or of a two-sided test, does not depend on the
# sign of the effect). Each tail is a probability of .nc_t_prob(), asked for
# directly rather than as 1 minus the other tail, so the upper tail is not
# rounded away when it is small.
.nc_t_power <- function(crit, df, ncp, two_sided = TRUE) {
  ncp <- abs(ncp)
  upper <- .nc_t_prob(crit, df, ncp, lower_tail = FALSE)
  if (!two_sided) return(upper)
  upper + .nc_t_prob(-crit, df, ncp)
}

# The density of the noncentral t ------------------------------------------
#
# Given S = s, T = (Z + ncp) / s has density s phi(x s - ncp), so the
# density of T is
#   f(x) = E[S phi(x S - ncp)] = integral over s of s phi(x s - ncp) g(s) ds,
# g the density of S (.nc_t_integral()). At x = 0 it is
# phi(ncp) E(S) = exp(-ncp^2 / 2) Gamma((df + 1) / 2) /
# (Gamma(df / 2) sqrt(pi df)), the form R's dnt() uses there.
#
# R's dt() with ncp (dnt()) computes
#   df / |x| * |pnt(x sqrt((df + 2) / df), df + 2, ncp) - pnt(x, df, ncp)|,
# a difference of two distribution functions that cancels, so the absolute
# error of pnt() is multiplied by df / |x|. It is therefore less accurate
# than pt(), and over a different region: near x = 0 at every df and at
# large df everywhere, inside the range in which pt() is accurate as well
# as outside it. Against the Poisson mixture of beta densities (Guenther,
# 1978), with the error divided by the largest density of the
# distribution (the analog of the absolute error used for pt()), over 25
# df from 1 to 4e5, 20 ncp inside pt()'s range, and x at 13 quantiles and
# at +/-1e-8 to 1 (02, 15,035 points; 03, 53,436 points with x from .01 to
# the 1e-7 and 1 - 1e-7 quantiles): with df <= 1,000 and |x| >= .1 the
# largest error was 1.1e-9 (at df = 1,000 and x = -.1), and 6e-10 with
# df <= 700; with smaller |x| it reached 5.5e-9 at |x| = .01, 2.2e-8 at
# .001, 4.1e-8 at 1e-4, and 1.8e-6 at 1e-6 (df = 200); with |x| >= .1 and
# more than 1,000 df it was 1.6e-9 at 1,500, 3.7e-9 at 2,000, 2.2e-8 at
# 5,000, 9.8e-8 at 1e4, 7.2e-6 at 1e5, and 1.2e-4 at 4e5, and at 4e5 it
# was 1.3e-2 near x = 0. dt() is therefore used only where
# df <= 1,000, |ncp| <= 37.62, and |x| >= .1 (.nc_t_dt_accurate()), and the
# integral elsewhere; ncp = 0 is the central density, dt(x, df), which is
# accurate everywhere.

# .nc_t_density_integral(x, df, ncp)
#
# The density by the integral over S, split where phi(x s - ncp) is more
# than 8 standard deviations from its peak (s = (ncp -/+ 8) / x), so that a
# peak narrow beside the spread of S is not missed. Over df from 0.5 to
# 1e6, ncp from -60 to 300, and x at nine quantiles and near 0 (07, 2,112
# points), it was within 3e-12 of the largest density of its distribution
# of the oracle that conditions on Z, at the 2,109 points that oracle
# reached, and within 5.6e-10 of the series, which is itself that far from
# the conditioning oracle at df = 1e6. The integral carries the absolute
# tolerance of .nc_t_integral(), so a density far below the largest one is
# accurate in absolute rather than relative terms: one 1e-10 of the
# largest was within 8.5e-5 of itself. R's dt() was off by as much as .15
# of the largest density at those points outside .nc_t_dt_accurate().
.nc_t_density_integral <- function(x, df, ncp) {
  if (is.infinite(df)) return(stats::dnorm(x, mean = ncp))
  if (!is.finite(x)) return(0)
  if (x == 0) {
    # Gamma((df + 1) / 2) / {Gamma(df / 2) sqrt(pi df)} = 1 / {B(df / 2,
    # 1/2) sqrt(df)}. R's lbeta() takes the large argument's log gamma
    # apart (Stirling's series and its correction term), so the log of the
    # ratio does not come from differencing two large log gammas, which
    # lost 3.8e-11 of the density at 1e5 degrees of freedom and 9.7e-9 at
    # 9.99e6. Against 256-bit values (work/round6/evidence/m4_nct,
    # dens0_compare.log) it was within 5.6e-16 below 1 degree of freedom,
    # 2.7e-15 from 1 to 1e16, and 6e-15 from 1e16 to 1e32.
    return(exp(-lbeta(df / 2, 0.5) - 0.5 * log(df) - ncp^2 / 2))
  }
  edges <- if (x > 0) c(-8, 8) else c(8, -8)
  value <- .nc_t_integral(
    function(s, d) s * stats::dnorm(.nc_t_affine(x, ncp, s, d)), df,
    (ncp + edges) / x, (edges - (x - ncp)) / x)
  max(value, 0)
}

# .nc_t_dt_accurate(x, df, ncp)
#
# Elementwise, whether R's dt() is accurate at (x, df, ncp): ncp = 0 (the
# central density), or 1 <= df <= 1,000, |ncp| <= 37.62, and |x| >= .1 (see
# the evidence above, whose grid starts at df = 1; dt() differences two
# values of pt(), which is not accurate below one degree of freedom).
.nc_t_dt_accurate <- function(x, df, ncp) {
  ncp == 0 | (df >= 1 & df <= 1000 & abs(ncp) <= 37.62 & abs(x) >= .1)
}

# .nc_t_dt(x, df, ncp)
#
# dt(x, df, ncp), vectorized over x and ncp, computed so that it does not
# warn. dnt() calls pnt() for the lower tail at x and at x sqrt((df + 2) /
# df), and pnt() warns that full precision may not have been achieved when
# that sum is within 1e-10 of 1, which happens for x > 0 in the upper tail
# (for x < 0 pnt() sums the other tail and does not warn). There the upper
# tails are asked for instead: pt(lower.tail = FALSE) returns 1 minus the
# sum, which is exact when the sum is at least 1/2, so 1 minus it gives
# back the sum dnt() would have used, bit for bit; a sum below 1/2 is far
# from 1 and is asked for directly. The difference is then put together as
# dnt() does, log(df) - log|x| + log|difference|, and exponentiated.
.nc_t_dt <- function(x, df, ncp) {
  n <- max(length(x), length(ncp))
  x <- rep_len(x, n)
  ncp <- rep_len(ncp, n)
  out <- numeric(n)
  # Where dnt() does not call pnt() for an upper-tail sum: x <= 0, x at most
  # sqrt(df DBL_EPSILON) (its closed form), and a central density.
  plain <- x <= sqrt(df * .Machine$double.eps) | ncp == 0 | !is.finite(x)
  if (any(plain)) out[plain] <- stats::dt(x[plain], df, ncp = ncp[plain])
  i <- which(!plain)
  if (length(i)) {
    xi <- x[i]
    nci <- ncp[i]
    x_2 <- xi * sqrt((df + 2) / df)
    lower_sum <- function(q, d, nc) {
      upper <- stats::pt(q, d, ncp = nc, lower.tail = FALSE)
      small <- upper > 0.5
      value <- 1 - upper
      if (any(small)) value[small] <- stats::pt(q[small], d, ncp = nc[small])
      value
    }
    difference <- lower_sum(x_2, df + 2, nci) - lower_sum(xi, df, nci)
    out[i] <- exp(log(df) - log(abs(xi)) + log(abs(difference)))
  }
  out
}

# .nc_t_density(x, df, ncp)
#
# The density of the noncentral t, vectorized over x and ncp (df a single
# value): dt() where .nc_t_dt_accurate() holds, which is fast, and the
# integral elsewhere.
.nc_t_density <- function(x, df, ncp) {
  n <- max(length(x), length(ncp))
  x <- rep_len(x, n)
  ncp <- rep_len(ncp, n)
  out <- numeric(n)
  fast <- .nc_t_dt_accurate(x, df, ncp)
  if (any(fast)) out[fast] <- .nc_t_dt(x[fast], df, ncp[fast])
  for (i in which(!fast)) out[i] <- .nc_t_density_integral(x[i], df, ncp[i])
  out
}
