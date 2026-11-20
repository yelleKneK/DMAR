# The studentized maximum modulus (SMM) distribution, shared by cv_smm() and
# ci_dunnett_t3(). The SMM variate is M = max_i |Z_i| / S, with Z_1, ..., Z_m
# independent standard normal and S = sqrt(U / df), U ~ chi^2_df, independent
# of them. Conditioning on S makes the m moduli independent, so
#   P(M <= q) = E_S[(2 * Phi(q * S) - 1)^m],
# a one-dimensional integral. Everything here is deterministic (no Monte
# Carlo), and df need not be an integer, which the Welch degrees of freedom of
# ci_dunnett_t3() never are.

# P(M <= q), or P(M > q) when lower_tail = FALSE, for a single q.
#
# The two tails are integrated in different coordinates because each needs
# something different.
#
# The lower tail is the integral over the chi squared variate U that cv_smm()
# has used since 1.0.0, kept exactly so that its critical values are unchanged.
# Its limits are the extreme quantiles of U rather than (0, Inf): for a large
# df the density of U is a narrow bump far from the origin, and an adaptive
# rule asked to cover the whole half-line samples its way past the mass and
# returns zero (at df = 200 it does exactly that).
#
# The upper tail is a p-value, so it needs relative accuracy far into the tail,
# where 1 minus the lower tail has none (the lower tail is accurate to about
# 1e-10 in absolute terms, so 1 minus it cannot resolve a p-value of 1e-12).
# It conditions instead on W = max_i |Z_i|, whose density is
# g(w) = m * (2 * Phi(w) - 1)^(m - 1) * 2 * phi(w):
#   P(M > q) = P(S < W / q) = E_W[F_S(W / q)],  F_S(s) = pchisq(df * s^2, df),
# with abs.tol = 0 so the rule works to rel.tol however small the integral.
# The integrand is the product of an increasing distribution function and a
# light-tailed density, so its mass is a single bump that, for a large q and a
# large df, sits far out where both factors are tiny; split at the bump's mode,
# both pieces have the bump at an endpoint, where the rule always samples.
# Checked against the closed forms (m = 1, df = Inf), against 2 * pt(-q, df) <=
# P(M > q) <= 1 - (1 - 2 * pt(-q, df))^m, and against a trapezoid rule on log S
# in tests/testthat/test-ci_dunnett_t3.R, and live against mvtnorm::pmvt() in
# tools/oracle_snippets/pmcmrplus_dunnett_t3.R.
.smm_cdf <- function(q, df, m, lower_tail = TRUE) {
  if (q <= 0) return(if (lower_tail) 0 else 1)
  if (is.infinite(q)) return(if (lower_tail) 1 else 0)
  # One comparison: the two-sided t (or, at df = Inf, normal) probability.
  if (m == 1) {
    upper <- 2 * stats::pt(-q, df)
    return(if (lower_tail) 1 - upper else upper)
  }
  # Known variance: the m moduli are independent normal moduli.
  if (is.infinite(df)) {
    log_lower <- m * log1p(-2 * stats::pnorm(-q))
    return(if (lower_tail) exp(log_lower) else -expm1(log_lower))
  }
  if (lower_tail) {
    u_lo <- stats::qchisq(1e-12, df)
    u_hi <- stats::qchisq(1e-12, df, lower.tail = FALSE)
    return(stats::integrate(function(u)
      (2 * stats::pnorm(q * sqrt(u / df)) - 1)^m * stats::dchisq(u, df),
      lower = u_lo, upper = u_hi, rel.tol = 1e-10)$value)
  }
  log_integrand <- function(w)
    stats::pchisq(df * (w / q)^2, df, log.p = TRUE) + log(m) +
      (m - 1) * log1p(-2 * stats::pnorm(-w)) + log(2) +
      stats::dnorm(w, log = TRUE)
  integrand <- function(w) exp(log_integrand(w))
  peak <- stats::optimize(log_integrand, c(0, 2 * q + 10),
                          maximum = TRUE)$maximum
  stats::integrate(integrand, 0, peak, rel.tol = 1e-10, abs.tol = 0,
                   subdivisions = 1000L)$value +
    stats::integrate(integrand, peak, Inf, rel.tol = 1e-10, abs.tol = 0,
                     subdivisions = 1000L)$value
}

# The upper-alpha_level SMM critical value: the q with P(M > q) = alpha_level.
#
# The root is found on the lower tail, as cv_smm() always has, and then
# checked against the upper tail. The lower-tail integrand rises from 0 to 1
# near u = df / q^2, and when df is small and q large (df below 1 at
# conventional levels, for example below about .8 at alpha_level = .05 and .6
# at .10, or df of about 2 or less with an alpha_level of .02 or less) that
# rise is a sliver by the origin, where the density of U is largest, which the
# adaptive rule can miss. It then either stops ("the integral is probably
# divergent") or returns a root whose tail area is wrong: at df = 2 and
# alpha_level = .001 the tail area was 1.5 to 5 times alpha_level, depending
# on m. The check accepts the root when its tail area is alpha_level to a
# relative 1e-6 (an accurate root meets it to about 1e-10), which on a grid of
# 600 cells was every case in which cv_smm() had been right, each returned bit
# for bit as before; otherwise it solves on the upper tail, on the log scale
# since alpha_level may be small.
.smm_upper_cv <- function(alpha_level, df, m) {
  if (m == 1)
    return(if (is.infinite(df)) stats::qnorm(1 - alpha_level / 2)
           else                 stats::qt(1 - alpha_level / 2, df))
  if (is.infinite(df))
    return(stats::qnorm((1 + (1 - alpha_level)^(1 / m)) / 2))
  target <- 1 - alpha_level
  # The SMM value is at least the two-sided single-t value (equality at
  # m = 1); use that as a lower bracket and grow the upper end until the CDF
  # clears 1 - alpha_level.
  lo <- stats::qt(1 - alpha_level / 2, df)
  root <- tryCatch({
    cdf <- function(q) .smm_cdf(q, df, m)
    hi <- lo
    while (cdf(hi) < target && hi < 1e4) hi <- hi * 1.5
    stats::uniroot(function(q) cdf(q) - target, c(lo, hi), tol = 1e-10)$root
  }, error = function(e) NA_real_)
  if (!is.na(root) &&
      abs(.smm_cdf(root, df, m, lower_tail = FALSE) / alpha_level - 1) < 1e-6)
    return(root)
  log_excess <- function(q)
    log(.smm_cdf(q, df, m, lower_tail = FALSE)) - log(alpha_level)
  hi <- lo
  while (log_excess(hi) > 0) hi <- hi * 1.5
  stats::uniroot(log_excess, c(lo, hi), tol = 1e-10)$root
}
