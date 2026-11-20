# The noncentral chi square distribution function for every caller that
# inverts it or reads a tail probability from it: ci_nc_chisq() (and so
# ci_rmsea(), ecvi(), and the RMSEA limits of sem_fit_measures() and
# equivalence_rmsea()), and the tests of close and not close fit.
#
# R's pchisq() with a noncentrality parameter sums the Poisson mixture of
# central chi square distributions from its first term upward, and stops
# after 10^6 terms. Its help page warns that noncentralities above about
# 1e5 may give inaccurate results, and above about 2e6 the series is cut
# off before it reaches the bulk of the Poisson weights, so the
# probability it returns is wrong (0 instead of .19 at a statistic and a
# noncentrality of 3e6 on 3000 df) and the warning that says so is easy to
# miss. A noncentrality that large needs a statistic of millions, as from a
# model of many variables fitted to a sample of 10^5 cases, and the RMSEA
# interval then collapsed to a point below its own estimate.
#
# Above 1e5 the mixture is summed here directly, over the Poisson weights
# that carry all but 1e-17 of the mass in each tail:
#   P(X <= q) = sum_j dpois(j, ncp/2) pchisq(q, df + 2j).
# Each term is a central chi square probability, which R computes to full
# precision for any df, so the sum is accurate to rounding: against the
# representation X = (Z + sqrt(ncp))^2 + chi square on df - 1, integrated
# numerically over Z, it agrees to 4e-15 at 1e6 and 2e-14 at 3e6, where
# R's pchisq() is off by 7e-9 and by .93. The window holds about
# 24 sqrt(ncp/2) terms, some 4,000 at 1e5 and 30,000 at 3e6, a few
# milliseconds. At and below 1e5 R's pchisq() is used (off by about 3e-10
# there), so no result in that range changes. The upper tail is summed as
# the upper tail, not as one minus the lower, so that a small p-value keeps
# its relative accuracy.
.pchisq_nc <- function(q, df, ncp, lower.tail = TRUE) {
  if (!is.finite(ncp) || ncp <= 1e5) {
    return(stats::pchisq(q, df = df, ncp = ncp, lower.tail = lower.tail))
  }
  mu <- ncp / 2
  j <- seq(stats::qpois(1e-17, mu), stats::qpois(1e-17, mu, lower.tail = FALSE))
  sum(stats::dpois(j, mu) * stats::pchisq(q, df = df + 2 * j, lower.tail = lower.tail))
}
