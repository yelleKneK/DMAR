# Oracles for the noncentral t, coded apart from the package
# (R/nc_t_internals.R integrates over S = sqrt(V / df)); shared by
# test-ci_nc_t.R and test-nc_t_internals.R.
#
# The distribution function conditions on the normal variate Z instead and
# takes the chi square probabilities from pchisq():
#   q > 0: P(T <= q) = Phi(-ncp) + E_Z[ P(V >= df ((Z + ncp) / q)^2); Z > -ncp ]
#   q < 0: P(T <= q) = E_Z[ P(V <= df ((Z + ncp) / q)^2); Z < -ncp ]
# with the other tail written the same way (no subtraction from 1).
.oracle_nct <- function(q, df, ncp, lower_tail = TRUE) {
  if (q == 0) return(stats::pnorm(-ncp, lower.tail = lower_tail))
  if (q > 0) {
    a <- max(-ncp, -40); b <- 40
    base <- if (lower_tail) stats::pnorm(-ncp) else 0
    upper_chisq <- lower_tail
  } else {
    a <- -40; b <- min(-ncp, 40)
    base <- if (lower_tail) 0 else stats::pnorm(ncp)
    upper_chisq <- !lower_tail
  }
  if (a >= b) return(base)
  g <- function(z) stats::dnorm(z) *
    stats::pchisq(df * ((z + ncp) / q)^2, df, lower.tail = !upper_chisq)
  spread <- abs(q) / sqrt(df / 2 + 1)
  edges <- sort(unique(pmin(pmax(c(a, b, -8, 0, 8, q - ncp, -q - ncp,
                                   q - ncp + c(-3, 3) * spread,
                                   -q - ncp + c(-3, 3) * spread), a), b)))
  total <- 0
  for (i in seq_len(length(edges) - 1L)) {
    if (edges[i + 1L] > edges[i])
      total <- total + stats::integrate(g, edges[i], edges[i + 1L],
                                        rel.tol = 1e-13, abs.tol = 1e-17,
                                        subdivisions = 5000L)$value
  }
  base + total
}

# The tail probabilities at the limits ci_nc_t() returns, by the oracle.
.oracle_tails <- function(res, ncp, df) {
  c(lower = .oracle_nct(ncp, df, res$value[1], lower_tail = FALSE),
    upper = .oracle_nct(ncp, df, res$value[2], lower_tail = TRUE))
}

# The density, by the same conditioning: given Z = z, T = (z + ncp) / S, so
# for x > 0 (z > -ncp) the density is g((z + ncp) / x) (z + ncp) / x^2, g
# the density of S, 2 df s dchisq(df s^2, df), and for x < 0 the same over
# z < -ncp with |z + ncp|; at x = 0 it is phi(ncp) E(S).
.oracle_nct_density <- function(x, df, ncp) {
  if (x == 0)
    return(exp(lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * log(pi * df) -
                 ncp^2 / 2))
  g <- function(s) ifelse(s > 0, 2 * df * s * stats::dchisq(df * s^2, df), 0)
  f <- function(z) stats::dnorm(z) * g((z + ncp) / x) * abs(z + ncp) / x^2
  a <- if (x > 0) max(-ncp, -40 - abs(ncp)) else -40 - abs(ncp)
  b <- if (x > 0) 40 + abs(ncp) else min(-ncp, 40 + abs(ncp))
  # S has its mass where (z + ncp) / x is near 1, i.e. z near x - ncp, over
  # a few multiples of |x| / sqrt(2 df).
  spread <- abs(x) * min(1, 1 / sqrt(2 * df))
  edges <- c(a, b, -8, 0, 8, x - ncp + c(-30, -8, -3, -1, 0, 1, 3, 8, 30) * spread)
  edges <- sort(unique(edges[edges >= a & edges <= b]))
  # QUADPACK can call a piece whose integral is near 0 divergent with
  # abs.tol = 0; such a piece is integrated again with an absolute
  # tolerance of 1e-15 of the integrand's largest value times its length.
  piece <- function(lo, hi) {
    tryCatch(stats::integrate(f, lo, hi, rel.tol = 1e-12, abs.tol = 0,
                              subdivisions = 5000L)$value,
             error = function(e) {
               scale <- max(f(seq(lo, hi, length.out = 401L))) * (hi - lo)
               stats::integrate(f, lo, hi, rel.tol = 1e-12,
                                abs.tol = 1e-15 * scale,
                                subdivisions = 5000L)$value
             })
  }
  total <- 0
  for (i in seq_len(length(edges) - 1L)) {
    if (edges[i + 1L] > edges[i]) total <- total + piece(edges[i], edges[i + 1L])
  }
  total
}

# A quantile, by bisection on .oracle_nct(), each tail matched directly.
.oracle_nct_quantile <- function(p, df, ncp, lo, hi) {
  f <- if (p <= .5) function(q) .oracle_nct(q, df, ncp) - p else
    function(q) (1 - p) - .oracle_nct(q, df, ncp, lower_tail = FALSE)
  stopifnot(f(lo) < 0, f(hi) > 0)
  while (hi - lo > 1e-11 * max(1, abs(lo))) {
    mid <- (lo + hi) / 2
    if (f(mid) > 0) hi <- mid else lo <- mid
  }
  (lo + hi) / 2
}

# A second oracle for the distribution function, the Poisson mixture of
# incomplete beta functions (Guenther, 1978; the series R's pnt() sums up to
# a noncentrality of 37.62, Lenth, 1989), written here in logs over every
# term that matters. For q > 0, with lambda = ncp^2 / 2,
#   p_j = exp(-lambda) lambda^j / j!,
#   r_j = ncp exp(-lambda) lambda^j / (sqrt(2) Gamma(j + 3/2)),
#   P(T > q)  = (1/2) sum_j [p_j I_y(df/2, j + 1/2) + r_j I_y(df/2, j + 1)],
#   P(T <= q) = Phi(-ncp) + (1/2) sum_j [p_j {1 - I_y(df/2, j + 1/2)} +
#                                        r_j {1 - I_y(df/2, j + 1)}],
# with y = df / (q^2 + df) formed directly, so a far tail keeps its
# relative accuracy; for q < 0, P(T <= q; ncp) = P(T > -q; -ncp). With a
# negative effective noncentrality the terms alternate in sign; the
# function then refuses (returns NA) when the rounding error, about 1e-16
# times the sum of the absolute values of the terms, exceeds 1e-9 of the
# value. It is independent of both the integral of R/nc_t_internals.R and
# .oracle_nct() above, and it reaches fewer than one degree of freedom.
.oracle_nct_series <- function(q, df, ncp, lower_tail = TRUE) {
  if (q == 0) return(stats::pnorm(-ncp, lower.tail = lower_tail))
  if (q < 0) return(.oracle_nct_series(-q, df, -ncp, !lower_tail))
  lam <- ncp^2 / 2
  j <- seq(max(0, floor(lam - 40 * sqrt(lam) - 60)), ceiling(lam + 40 * sqrt(lam) + 60))
  y <- df / (q^2 + df)
  pj <- if (lam == 0) as.numeric(j == 0) else exp(-lam + j * log(lam) - lgamma(j + 1))
  rj <- if (ncp == 0) 0 * j else
    sign(ncp) * exp(log(abs(ncp)) - 0.5 * log(2) - lam + j * log(lam) - lgamma(j + 1.5))
  # I_y(df/2, j + 1/2) is 1 - I_x(j + 1/2, df/2): the upper tail's term.
  a <- pj * stats::pbeta(y, df / 2, j + 0.5, lower.tail = !lower_tail)
  b <- rj * stats::pbeta(y, df / 2, j + 1, lower.tail = !lower_tail)
  base <- if (lower_tail) stats::pnorm(-ncp) else 0
  value <- base + 0.5 * sum(a + b)
  scale <- base + 0.5 * sum(abs(a) + abs(b))
  if (value <= 1e7 * .Machine$double.eps * scale) NA_real_ else value
}

# The far tail, P(T > q) for q > 0 and P(T <= q) for q < 0, for a quantile
# as large as the largest double: .oracle_nct() squares (z + ncp) / q, which
# underflows to 0 for |q| beyond about 1e154. Conditioning on Z again,
#   q > 0: P(T > q)  = E_Z[ P(S <= (Z + ncp) / q); Z > -ncp ],
#   q < 0: P(T <= q) = E_Z[ P(S <= -(Z + ncp) / |q|); Z < -ncp ],
# with P(S <= x) = P(V / 2 <= y), y = df x^2 / 2, taken in logs
# (log y = log(df / 2) + 2 log x): pgamma(y, df / 2) when log y >= -600,
# and otherwise its leading term exp(a log y - lgamma(a + 1)), a = df / 2,
# which is the probability to within a factor 1 - a y / (a + 1), that is,
# exactly in double precision there.
.oracle_nct_far <- function(q, df, ncp) {
  a <- df / 2
  s_le <- function(x) {
    ly <- log(df / 2) + 2 * log(x)
    out <- numeric(length(x))
    far <- ly < -600
    out[far] <- exp(a * ly[far] - lgamma(a + 1))
    out[!far] <- stats::pgamma(exp(ly[!far]), a)
    out
  }
  g <- if (q > 0) function(z) stats::dnorm(z) * s_le(exp(log(pmax(z + ncp, 1e-300)) - log(q))) else
    function(z) stats::dnorm(z) * s_le(exp(log(pmax(-(z + ncp), 1e-300)) - log(-q)))
  lo <- if (q > 0) max(-ncp, -40) else -40
  hi <- if (q > 0) 40 else min(-ncp, 40)
  if (lo >= hi) return(0)
  cuts <- unique(c(lo, seq(ceiling(lo), floor(hi)), hi))
  total <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    if (cuts[i + 1L] > cuts[i])
      total <- total + stats::integrate(g, cuts[i], cuts[i + 1L], rel.tol = 1e-12,
                                        abs.tol = 0, subdivisions = 1000L)$value
  }
  total
}
