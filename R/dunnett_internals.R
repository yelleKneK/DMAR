# Shared internals for Dunnett's many-to-one comparisons: the exact CDF of the
# maximum (or maximum modulus) of the balanced-design statistics, and the
# critical-value inversion built on it. Both cv_dunnett() (critical value) and
# ci_dunnett() (adjusted p-values) use these, so they live in one internals file
# rather than being duplicated across the two public files.
#
# The statistics T_1, ..., T_m are the treatment-versus-control comparisons of a
# balanced one-way design; they are a multivariate t with common correlation
# rho = 1/2 (each shares the one control mean) on df error degrees of freedom.
# A single common correlation admits the one-factor representation
#   Z_i = sqrt(rho) W + sqrt(1 - rho) U_i,   T_i = Z_i / S,
# with a shared standard-normal factor W, independent standard-normal U_i, and a
# shared scale S = sqrt(chi^2_df / df). Conditioning on W and S makes the U_i
# independent, so the m-dimensional probability factorizes and collapses to two
# nested one-dimensional integrals, evaluated deterministically here (no Monte
# Carlo). With rho = 1/2 both sqrt(rho) and sqrt(1 - rho) equal 1/sqrt(2).

# P(max_i T_i <= q)  (two_sided = FALSE)  or  P(max_i |T_i| <= q)  (two_sided
# = TRUE), for m equicorrelated (rho = 1/2) balanced-design Dunnett statistics
# on df error degrees of freedom.  df = Inf gives the known-variance limit.
.dunnett_cdf <- function(q, df, m, two_sided) {
  a <- sqrt(0.5)                                   # sqrt(rho) = sqrt(1 - rho)
  # P(all comparisons on the correct side | shared factor w, scale s).
  inner <- function(s) stats::integrate(function(w) {
      p <- if (two_sided)
             stats::pnorm((q * s - a * w) / a) - stats::pnorm((-q * s - a * w) / a)
           else
             stats::pnorm((q * s - a * w) / a)
      p[p < 0] <- 0
      p^m * stats::dnorm(w)
    }, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value
  if (is.infinite(df)) return(inner(1))
  # Integrate the chi squared error variate between its extreme quantiles rather
  # than over (0, Inf). At a large df its density is a narrow bump far from the
  # origin, and an adaptive rule asked to cover the whole half-line samples its
  # way past the mass and returns zero (at df = 200 it does exactly that), which
  # leaves the root finder in .dunnett_upper_cv() with no bracket.
  stats::integrate(function(u)
      vapply(sqrt(u / df), function(si) inner(si), numeric(1)) * stats::dchisq(u, df),
      lower = stats::qchisq(1e-12, df),
      upper = stats::qchisq(1e-12, df, lower.tail = FALSE),
      rel.tol = 1e-8)$value
}

# Upper-alpha Dunnett critical value: the quantile d with .dunnett_cdf(d) = 1 -
# alpha. Deterministic; m = 1 is the ordinary one- or two-sided t (or normal).
.dunnett_upper_cv <- function(alpha, df, m, two_sided) {
  if (m == 1)
    return(if (two_sided) (if (is.infinite(df)) stats::qnorm(1 - alpha / 2)
                           else                 stats::qt(1 - alpha / 2, df))
           else           (if (is.infinite(df)) stats::qnorm(1 - alpha)
                           else                 stats::qt(1 - alpha, df)))
  target <- 1 - alpha
  # Bracket: the value lies between the single-test t (lower) and the Bonferroni
  # t (upper); pad slightly so the m = 1 equalities do not sit on a bracket end.
  if (is.infinite(df)) {
    lo <- if (two_sided) stats::qnorm(1 - alpha / 2)       else stats::qnorm(1 - alpha)
    hi <- if (two_sided) stats::qnorm(1 - alpha / (2 * m)) else stats::qnorm(1 - alpha / m)
  } else {
    lo <- if (two_sided) stats::qt(1 - alpha / 2, df)       else stats::qt(1 - alpha, df)
    hi <- if (two_sided) stats::qt(1 - alpha / (2 * m), df) else stats::qt(1 - alpha / m, df)
  }
  stats::uniroot(function(d) .dunnett_cdf(d, df, m, two_sided) - target,
                 c(lo * 0.99, hi * 1.01), tol = 1e-9)$root
}
