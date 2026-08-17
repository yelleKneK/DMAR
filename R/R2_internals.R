# =============================================================================
#  Internal helpers for the squared multiple correlation coefficient family
# =============================================================================
#
# Helpers shared by exported R^2 functions:
#
#   ci_R2()                     -- CI for rho^2 (Lee, 1971 for random predictors;
#                                  noncentral F for fixed)
#   ss_power_R2()               -- sample size for the omnibus test of R^2
#   ss_power_R2_sensitivity()   -- Monte Carlo sensitivity analysis for power
#
# Nothing in this file is exported. Names are dot-prefixed by convention to
# mark them as internal; a knowledgeable user can reach them via DMAR:::.name
# for inspection, but they are not part of the package's public API and their
# signatures may change.
#
# Layout:
#   (A) Patnaik (1949) / Lee (1971) two-moment approximation to the CDF of
#       the sample R^2 under random predictors:   .lee_random_R2_cdf
#   (B) Bisection inversion of (A) used for CI construction: .lee_random_R2_bisect
#   (C) Gauss hypergeometric 2F1 used by the exact R^2 moment functions
#       (expected_R2, var_R2, unbiased_R2, ss_aipe_R2, ss_aipe_reg_coef):
#       .hyperg_2F1

# --- (C) Gauss hypergeometric 2F1 (base-R replacement for gsl::hyperg_2F1) ----
# Sums the Gauss hypergeometric series 2F1(a, b; c; x) with the term-ratio
# recurrence  t_{k+1} / t_k = (a + k)(b + k) / ((c + k)(k + 1)) * x.
#
# This replaces gsl::hyperg_2F1() so DMAR carries no GSL system-library
# dependency. DMAR only ever calls it with a = b in {1, 2}, c = (N + const)/2
# (large, positive, > a), and x in [0, 1) (the squared multiple correlation R^2,
# or 1 - R^2), where the series converges. It matches gsl to machine precision
# across that range and is in fact MORE accurate than gsl as x -> 1 (gsl loses
# precision near the boundary; the series was verified term by term against the
# independent Euler integral representation), besides being several times
# faster. max_terms is a generous ceiling: the slowest corner DMAR can reach
# (c near 1.5 with x near 1) converges in about 28,000 terms, and typical calls
# converge in a few dozen. Not exported.
.hyperg_2F1 <- function(a, b, c, x, tol = 1e-15, max_terms = 1e6L) {
  if (x == 0) return(1)
  term  <- 1
  total <- 1
  for (k in 0:(max_terms - 1L)) {
    term  <- term * (a + k) * (b + k) / ((c + k) * (k + 1)) * x
    total <- total + term
    if (abs(term) <= tol * abs(total)) return(total)
  }
  warning("'.hyperg_2F1' did not reach tolerance in ", max_terms,
          " terms; returning the last partial sum.", call. = FALSE)
  total
}

# --- (A) CDF of sample R^2 under random predictors ---------------------------
#
# Lee (1971, JRSS-B, Section 4) gives a two-moment Patnaik (1949) approximation
# to the sampling distribution of R^2 under joint multivariate normality with
# random predictors. Algina and Olejnik (2000) implemented the bisection in
# SAS (ci.smcc.bisec.sas); this is the R port.
#
# Returns P(R^2_sample <= R2_obs | rho^2 = rho2, N, p, random predictors).
#
# Arguments
#   R2_obs : observed sample R^2 (a scalar in (0, 1))
#   rho2   : assumed population rho^2 (a scalar in (0, 1))
#   N      : total sample size
#   p      : number of predictors
#
# The internal variables (PHI_1, PHI_2, PHI_3, GAMMA, g, nu, LAMBDA_U) are
# named to match Lee's (1971) own notation.

.lee_random_R2_cdf <- function(R2_obs, rho2, N, p) {
  R2_tilde <- R2_obs / (1 - R2_obs)
  df_1     <- N - 1
  df_2     <- N - p - 1

  yy    <- rho2 / (1 - rho2)
  GAMMA <- sqrt(1 + yy)
  PHI_1 <- df_1 * (GAMMA^2 - 1) + p
  PHI_2 <- df_1 * (GAMMA^4 - 1) + p
  PHI_3 <- df_1 * (GAMMA^6 - 1) + p
  g     <- (PHI_2 - sqrt(PHI_2^2 - PHI_1 * PHI_3)) / PHI_1
  nu    <- (PHI_2 - 2 * yy * GAMMA * sqrt(df_1 * df_2)) / (g^2)
  LAMBDA_U <- yy * GAMMA * sqrt(df_1 * df_2) / (g^2)
  limit <- df_2 * R2_tilde / (nu * g)
  pf(limit, nu, df_2, ncp = LAMBDA_U)
}

# --- (B) Lee (1971) bisection for CI inversion ------------------------------
#
# Find the population rho^2 value such that
#   P(R^2_sample <= R2_obs | rho^2, N, p) = p_target.
# Used by ci_R2() for both the lower and upper random-predictors limits.
#
# Bisection is over rho^2 on (1e-6, 1 - 1e-6); convergence tolerance and the
# guard against rho^2 reaching the boundary follow the original Algina and
# Olejnik (2000) implementation.

.lee_random_R2_bisect <- function(R2_obs, p_target, N, p,
                                  tol = 1e-5, lower = 1e-6, upper = 1 - 1e-6) {
  x1 <- lower
  x2 <- upper
  x3 <- 0.5
  rhosq <- 0.5

  # Cache the lower-bound CDF value; under bisection it only changes
  # when x1 changes, so recomputing it every iteration (as the
  # previous implementation did) doubles the per-bisection cost.
  diff1 <- .lee_random_R2_cdf(R2_obs, x1, N, p) - p_target
  diff3 <- .lee_random_R2_cdf(R2_obs, x3, N, p) - p_target

  while ((abs(diff3) > tol) && (round(rhosq, 5) < 1)) {
    x3 <- (x1 + x2) / 2
    diff3 <- .lee_random_R2_cdf(R2_obs, x3, N, p) - p_target
    if (diff1 * diff3 < 0) {
      x2 <- x3
    } else {
      x1 <- x3
      diff1 <- diff3
    }
    rhosq <- x3
  }
  rhosq
}


# --- (C) Inner-loop fast path for ci_R2 with random predictors -------------
#
# Returns just the (lower, upper) numeric vector of the R^2 CI under the
# Lee (1971) random-predictors method, skipping the data.frame
# construction and class tagging that ci_R2() does. Used by ss_aipe_R2
# (and other iterative sample size searches) where building thousands
# of data.frames in the inner loop dominates the runtime.
#
# Returns NA boundaries for the degenerate alpha = 0 / 1 cases (rare in
# the iterative search) so callers can detect and dispatch to ci_R2()
# for the full handling.

.ci_R2_random_limits_fast <- function(R2, N, p, alpha_lower, alpha_upper) {
  pul <- alpha_upper
  pll <- 1 - alpha_lower
  if (pul == 0 || pll == 1) return(c(NA_real_, NA_real_))
  # See ci_R2(): the CDF's supremum over rho^2 is at rho^2 = 0, so a target
  # probability above it has no interior solution and the limit is 0.
  C0 <- .lee_random_R2_cdf(R2_obs = R2, rho2 = 0, N = N, p = p)
  ll <- if (C0 <= pll) 0 else .lee_random_R2_bisect(R2_obs = R2, p_target = pll, N = N, p = p)
  ul <- if (C0 <= pul) 0 else .lee_random_R2_bisect(R2_obs = R2, p_target = pul, N = N, p = p)
  if (ll > ul) {
    return(c(NA_real_, NA_real_))
  }
  c(ll, ul)
}


# --- (D) Inner-loop fast path for ci_R2 with fixed predictors --------------
#
# Returns just the (lower, upper) numeric vector of the R^2 CI under
# the fixed-predictors method (inversion of the noncentral F), skipping
# the data.frame construction and class tagging that ci_R2() does.
# Used by ss_aipe_R2's fixed-predictors inner loop.
#
# The bulk of the cost in this path is the conf_limits_ncf root finding
# itself; this helper saves the per-iteration data.frame allocations
# that wrap the inputs and outputs to conf_limits_ncf.

.ci_R2_fixed_limits_fast <- function(R2, N, p, alpha_lower, alpha_upper,
                                     tol = 1e-9) {
  df_1 <- p
  df_2 <- N - p - 1
  F_val <- .convert_R2_f_fast(R2, df_1, df_2)
  Limits <- conf_limits_ncf(F_value = F_val, df_1 = df_1, df_2 = df_2,
                            conf_level = NULL, tol = tol,
                            alpha_lower = alpha_lower,
                            alpha_upper = alpha_upper)
  lower_lambda <- Limits$value[Limits$term == "lower_limit"]
  upper_lambda <- Limits$value[Limits$term == "upper_limit"]
  ll <- .convert_lambda_R2_fast(lower_lambda, N)
  ul <- if (is.infinite(upper_lambda)) 1
        else .convert_lambda_R2_fast(upper_lambda, N)
  c(ll, ul)
}
