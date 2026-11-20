# =============================================================================
#  Internal helpers for the squared multiple correlation coefficient family
# =============================================================================
#
# Helpers shared by exported R^2 functions:
#
#   ci_R2()                     -- CI for rho^2 (Lee, 1971 for random predictors;
#                                  noncentral F for fixed)
#   ss_power_R2()               -- sample size for the omnibus test of R^2
#   ss_power_R2_change()        -- sample size for the test of a change in R^2
#   ss_power_R2_sensitivity()   -- Monte Carlo sensitivity analysis for power
#
# Nothing in this file is exported. Names are dot-prefixed by convention to
# mark them as internal; a knowledgeable user can reach them via DMAR:::.name
# for inspection, but they are not part of the package's public API and their
# signatures may change.
#
# Layout:
#   (A) Lee (1971, Section 5.1, Equation 25): three-moment noncentral F
#       approximation to the distribution of the sample R^2 under random
#       predictors:   .lee_random_R2_parms, .lee_random_R2_cdf,
#       .lee_random_R2_quantile
#   (B) Bisection inversion of (A) used for CI construction: .lee_random_R2_bisect
#   (C) Gauss hypergeometric 2F1 used by the exact R^2 moment functions
#       (expected_R2, var_R2, unbiased_R2, ss_aipe_R2, ss_aipe_reg_coef):
#       .hyperg_2F1
#   (D), (E) Inner-loop fast paths for the limits of ci_R2()
#   (F) Power of the model comparison F test of a change in R^2 and the
#       sample size search, shared by ss_power_R2() and ss_power_R2_change():
#       .R2_change_power, .R2_change_power_guarded, .R2_imprecise_message,
#       .R2_power_search_N
#   (G) The half-width search of ss_aipe_R2(): .R2_half_width_search_N
#   (H) ci_R2() as the planner ss_aipe_R2() uses it: .ci_R2_planning

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
  .dmar_warn(c(
    "'.hyperg_2F1' did not reach tolerance in ", max_terms,
    " terms; returning the last partial sum."), "result_qualified")
  total
}

# --- (A) CDF of sample R^2 under random predictors ---------------------------
#
# Lee (1971, Section 5.1, Equation 25): three-moment noncentral F approximation
# to the sampling distribution of R^2 under joint multivariate normality with
# random predictors. Algina and Olejnik (2000) implemented the bisection in
# SAS (ci.smcc.bisec.sas); this is the R port.
#
# With n = N - 1, gamma^2 = 1 / (1 - rho^2), and
#   theta_j = n (gamma^(2j) - 1) + p,   j = 1, 2, 3,
# the approximation is
#   R^2 / (1 - R^2)  ~  {nu g / (n - p)} F'(nu, n - p; lambda),
# where F' is a noncentral F variate and g, nu and lambda are functions of
# the theta_j. In the code PHI_j is theta_j, GAMMA is gamma, and LAMBDA_U is
# lambda; g and nu keep Lee's names.

# .lee_random_R2_parms(rho2, N, p)
#
# The parameters of Lee's approximation at a population rho^2.
#
# The square root of PHI_2^2 - PHI_1 PHI_3, with PHI_3 = df_1 * (GAMMA^6 - 1) + p,
# is not computed from that difference. Writing
# x = gamma^2 - 1 = rho^2 / (1 - rho^2) and n = N - 1,
#   PHI_2^2 - PHI_1 PHI_3 = n (n - p) x^2 (1 + x)
#                         = df_1 * df_2 * yy^2 * GAMMA^2,
# so the root is yy * GAMMA * sqrt(df_1 * df_2) exactly. The subtraction
# cancels for a small rho^2 and came out negative near rho^2 = 1e-8, where
# the square root was NaN; the closed form has no subtraction.
#
# Arguments
#   rho2 : population rho^2 (a scalar in [0, 1))
#   N    : total sample size
#   p    : number of predictors
# Returns
#   A list with g, nu, LAMBDA_U and df_2 = N - p - 1. At rho2 = 0 these are
#   g = 1, nu = p and LAMBDA_U = 0, the central F of the usual omnibus test.

.lee_random_R2_parms <- function(rho2, N, p) {
  df_1 <- N - 1
  df_2 <- N - p - 1

  yy    <- rho2 / (1 - rho2)
  GAMMA <- sqrt(1 + yy)
  PHI_1 <- df_1 * (GAMMA^2 - 1) + p
  PHI_2 <- df_1 * (GAMMA^4 - 1) + p
  root  <- yy * GAMMA * sqrt(df_1 * df_2)
  g     <- (PHI_2 - root) / PHI_1
  nu    <- (PHI_2 - 2 * root) / (g^2)
  LAMBDA_U <- root / (g^2)
  list(g = g, nu = nu, LAMBDA_U = LAMBDA_U, df_2 = df_2)
}

# .lee_random_R2_cdf(R2_obs, rho2, N, p)
#
# Returns P(R^2_sample <= R2_obs | rho^2 = rho2, N, p, random predictors).
#
# Arguments
#   R2_obs : observed sample R^2 (a scalar in (0, 1))
#   rho2   : assumed population rho^2 (a scalar in [0, 1))
#   N      : total sample size
#   p      : number of predictors

.lee_random_R2_cdf <- function(R2_obs, rho2, N, p) {
  parms    <- .lee_random_R2_parms(rho2 = rho2, N = N, p = p)
  R2_tilde <- R2_obs / (1 - R2_obs)
  limit    <- parms$df_2 * R2_tilde / (parms$nu * parms$g)
  pf(limit, parms$nu, parms$df_2, ncp = parms$LAMBDA_U)
}

# .lee_random_R2_quantile(prob, rho2, N, p)
#
# The inverse of .lee_random_R2_cdf() in its first argument: the value of the
# sample R^2 below which a proportion prob of its sampling distribution lies
# at rho^2 = rho2. The CDF is pf() of an increasing function of R^2, so the
# quantile is closed form in qf(). At rho2 = 0 the central F quantile is
# used, which is the critical value of the usual omnibus test.
# .R2_change_power() takes its critical value from here at rho2 = 0, for
# ss_power_R2() and ss_power_R2_change(). The branch for a nonzero rho2 has
# no exported caller; it is the dual of the one-sided limit of ci_R2().
#
# Arguments
#   prob : cumulative probability (a scalar in (0, 1))
#   rho2 : population rho^2 (a scalar in [0, 1))
#   N    : total sample size
#   p    : number of predictors

.lee_random_R2_quantile <- function(prob, rho2, N, p) {
  parms <- .lee_random_R2_parms(rho2 = rho2, N = N, p = p)
  F_q <- if (rho2 == 0) {
    qf(prob, p, parms$df_2)
  } else {
    qf(prob, parms$nu, parms$df_2, ncp = parms$LAMBDA_U)
  }
  R2_tilde <- parms$nu * parms$g * F_q / parms$df_2
  R2_tilde / (1 + R2_tilde)
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


# --- (D) Inner-loop fast path for ci_R2 with random predictors -------------
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
  # probability above it has no interior solution. The lower limit is then
  # 0, and the upper limit does not exist and is NA, as in ci_R2(). (At the
  # expected sample R^2 that ss_aipe_R2() passes, the CDF at rho^2 = 0 is
  # above one half, so with tails of at most one half neither case arises
  # there.)
  C0 <- .lee_random_R2_cdf(R2_obs = R2, rho2 = 0, N = N, p = p)
  ll <- if (C0 <= pll) 0 else .lee_random_R2_bisect(R2_obs = R2, p_target = pll, N = N, p = p)
  ul <- if (C0 <= pul) NA_real_ else .lee_random_R2_bisect(R2_obs = R2, p_target = pul, N = N, p = p)
  if (!is.na(ul) && ll > ul) {
    return(c(NA_real_, NA_real_))
  }
  c(ll, ul)
}


# --- (E) Inner-loop fast path for ci_R2 with fixed predictors --------------
#
# Returns just the (lower, upper) numeric vector of the R^2 CI under
# the fixed-predictors method (inversion of the noncentral F), skipping
# the data.frame construction and class tagging that ci_R2() does.
# Used by ss_aipe_R2's fixed-predictors inner loop.
#
# The bulk of the cost in this path is the ci_nc_F root finding
# itself; this helper saves the per-iteration data.frame allocations
# that wrap the inputs and outputs to ci_nc_F.

.ci_R2_fixed_limits_fast <- function(R2, N, p, alpha_lower, alpha_upper,
                                     tol = 1e-9) {
  df_1 <- p
  df_2 <- N - p - 1
  F_val <- .convert_R2_f_fast(R2, df_1, df_2)
  Limits <- ci_nc_F(F_value = F_val, df_1 = df_1, df_2 = df_2,
                            conf_level = NULL, tol = tol,
                            alpha_lower = alpha_lower,
                            alpha_upper = alpha_upper)
  lower_lambda <- Limits$value[Limits$term == "ci_lower"]
  upper_lambda <- Limits$value[Limits$term == "ci_upper"]
  ll <- .convert_lambda_R2_fast(lower_lambda, N)
  ul <- if (is.infinite(upper_lambda)) 1
        else .convert_lambda_R2_fast(upper_lambda, N)
  c(ll, ul)
}


# --- (F) Power and sample size for the test of a change in R^2 ---------------
#
# The model comparison F test of u = p_full - p_restricted predictors added
# to a model that has p_restricted predictors. With p_restricted = 0 it is
# the omnibus test of all p_full predictors.

# .R2_change_power(N, f2, p_full, p_restricted, alpha_level, random_predictors)
#
# Power of the test at total sample size N.
#
# Fixed predictors: F is noncentral F with u and N - p_full - 1 degrees of
# freedom and noncentrality N f2. The numerator sum of squares is sigma^2
# times a noncentral chi square on u degrees of freedom with noncentrality
# mu'(P_full - P_restricted) mu / sigma^2, P the projection onto a model's
# columns and mu the mean vector, which is N f2 when the two R^2 values are
# those of the design itself (sums of squares over N). For one added
# predictor this is Kelley and Maxwell (2008, Equations 19, 21, and 31).
# (Equation 22 of the chapter prints the denominator of f2 as the full
# model's R^2; Equations 19, 31, and 32 and the worked example of pp. 181 to
# 182 use one minus it, as the code does.)
#
# Random predictors (multivariate normal): the sample squared multiple
# partial correlation of the outcome with the u added predictors given the
# other p_restricted has the distribution of a sample squared multiple
# correlation with u predictors and a sample size of N - p_restricted
# (Gatsonis and Sampson, 1989, p. 518). That distribution is evaluated with
# Lee's approximation, (A) above, at rho^2 = f2 / (1 + f2), the population
# squared multiple partial correlation; the critical value is the central F
# quantile, which the quantile function returns at rho2 = 0.
#
# Arguments
#   N                 : total sample size
#   f2                : Cohen's f^2 for the change, (rho^2_full -
#                       rho^2_restricted) / (1 - rho^2_full)
#   p_full            : predictors in the full model
#   p_restricted      : predictors in the restricted model (0 for none)
#   alpha_level       : Type I error rate
#   random_predictors : TRUE or FALSE
# Returns
#   The power; 0 when N leaves no error degrees of freedom.

.R2_change_power <- function(N, f2, p_full, p_restricted, alpha_level,
                             random_predictors) {
  u <- p_full - p_restricted
  df_2 <- N - p_full - 1
  if (df_2 < 1) return(0)
  if (random_predictors) {
    N_star  <- N - p_restricted
    R2_crit <- .lee_random_R2_quantile(1 - alpha_level, rho2 = 0,
                                       N = N_star, p = u)
    1 - .lee_random_R2_cdf(R2_obs = R2_crit, rho2 = f2 / (1 + f2),
                           N = N_star, p = u)
  } else {
    f_crit <- qf(1 - alpha_level, df1 = u, df2 = df_2)
    1 - pf(f_crit, df1 = u, df2 = df_2, ncp = N * f2)
  }
}

# .R2_change_power_guarded(...)
#
# .R2_change_power() with R's two noncentral beta warnings turned into a
# missing value. The noncentral beta routine, on which pf() rests, loses
# accuracy once the noncentrality passes about one million and says so with
# one of two warnings ("full precision may not have been achieved in
# 'pnbeta'", "convergence failed in 'pnbeta'"). A power computed through
# either can be wrong by an amount that is not small, so the warnings are
# muffled here and NA is returned in place of a number; the callers refuse a
# missing power. Other warnings pass through.

.R2_change_power_guarded <- function(N, f2, p_full, p_restricted, alpha_level,
                                     random_predictors) {
  state <- new.env(parent = emptyenv())
  state$imprecise <- FALSE
  power_now <- withCallingHandlers(
    .R2_change_power(N, f2, p_full, p_restricted, alpha_level,
                     random_predictors),
    warning = function(w) {
      if (grepl("pnbeta", conditionMessage(w), fixed = TRUE)) {
        state$imprecise <- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  if (state$imprecise) NA_real_ else power_now
}

# .R2_imprecise_message(N, f2)
#
# The reason given when .R2_change_power_guarded() returned NA at N.

.R2_imprecise_message <- function(N, f2) {
  paste0("At a sample size of ",
         format(N, scientific = FALSE, big.mark = ","),
         " the noncentrality parameter of the test is about ",
         format(signif(N * f2, 3), scientific = FALSE, big.mark = ","),
         ", beyond the point (about one million) to which R computes ",
         "the noncentral F distribution accurately, so no power is ",
         "returned. This happens only at very large sample sizes.")
}

# .R2_power_search_N(power_at_N, N_min, desired_power, imprecise_message,
#                    N_max = 1e7)
#
# The smallest whole N >= N_min at which power_at_N(N) is no less than
# desired_power.
#
# The smallest admissible size is evaluated first, so a target the minimum
# already attains returns that minimum. Otherwise the size is doubled until
# the target is met and the smallest such size is then found by bisection on
# whole numbers; power is nondecreasing in N for these tests.
#
# A size at which the power cannot be computed accurately (NA from
# power_at_N) is treated as lying above the sizes at which it can, since the
# noncentrality grows with the size. Such a size closes the bracket from
# above without counting as having met the target; the bisection then either
# finds a size below it that meets the target, or ends with the bracket still
# closed by a size without a power, and the search stops with the reason.
#
# Arguments
#   power_at_N        : function of N returning the power, or NA
#   N_min             : smallest admissible total sample size
#   desired_power     : the target
#   imprecise_message : function of N returning the reason for an NA
#   N_max             : the search does not go beyond this size
# Returns
#   list(N = , power = ).

.R2_power_search_N <- function(power_at_N, N_min, desired_power,
                               imprecise_message, N_max = 1e7) {
  N_low <- N_min
  power_low <- power_at_N(N_low)
  if (is.na(power_low)) {
    .dmar_refuse(imprecise_message(N_low), "computation_failed")
  }
  if (power_low >= desired_power) {
    return(list(N = N_low, power = power_low))
  }
  repeat {
    if (N_low >= N_max) {
      .dmar_refuse(
        "Failed to reach 'desired_power' within a reasonable sample size.",
        "computation_failed")
    }
    N_high <- min(2 * N_low, N_max)
    power_high <- power_at_N(N_high)
    if (is.na(power_high) || power_high >= desired_power) break
    N_low <- N_high
  }
  # Here power at N_low is below the target and power at N_high meets it or
  # could not be computed.
  while (N_high - N_low > 1) {
    N_mid <- floor((N_low + N_high) / 2)
    power_mid <- power_at_N(N_mid)
    if (is.na(power_mid) || power_mid >= desired_power) {
      N_high <- N_mid
      power_high <- power_mid
    } else {
      N_low <- N_mid
    }
  }
  if (is.na(power_high)) {
    .dmar_refuse_argument(c(
      "The sample size that attains 'desired_power' is above ",
      format(N_low, scientific = FALSE, big.mark = ","), ". ",
      imprecise_message(N_high)), arg = "desired_power")
  }
  list(N = N_high, power = power_high)
}


# --- (G) The half-width search of ss_aipe_R2() -------------------------------
#
# .R2_half_width_search_N(distances_at_N, N_min, width, tol, N_max = 1e7)
#
# The smallest whole N >= N_min such that the distance from the expected
# sample R^2 to one confidence limit is no larger than width at N and at
# every larger sample size: one more than the largest N at which the
# distance is larger than width, and N_min when there is no such N. A
# distance is within width when it exceeds it by no more than tol, and a
# distance that cannot be computed (NA) is not within it.
#
# The first N of an upward search at which the distance is within width is
# not that sample size in general. At a small N with several predictors the
# expected sample R^2 is far above the population value and the upper limit
# is close to that expected value or below it, so the distance to the upper
# limit is small or negative there; it rises to one maximum as N grows and
# falls from then on.
# The distance to the lower limit falls from the smallest N on (at 99
# percent with one predictor and a population R^2 of .70 to .90 it first
# rises, by less than .005, over the first few N).
#
# The search rests on three properties of these distances as functions of
# N. They were checked numerically and are not proved: on 1,440 profiles
# (random and fixed predictors; 12 values of the population R^2 from .02 to
# .95; 1 to 20 predictors; 90, 95, and 99 percent), at every N from 2p + 3
# to 1,200 and on a geometric grid from there to 60,000.
#
#   (a) The covering width, the largest of the full width of the interval
#       and the two distances, does not rise with N. Neither distance is
#       larger than it, so at and beyond a sample size at which it is
#       within width both distances are within width. (The full width did
#       not rise in any profile either, and the lower limit was never above
#       the expected sample R^2, so the full width bounds the distance to
#       the upper limit. It does not bound the distance to the lower limit:
#       in 437 of the profiles, all with 7 or more predictors, the upper
#       limit is below the expected sample R^2 at some N, and there that
#       distance is larger than the full width.)
#   (b) The distance rises to at most one maximum, so the sample sizes at
#       which it is larger than width are consecutive.
#   (c) The limits are computed to a tolerance, so (b) holds only up to
#       their rounding. With random predictors a bisection stops when the
#       distribution function is within 1e-5 of its target, and where the
#       distance passes width on its way down it can pass it more than
#       once. In windows of consecutive N beginning at 300 to 30,000 the
#       sample sizes at which that happened were at most 3 apart for N up
#       to 3,000 and N / 800 beyond at 90 to 99 percent, and at most 8
#       apart for N up to 1,000 and N / 130 beyond at 99.9 percent. With
#       fixed predictors it did not happen. Rounding also limits step (2)
#       below: when width is within about 1e-5 of the maximum of the
#       distance, a single sample size beside the maximum can miss width
#       by rounding while every sample size the search evaluates meets it
#       (1 of 6,560 plans whose width was set 1e-5 to 5e-3 below the
#       maximum), and the search then returns N_min.
#
# The steps: (1) double N from N_min until the covering width is within
# width, which bounds the answer from above by (a); (2) find a sample size
# below that bound at which the distance is larger than width, N_min itself
# when it is one, and otherwise by a ternary search for the maximum of the
# distance that stops at the first such size it meets, and return N_min
# when there is none; (3) bisect between that sample size and the bound,
# which by (b) ends one above the largest sample size that misses, up to
# the rounding of (c); (4) step upward from there until the distance has
# been within width at max(10, N / 50) consecutive sample sizes, and return
# the first of them. The search evaluates the limits at a few dozen sample
# sizes, and at about a hundred for a plan of several thousand.
#
# The search stops with an error for width when step (1) reaches N_max
# without a bound, or when step (4) reaches N_max. Neither means that no
# sample size meets width, and each message says what the search did not
# find. The bound of step (1) is about four times the sample size the
# distance needs (the full width is about twice the distance, and both
# fall as 1 / sqrt(N)), and up to twice that by the doubling; and at a
# very large N the limits are beyond what R computes accurately. R's
# documentation of pf() (?FDist) says that its code for a nonzero ncp "is
# principally intended to be used for moderate values of ncp: it will not
# be highly accurate, especially in the tails, for large values", and the
# noncentral pbeta() it calls is Lenth's (1987) algorithm with the remarks
# of Frick (1990) and Lam (1995). With 5 predictors, 95 percent, and a
# population R^2 of .30 to .98, the lower limit for fixed predictors goes
# wrong without a warning once N * rho^2 / (1 - rho^2) passes about two
# million (at .98 the full width is .00068 at N = 26,624, a noncentrality
# of 1.3 million, and .00117 at N = 53,248, 2.6 million, where it should
# be smaller), and for random predictors pf() warns that the noncentral
# beta distribution did not converge from an N of one and a half to four
# million. So the message of step (1) says that the width is narrower
# than can be planned with R's distribution functions. A call whose
# plan is far below N_max can therefore be refused: with fixed predictors,
# a population R^2 of .98, 5 predictors, and a width of .0005 the distance
# to the upper limit is within width from N = 11,351, and step (1) finds
# no bound, because the full width is .00068 at N = 26,624 and the lower
# limit is wrong at N = 53,248. Of 600 calls with a width of .002 to .005
# (random and fixed predictors; a population R^2 of .10 to .90; 1 to 20
# predictors; 90, 95, and 99 percent) 21 were refused: 17 of the 120 with
# a width of .002, 4 of the 120 with .0025, and none with .003 or more.
#
# Arguments
#   distances_at_N : function of N returning c(the distance planned for,
#                    the covering width)
#   N_min          : smallest sample size the search considers
#   width          : the target for the distance
#   tol            : a distance above width by no more than this is within it
#   N_max          : the search does not go beyond this size, and stops
#                    with an error for width when it would have to
# Returns
#   The sample size.

.R2_half_width_search_N <- function(distances_at_N, N_min, width, tol,
                                    N_max = 1e7) {
  within <- function(x) !is.na(x) && x - width <= tol
  N_max_text <- format(N_max, scientific = FALSE, big.mark = ",")
  refuse_width <- function(what) {
    .dmar_refuse_argument(c(
      "'width' is ", format(width), ", and the search ", what,
      " Give a larger 'width'."),
      arg = "width", value = width)
  }
  at_min <- distances_at_N(N_min)

  # (1) A sample size from which both distances stay within width.
  N_high <- N_min
  covering <- at_min[2L]
  while (!within(covering)) {
    if (N_high >= N_max) .dmar_refuse_argument(c(
      "'width' is ", format(width), ", narrower than can be planned with ",
      "R's distribution functions. To confirm that the distance from the ",
      "expected sample R2 to the limit stays within 'width', the search ",
      "needs the confidence limits at sample sizes several times the plan: ",
      "either where the noncentrality parameter is in the millions, beyond ",
      "the range in which R's noncentral F and beta distribution functions ",
      "(pf() and pbeta()) are accurate (see ?pf), or above ", N_max_text,
      ", the largest sample size the search tries. This does not mean that ",
      "no sample size meets 'width'. Give a larger 'width'."),
      arg = "width", value = width)
    N_high <- min(2 * N_high, N_max)
    covering <- distances_at_N(N_high)[2L]
  }
  if (N_high == N_min) return(N_min)

  # (2) A sample size below N_high at which the distance is larger than
  # width.
  if (!within(at_min[1L])) {
    N_low <- N_min
  } else {
    N_low <- NA_real_
    left <- N_min
    right <- N_high
    while (is.na(N_low) && right - left > 2) {
      third <- (right - left) %/% 3
      N_a <- left + third
      N_b <- right - third
      d_a <- distances_at_N(N_a)[1L]
      d_b <- distances_at_N(N_b)[1L]
      if (!within(d_a)) {
        N_low <- N_a
      } else if (!within(d_b)) {
        N_low <- N_b
      } else if (d_a < d_b) {
        left <- N_a
      } else {
        right <- N_b
      }
    }
    # The bracket has closed on the maximum without meeting a distance
    # larger than width (its one interior point was N_a or N_b of the last
    # pass), so the distance is within width from N_min on.
    if (is.na(N_low)) return(N_min)
  }

  # (3) The distance is larger than width at N_low and within it at N_high.
  while (N_high - N_low > 1) {
    N_mid <- (N_low + N_high) %/% 2
    if (within(distances_at_N(N_mid)[1L])) N_high <- N_mid else N_low <- N_mid
  }

  # (4) N_high is within width; so must be the sample sizes after it.
  run_needed <- max(10, ceiling(N_high / 50))
  run <- 1
  N <- N_high
  while (run < run_needed) {
    if (N >= N_max) refuse_width(c(
      "reached ", N_max_text, ", the largest sample size it tries, before ",
      "it had confirmed that the distance from the expected sample R2 to ",
      "the limit stays within 'width'."))
    N <- N + 1
    if (within(distances_at_N(N)[1L])) {
      run <- run + 1
    } else {
      N_high <- N + 1
      run <- 0
    }
  }
  N_high
}


# --- (H) ci_R2() as the planner ss_aipe_R2() uses it -------------------------
#
# .ci_R2_planning(...) takes the arguments of ci_R2() and returns what it
# returns, except in one case. The assurance step of Kelley (2008), as
# MBESS and DMAR 1.0.0 carry it out, evaluates the width of the interval at
# a confidence limit of R^2 used as a planning value, and with a small
# population R^2 that limit can be 0 or close to it. With random predictors
# ci_R2() reports the upper limit there as undefined (NA, with a warning,
# when the planning value is at or below the alpha_upper quantile of its
# distribution at rho^2 = 0). The procedure takes the interval there as
# running from 0 to 0, of width 0, and so does this function, without the
# warning, so that the plans are those of the procedure. Every other result
# of ci_R2(), and every other warning, passes through unchanged.
.ci_R2_planning <- function(...) {
  # The handler records the case in an environment made for it, so that it
  # updates it by ordinary assignment.
  state <- new.env(parent = emptyenv())
  state$undefined <- FALSE
  r <- withCallingHandlers(
    ci_R2(...),
    warning = function(w) {
      if (grepl("population squared multiple correlation of 0 the probability",
                conditionMessage(w), fixed = TRUE)) {
        state$undefined <- TRUE
        invokeRestart("muffleWarning")
      }
    })
  if (state$undefined) {
    up <- r$term == "ci_upper"
    a <- list(...)
    alpha_upper <- if (!is.null(a$alpha_upper)) a$alpha_upper else
      (1 - (if (is.null(a$conf_level)) .95 else a$conf_level)) / 2
    r$value[up] <- 0
    r$prob_less[up] <- 1 - alpha_upper
    r$prob_greater[up] <- alpha_upper
  }
  r
}
