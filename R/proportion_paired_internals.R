# Internal helpers shared by the paired proportion family.
#
# Two public functions live in this family: ci_proportion_paired()
# (R/ci_proportion_paired.R), the score interval of Tango (1998) for the
# difference between two proportions measured on the same units, and
# equivalence_proportion_paired() (R/equivalence_proportion_paired.R), the
# score tests of Nam (1997) for noninferiority and equivalence on that
# difference. Both rest on one statistic, so it lives here once.
#
# Notation. Each of N units gives a pair of binary outcomes (x, y). The
# discordant counts are n_10 (x = 1, y = 0) and n_01 (x = 0, y = 1), and
# the estimand is psi = p_10 - p_01 = P(x = 1) - P(y = 1). The concordant
# counts enter only through N, which is why the counts interface takes
# n_10, n_01, and N and nothing else.
#
# Nothing in this file is exported.


# ---- Counts from raw outcomes or from the user -------------------------
#
# Returns c(n_10, n_01, N) from either the two 0/1 vectors or the counts,
# never both. Pairs with a missing value on either member are dropped, as
# the help pages say, because a unit observed under one condition only
# carries no information about the paired difference.
.paired_rd_counts <- function(x, y, n_10, n_01, N) {
  raw    <- !is.null(x) || !is.null(y)
  counts <- !is.null(n_10) || !is.null(n_01) || !is.null(N)
  if (raw && counts)
    .dmar_refuse_argument(c(
      "Supply either 'x' and 'y' or the counts 'n_10', 'n_01', and ",
      "'N', not both."), arg = c("x", "y", "n_10"))
  if (raw) {
    if (is.null(x) || is.null(y))
      .dmar_refuse_argument(c(
        "Supply both 'x' and 'y', or the counts 'n_10', 'n_01', and ",
        "'N'."), arg = c("x", "y", "n_10"))
    if (!(is.numeric(x) || is.logical(x)) || !(is.numeric(y) || is.logical(y)))
      .dmar_refuse_argument(c(
        "'x' and 'y' must be numeric or logical vectors coded 1 and 0 ",
        "(or TRUE and FALSE)."), arg = c("x", "y"))
    if (length(x) != length(y))
      .dmar_refuse_argument(
        "'x' and 'y' must have the same length (one element per unit).",
        arg = c("x", "y"))
    ok <- !(is.na(x) | is.na(y))
    x <- as.numeric(x[ok])
    y <- as.numeric(y[ok])
    if (!all(x %in% c(0, 1)) || !all(y %in% c(0, 1)))
      .dmar_refuse_argument(
        "'x' and 'y' must contain only 0 and 1 (or FALSE and TRUE).",
        arg = c("x", "y"))
    n_10 <- sum(x == 1 & y == 0)
    n_01 <- sum(x == 0 & y == 1)
    N    <- length(x)
    if (N < 1L)
      .dmar_refuse_argument(
        "At least one unit observed under both conditions is required.",
        arg = NULL)
  } else {
    if (is.null(n_10) || is.null(n_01) || is.null(N))
      .dmar_refuse_argument(c(
        "Supply both 'x' and 'y', or the counts 'n_10', 'n_01', and ",
        "'N'."), arg = c("x", "y", "n_10"))
    for (nm in c("n_10", "n_01", "N")) {
      v <- get(nm)
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 ||
          v != round(v))
        .dmar_refuse_argument(c(
          "'", nm, "' must be a single non-negative integer."), arg = NULL)
    }
    if (N < 1 || n_10 + n_01 > N)
      .dmar_refuse_argument(
        "'N' must be positive and at least 'n_10' + 'n_01'.",
        arg = c("N", "n_10", "n_01"))
  }
  c(n_10 = unname(n_10), n_01 = unname(n_01), N = unname(N))
}


# ---- Restricted maximum likelihood estimate of p_01 --------------------
#
# Under psi = psi0 the multinomial log likelihood, in the nuisance
# parameter q = p_01, is
#   n_10 log(q + psi0) + n_01 log(q) + (N - n_10 - n_01) log(1 - 2q - psi0),
# and its maximizer is the larger root of A q^2 + B q + C = 0 with
#   A = 2N, B = (2N - n_10 + n_01) psi0 - (n_10 + n_01),
#   C = -n_01 psi0 (1 - psi0)
# (Nam, 1997, equation 1). The root is taken in the form that avoids
# cancellation: when B > 0 the textbook (-B + sqrt(D)) / (2A) subtracts
# two nearly equal numbers whenever 4AC is small beside B^2, so the
# equivalent -2C / (B + sqrt(D)), from the product of the roots, is used
# there. Vectorized over psi0.
.paired_rd_p01 <- function(n_10, n_01, N, psi0) {
  A <- 2 * N
  B <- (2 * N - n_10 + n_01) * psi0 - (n_10 + n_01)
  C <- -n_01 * psi0 * (1 - psi0)
  s <- sqrt(pmax(B^2 - 4 * A * C, 0))
  ifelse(B > 0, -2 * C / (B + s), (-B + s) / (2 * A))
}


# ---- The score statistic -----------------------------------------------
#
#   Z(psi0) = (n_10 - n_01 - N psi0) / sqrt(N {2 p01 + psi0 (1 - psi0)}),
#
# with p01 the restricted estimate above: Nam's (1997) equation 3, which
# his Appendix 2 shows equals the efficient score statistic, and the
# statistic whose inversion is Tango's (1998) interval. The denominator is
# N times the variance of the per-unit difference, p_10 + p_01 - psi^2,
# evaluated at the restricted estimates, since p_10 = p_01 + psi. At
# psi0 = 0 it reduces to McNemar's statistic, (n_10 - n_01) / sqrt(n_10 +
# n_01). Z decreases in psi0. Vectorized over psi0.
.paired_rd_score <- function(n_10, n_01, N, psi0) {
  p01 <- .paired_rd_p01(n_10, n_01, N, psi0)
  num <- n_10 - n_01 - N * psi0
  # At psi0 equal to the estimate the numerator is zero; with no
  # discordant pair the denominator is zero there as well, and the
  # statistic is 0 (the data sit exactly at the hypothesized value).
  ifelse(num == 0, 0, num / sqrt(N * (2 * p01 + psi0 * (1 - psi0))))
}


# ---- Score limits ------------------------------------------------------
#
# The lower limit is the psi0 at which Z(psi0) = +q and the upper limit the
# psi0 at which Z(psi0) = -q, each found by uniroot() between the estimate,
# where Z is 0, and the edge of the parameter space. When every pair is
# (0, 1) the estimate is -1 and so is the lower limit; when every pair is
# (1, 0) the estimate and the upper limit are 1. Away from those cases Z
# grows without bound toward the edges, so each bracket holds a root; the
# fallback to -1 or 1 only guards against a statistic that floating point
# arithmetic could not separate from q at the edge.
.paired_rd_score_limits <- function(n_10, n_01, N, q) {
  est <- (n_10 - n_01) / N
  eps <- 1e-12
  f <- function(p, target) .paired_rd_score(n_10, n_01, N, p) - target
  lower <- if (n_01 == N || f(-1 + eps, q) <= 0) {
    -1
  } else {
    stats::uniroot(f, c(-1 + eps, est), target = q, tol = 1e-12)$root
  }
  upper <- if (n_10 == N || f(1 - eps, -q) >= 0) {
    1
  } else {
    stats::uniroot(f, c(est, 1 - eps), target = -q, tol = 1e-12)$root
  }
  c(lower, upper)
}
