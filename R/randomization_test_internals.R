# Internal helpers shared by the randomization test family.
#
# Two public functions live in this family: randomization_test() for two
# independent groups (R/randomization_test.R) and
# randomization_test_paired() for paired observations
# (R/randomization_test_paired.R). The counting rule that turns a
# reference distribution into a p-value is identical for both, so it lives
# here rather than being written twice. The vectorized permutation
# machinery below is what makes the independent-groups test fast enough to
# invert for a confidence interval, which requires evaluating the whole
# reference distribution once per candidate shift.
#
# Nothing in this file is exported.


# ---- Thresholds for exact enumeration ----------------------------------
#
# The independent-groups test enumerates choose(N, n_1) reassignments. Two
# constants govern the choice between enumeration and Monte Carlo, and both
# are documented on ?randomization_test so a user can predict which branch
# a given design takes.
#
#   .rt_exact_threshold  the largest choose(N, n_1) for which the default
#                        exact = NULL enumerates rather than samples. At
#                        50,000 a balanced design is enumerated through 9
#                        observations per group (choose(18, 9) = 48,620),
#                        and the whole reference distribution still fits
#                        comfortably in memory for the confidence interval
#                        inversion.
#   .rt_exact_cap        the largest choose(N, n_1) accepted when exact =
#                        TRUE is forced. Beyond a million reassignments the
#                        index matrix alone runs to tens of megabytes and
#                        the inversion becomes impractical, so the request
#                        is refused rather than silently thrashing.

.rt_exact_threshold <- 50000
.rt_exact_cap       <- 1000000


# ---- Counting rule (shared by both public functions) -------------------
#
# Counts how many members of the reference distribution T_null are at least
# as extreme as the observed statistic T_obs, in the direction named by
# 'alternative'. The tolerance exists because a reassignment that is a
# genuine tie with the observed one can miss equality by a few units in the
# last binary place. It is scaled by the size of the reference distribution
# itself, never floored at an absolute constant, because a randomization
# test must be exactly equivariant under rescaling the response: multiply
# every score by c > 0 and the statistic and the whole null distribution
# scale by c, so the set of reassignments counted as at least as extreme,
# and hence the p-value, must not change. A floor of 1 would make the
# tolerance absolute, so a mean difference measured in millionths would
# count the entire distribution as tied while the same data in original
# units would not. Real neighboring values of a permutation distribution
# are separated by far more than this.

.rt_extreme <- function(T_null, T_obs, alternative,
                        tol = sqrt(.Machine$double.eps) *
                          max(abs(T_obs), max(abs(T_null)), 0)) {
  switch(alternative,
    two_sided = sum(abs(T_null) >= abs(T_obs) - tol),
    less      = sum(T_null <=  T_obs + tol),
    greater   = sum(T_null >=  T_obs - tol)
  )
}


# ---- Independent-groups permutation workspace --------------------------
#
# Everything the reference distribution needs that does not depend on the
# candidate shift is computed once. 'perm' is an n_1 by m integer matrix
# whose columns are the observation indices assigned to the first group by
# each reassignment; 'mem' is a 0/1 indicator of membership in the observed
# first group, which is what a shift is applied to.
#
# The response is centered at its own mean before anything else. Both
# statistics offered by randomization_test() are invariant to adding a
# constant to every observation, so centering changes no result, but it
# keeps the sum of squares used by the studentized statistic away from the
# catastrophic cancellation that q - n * m^2 suffers when the mean is large
# relative to the spread.

.rt_ind_workspace <- function(y, mem, perm, statistic) {
  n_1 <- nrow(perm)
  N   <- length(y)
  y_c <- y - mean(y)
  list(
    statistic = statistic,
    n_1       = n_1,
    n_2       = N - n_1,
    y1        = matrix(y_c[perm], nrow = n_1),
    k1        = matrix(mem[perm], nrow = n_1),
    sum_y     = sum(y_c),
    sum_y2    = sum(y_c^2),
    sum_ym    = sum(y_c * mem)
  )
}


# ---- Statistics over the reassignments ---------------------------------
#
# Returns the test statistic for every column of the workspace's
# reassignment matrix, computed on the data shifted by 'delta' on the
# observed first group. delta = 0 gives the reference distribution for the
# null of no effect; a nonzero delta gives the reference distribution for
# the null that the treatment adds delta to every score, which is what the
# inverted confidence interval needs.
#
# The observed statistic is obtained by calling this on a workspace whose
# reassignment matrix has the single observed column, so the observed value
# and the reference distribution travel through identical arithmetic and
# ties between them are exact rather than approximate.

.rt_ind_null <- function(w, delta = 0) {
  z1 <- if (delta == 0) w$y1 else w$y1 - delta * w$k1
  s1 <- colSums(z1)
  s  <- w$sum_y - delta * w$n_1
  m_1 <- s1 / w$n_1
  m_2 <- (s - s1) / w$n_2

  if (w$statistic == "mean") return(m_1 - m_2)

  # Studentized (Welch) statistic. The total sum of squares of the shifted
  # data is available in closed form from the precomputed scalars, so only
  # the per-reassignment first-group sums of squares need the matrix.
  q1 <- colSums(z1 * z1)
  q  <- w$sum_y2 - 2 * delta * w$sum_ym + delta^2 * w$n_1
  v_1 <- (q1 - w$n_1 * m_1^2) / (w$n_1 - 1)
  v_2 <- ((q - q1) - w$n_2 * m_2^2) / (w$n_2 - 1)
  v_1[v_1 < 0] <- 0                    # rounding can push a zero variance below 0
  v_2[v_2 < 0] <- 0
  out <- (m_1 - m_2) / sqrt(v_1 / w$n_1 + v_2 / w$n_2)
  out[!is.finite(out)] <- 0            # a reassignment with no within-group spread
  out
}


# ---- Inverting the test for a confidence interval ----------------------
#
# Finds one endpoint of the set of shift values the test does not reject.
# 'p_fun' returns the randomization p-value for a candidate shift, 'center'
# is a shift the test certainly accepts (the observed mean difference, at
# which the shifted groups have identical means), 'side' is +1 for the
# upper endpoint and -1 for the lower, and 'scale' sets the first bracketing
# step.
#
# The p-value is a step function of the shift, so bisection rather than a
# derivative-based root finder is used: bisection only ever needs the sign
# of p_fun(delta) - alpha and keeps a valid bracket at every iteration. The
# returned endpoint is the last shift on the accepted side, so the reported
# interval is contained in the acceptance set rather than straddling it.
#
# An interval endpoint of Inf (or -Inf) is a real answer, not a failure: a
# design small enough that the smallest attainable p-value exceeds alpha
# cannot reject any shift, and its randomization interval is the whole line.

.rt_invert <- function(p_fun, alpha, center, scale, side, tol,
                       max_doubling = 60L, max_bisect = 200L) {
  if (p_fun(center) <= alpha) return(NA_real_)

  step <- scale
  lo   <- center
  hi   <- center + side * step
  k    <- 0L
  while (p_fun(hi) > alpha) {
    k <- k + 1L
    if (k > max_doubling) return(side * Inf)
    step <- step * 2
    hi   <- center + side * step
  }

  for (i in seq_len(max_bisect)) {
    if (abs(hi - lo) <= tol) break
    mid <- (lo + hi) / 2
    if (p_fun(mid) > alpha) lo <- mid else hi <- mid
  }
  lo
}
