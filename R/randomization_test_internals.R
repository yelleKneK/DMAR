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
# .rt_is_extreme() marks the members of the reference distribution T_null
# that are at least as extreme as the observed statistic T_obs, in the
# direction named by 'alternative'; .rt_extreme() counts them. A
# reassignment that is mathematically tied with the observed one can miss
# equality by the rounding of the arithmetic that produced the two values,
# and a tie counts as at least as extreme, so the comparison allows 'tol'.
#
# 'tol' is the rounding bound of the statistic being compared, supplied by
# the caller, and it has no default. It is a bound on the error of the
# arithmetic, so it scales with the scores that entered the statistic and
# not with the statistic itself or with the largest reference value: the
# difference of two means of scores near a million carries the rounding of
# the scores however small the difference is. Because every bound used here
# is a multiple of the size of the scores, multiplying the response by
# c > 0 multiplies the statistic, the reference distribution, and the
# tolerance by c, and the set of reassignments counted does not change.
# 'tol' is a scalar or a vector with one entry per reference value. An
# entry that is not finite is taken as 0, so that an infinite statistic is
# compared directly.

.rt_is_extreme <- function(T_null, T_obs, alternative, tol) {
  tol[!is.finite(tol)] <- 0
  switch(alternative,
    two_sided = abs(T_null) >= abs(T_obs) - tol,
    less      = T_null <=  T_obs + tol,
    greater   = T_null >=  T_obs - tol
  )
}

.rt_extreme <- function(T_null, T_obs, alternative, tol) {
  sum(.rt_is_extreme(T_null, T_obs, alternative, tol))
}

# Rounding bound for comparing two signed sums of the paired differences
# d = y - x. With u = eps / 2, each difference carries at most
# u (|x_i| + |y_i|) + u |d_i| from the rounding of the scores and of the
# subtraction, and a sum of n terms adds at most (n - 1) u sum|d_i|. Two
# sums are compared, which gives eps (sum(|x| + |y|) + n sum|d|); the
# function returns twice that.

.rt_paired_tol <- function(x, y, d) {
  2 * .Machine$double.eps *
    (sum(abs(x) + abs(y)) + length(d) * sum(abs(d)))
}


# ---- Text of the refusal of an undefined studentized statistic ----------
#
# Both public functions refuse statistic = "t" when the observed standard
# error is 0 and the observed mean difference is not (see "When the
# studentized statistic is undefined" on either help page). The refusal
# states the exact p-value of the mean statistic as a fraction of the
# number of arrangements, and these three helpers format its parts.
#
# .rt_count_text() writes a number of arrangements in full, or, above
# 1e15, where a double no longer holds every whole number, as "about
# 1.18e+17". .rt_number_text() writes a score or a mean difference to
# seven significant digits. .rt_p_text() writes the fraction with its
# value, "2/70 = 0.02857", or "1" when the numerator is the whole
# reference set.

.rt_count_text <- function(count) {
  if (count > 1e15)
    return(paste("about", format(count, digits = 3, scientific = TRUE)))
  format(count, scientific = FALSE, trim = TRUE)
}

.rt_number_text <- function(x) format(x, digits = 7, trim = TRUE)

.rt_p_text <- function(numerator, count) {
  if (numerator >= count) return("1")
  denominator <- if (count > 1e15)
    format(count, digits = 3, scientific = TRUE) else .rt_count_text(count)
  paste0(numerator, "/", denominator, " = ",
         format(numerator / count, digits = 4, trim = TRUE))
}

# The sentence the refusal adds when the call would have been run by Monte
# Carlo sampling: the run estimates the exact value quoted and cannot
# report less than 1 / (n_resamples + 1). 'what' is "30 pairs" or
# "60 scores"; 'forced' says that the user asked for sampling with
# exact = FALSE where the default would have enumerated.

.rt_monte_carlo_text <- function(what, forced, n_resamples) {
  paste0(" With ", if (forced) "exact = FALSE" else what,
         " the test is run by Monte Carlo sampling, which estimates this ",
         "value and cannot report less than 1/(n_resamples + 1) = ",
         format(1 / (n_resamples + 1), digits = 4, scientific = FALSE,
                trim = TRUE), ".")
}


# ---- Independent-groups permutation workspace --------------------------
#
# Everything the reference distribution needs that does not depend on the
# candidate shift is computed once. 'perm' is an n_1 by m integer matrix
# whose columns are the observation indices assigned to the first group by
# each reassignment; 'mem' is a 0/1 indicator of membership in the observed
# first group, which is what a shift is applied to.
#
# The scores of both groups are kept for every reassignment (y1 and y2,
# with the membership indicators k1 and k2), so that each group's mean and
# sum of squares come from its own scores. The second group's quantities
# are not recovered by subtraction from totals, and no variance is formed
# as a sum of squares less n times a squared mean: both cancel to rounding
# error when a group has little spread. The indices of the second group,
# the complement of each column of 'perm', are found in one vectorized
# step.
#
# The response is centered at its own mean before anything else. Both
# statistics offered by randomization_test() are invariant to adding a
# constant to every observation, so centering changes no result, and it
# keeps the scores small relative to their spread. max_abs_y, the largest
# absolute value of the response before centering, is kept because the
# rounding of the centering itself is proportional to it.

.rt_ind_workspace <- function(y, mem, perm, statistic) {
  n_1 <- nrow(perm)
  m   <- ncol(perm)
  N   <- length(y)
  n_2 <- N - n_1
  y_c <- y - mean(y)

  in_first <- matrix(FALSE, nrow = N, ncol = m)
  in_first[cbind(as.vector(perm), rep(seq_len(m), each = n_1))] <- TRUE
  comp <- matrix((which(!in_first) - 1L) %% N + 1L, nrow = n_2)

  mem_int <- as.integer(mem)
  list(
    statistic = statistic,
    n_1       = n_1,
    n_2       = n_2,
    y_c       = y_c,
    mem       = mem,
    y1        = matrix(y_c[perm], nrow = n_1),
    k1        = matrix(mem_int[perm], nrow = n_1),
    y2        = matrix(y_c[comp], nrow = n_2),
    k2        = matrix(mem_int[comp], nrow = n_2),
    max_abs_y = max(abs(y))
  )
}


# ---- Statistics over the reassignments ---------------------------------
#
# Returns, for every column of the workspace's reassignment matrix, the
# test statistic computed on the data shifted by 'delta' on the observed
# first group, with the rounding bound of each value: a list with elements
# 'statistic' and 'tol'. delta = 0 gives the reference distribution for
# the null of no effect; a nonzero delta gives the reference distribution
# for the null that the treatment adds delta to every score, which is what
# the inverted confidence interval needs.
#
# The observed statistic is obtained by calling this on a workspace whose
# reassignment matrix has the single observed column, so the observed value
# and the reference distribution travel through the same arithmetic. The
# tolerance for comparing a reference value with the observed one is the
# sum of the two bounds.
#
# The bounds. Write z for a centered score less delta times its
# membership indicator, N for the total sample size, and eps for
# .Machine$double.eps.
#   a = eps (max|y| + max|z| + |delta|) bounds the error one shifted score
#       carries from the centering and the shift.
#   b = 2a + N eps max|z| bounds the error of a group mean, of a
#       difference between two group means, and of a residual about a
#       group mean.
# Mean statistic: the bound of each value is b.
# Studentized statistic: a group whose scores span no more than 2a is
# constant, and its sum of squares is exactly 0; otherwise the sum of
# squares is taken from the residuals about the group's own mean. A mean
# difference D with |D| <= b is 0. With se the separate-variances standard
# error, the statistic is D / se; when se is 0 it is Inf with the sign of
# D, or 0 when D is 0 as well (every score equal). An infinite value is the
# order of a reassignment inside the reference set only: when the observed
# assignment itself has se = 0 and D not 0, randomization_test() refuses
# statistic = "t" before the reference distribution is formed. The error
# of se^2 is at
# most 4 b se + 2 b^2, from which, with r = b / se, the error of one
# statistic t is at most b / se + |t| (2r + r^2) + 8 eps |t|. A bound that
# is not finite is returned as 0, so that infinite statistics are compared
# directly.
#
# The studentized comparison is exact for the data on which it was checked
# against rational arithmetic (see tests/testthat/test-randomization_test.R)
# and errs toward counting a tie otherwise: two reassignments whose
# statistics differ by less than the bound are counted as tied, so the
# p-value can be above the exact one and is not below it.

.rt_ind_null <- function(w, delta = 0) {
  eps <- .Machine$double.eps
  N   <- w$n_1 + w$n_2
  if (delta == 0) {
    z1 <- w$y1
    z2 <- w$y2
    max_abs_z <- max(abs(w$y_c))
  } else {
    z1 <- w$y1 - delta * w$k1
    z2 <- w$y2 - delta * w$k2
    max_abs_z <- max(abs(w$y_c - delta * w$mem))
  }
  a <- eps * (w$max_abs_y + max_abs_z + abs(delta))
  b <- 2 * a + N * eps * max_abs_z

  m_1 <- colSums(z1) / w$n_1
  m_2 <- colSums(z2) / w$n_2
  D   <- m_1 - m_2

  if (w$statistic == "mean") return(list(statistic = D, tol = b))

  # The studentized statistic does not depend on the unit of the scores, and
  # its squared residuals underflow or overflow for scores near 1e-200 or
  # 1e200. The scores, their means, and the two bounds are divided by the
  # power of 2 at or below the largest absolute score, which is exact in
  # binary arithmetic, so no statistic changes for scores of ordinary size.
  if (is.finite(max_abs_z) && max_abs_z > 0) {
    unit <- 2^floor(log2(max_abs_z))
    z1 <- z1 / unit
    z2 <- z2 / unit
    m_1 <- m_1 / unit
    m_2 <- m_2 / unit
    D <- D / unit
    a <- a / unit
    b <- b / unit
  }

  ss_1 <- .rt_group_ss(z1, m_1, 2 * a)
  ss_2 <- .rt_group_ss(z2, m_2, 2 * a)
  se <- sqrt(ss_1 / (w$n_1 * (w$n_1 - 1)) + ss_2 / (w$n_2 * (w$n_2 - 1)))

  D[abs(D) <= b] <- 0
  no_spread <- se == 0
  out <- D / se
  out[no_spread] <- sign(D[no_spread]) * Inf
  out[no_spread & D == 0] <- 0

  r   <- b / se
  tol <- b / se + abs(out) * (2 * r + r^2) + 8 * eps * abs(out)
  tol[!is.finite(tol)] <- 0
  list(statistic = out, tol = tol)
}

# Sum of squares of each column of 'z' about its own mean 'm', set to
# exactly 0 for a column whose range is at most 'constant_range', the
# rounding bound within which the scores of a group are one value.

.rt_group_ss <- function(z, m, constant_range) {
  hi <- z[1L, ]
  lo <- hi
  for (i in seq_len(nrow(z))[-1L]) {
    hi <- pmax(hi, z[i, ])
    lo <- pmin(lo, z[i, ])
  }
  res <- z - rep(m, each = nrow(z))
  ss  <- colSums(res * res)
  ss[hi - lo <= constant_range] <- 0
  ss
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
