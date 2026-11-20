# Internal helpers shared by the equivalence family.
#
# The verdict of an equivalence analysis of a contrast is read from its
# 100(1 - 2 alpha)% confidence interval against the bounds
# (-lower_bound, upper_bound), which equivalence_c() takes as
# contrast_lower and contrast_upper and plot_equivalence() as effect_lower
# and effect_upper. equivalence_c() attaches the verdict to its table and
# plot_equivalence() colors each interval by it. Both call
# .equivalence_verdict(), so the plot cannot disagree with the table.

# The six verdicts, in the order the labels are checked: a label applies only
# when none before it does, so the more informative claim wins. Equivalence
# implies both one-sided claims (noninferiority and nonsuperiority),
# superiority implies noninferiority, and inferiority implies nonsuperiority.
# What is left when neither one-sided claim holds is an interval that
# contains the whole region [-lower_bound, upper_bound]: inconclusive.
.EQUIVALENCE_VERDICTS <- c("Equivalent", "Superior", "Inferior",
                           "Noninferior only",
                           "Not relevantly superior only",
                           "Inconclusive")

# lower, upper : the confidence limits (vectors of one length).
# lower_bound  : the noninferiority margin as a positive magnitude, so the
#                lower bound of the region is -lower_bound, as
#                contrast_lower of equivalence_c() and effect_lower of
#                plot_equivalence() are given. This differs from the
#                arguments of the same name of power_equivalence_md()
#                and its planners, which are the signed ends of the region
#                (-1 and 1, or 0.8 and 1.25 on the ratio scale); the
#                internal never receives those.
# upper_bound  : the upper bound, positive.
# The callers pass the four by position. Returns a character vector of
# verdicts, one per interval, drawn from .EQUIVALENCE_VERDICTS. The
# inequalities are strict, as the tests are: a limit exactly on a bound
# establishes nothing about that bound.
.equivalence_verdict <- function(lower, upper, lower_bound, upper_bound) {
  noninferior <- lower > -lower_bound
  nonsuperior <- upper <  upper_bound
  ifelse(noninferior & nonsuperior, "Equivalent",
  ifelse(lower >  upper_bound,      "Superior",
  ifelse(upper < -lower_bound,      "Inferior",
  ifelse(noninferior,               "Noninferior only",
  ifelse(nonsuperior,               "Not relevantly superior only",
                                    "Inconclusive")))))
}


# The interval of equivalence_smd(), and the plan for it.
#
# equivalence_smd() reports, beside its two one-sided tests, the
# 100(1 - 2 alpha)% interval for the population standardized mean difference
#   d -/+ t_{1 - alpha, n_1 + n_2 - 2} SE(d),
#   SE(d) = sqrt((n_1 + n_2) / (n_1 n_2) + d^2 / (2 (n_1 + n_2))),
# the large-sample standard error of d with a Student t critical value. The
# interval lies inside (-smd_lower, smd_upper) exactly when both
# one-sided tests reject at alpha. equivalence_smd(), its planner
# ss_aipe_equivalence_smd(), and the simulation of
# ss_aipe_equivalence_smd_sensitivity() all form the interval here, so the
# interval a study is planned for is the interval the analysis reports and
# the interval the simulation scores.
#
# smd         : sample standardized mean difference(s), a numeric vector.
# n_1, n_2    : the two group sizes.
# alpha_level : the one-sided level; the interval is at 1 - 2 alpha_level.
# Returns a list of df, crit (the t quantile), se, lower, upper, and width
# (upper - lower, the width the analysis reports), the last four the length
# of smd.
.equivalence_smd_interval <- function(smd, n_1, n_2, alpha_level) {
  # Doubles, so the product n_1 * n_2 cannot overflow R's integer range
  # (it would once it passed 2^31 - 1, as at 46,341 per group, when the
  # sizes arrive as integers).
  n_1 <- as.numeric(n_1)
  n_2 <- as.numeric(n_2)
  df    <- n_1 + n_2 - 2
  se    <- sqrt((n_1 + n_2) / (n_1 * n_2) + smd^2 / (2 * (n_1 + n_2)))
  crit  <- stats::qt(1 - alpha_level, df)
  lower <- smd - crit * se
  upper <- smd + crit * se
  list(df = df, crit = crit, se = se, lower = lower, upper = upper,
       width = upper - lower)
}

# .equivalence_smd_plan_n(true_smd, width, alpha_level, n_max)
#
# The smallest per-group size n, at least 2 (the smallest equivalence_smd()
# accepts), at which the interval above, formed at a sample standardized mean
# difference equal to true_smd with n in each group, is no wider than
# `width`. With n per group the width is
#   t_{1 - alpha, 2n - 2} sqrt((8 + true_smd^2) / n),
# which decreases in n, so the answer is unique and at n - 1 the interval is
# wider. The search starts from z_{1 - alpha}^2 (8 + true_smd^2) / width^2,
# the answer when the normal quantile stands in for the t quantile; the t
# quantile is the larger, so the answer is not below it. A plan above
# n_max per group is refused, as it was when the search counted up from 4.
.equivalence_smd_plan_n <- function(true_smd, width, alpha_level,
                                    n_max = 1e6) {
  w <- function(n) .equivalence_smd_interval(true_smd, n, n, alpha_level)$width
  n <- max(2, floor(stats::qnorm(1 - alpha_level)^2 * (8 + true_smd^2) /
                      width^2))
  if (n > n_max) .equivalence_smd_refuse_n_max(FALSE)
  # The start is below the answer in exact arithmetic; the step down guards
  # against the rounding of the start only.
  while (n > 2 && w(n - 1) <= width) n <- n - 1
  while (w(n) > width) {
    n <- n + 1
    if (n > n_max) .equivalence_smd_refuse_n_max(FALSE)
  }
  n
}

# .equivalence_smd_check_smd(x, arg)
#
# The check of a population or planning standardized mean difference in
# ss_aipe_equivalence_smd() and ss_aipe_equivalence_smd_sensitivity(): a
# single finite number. .dmar_check_numbers() runs first and refuses NA,
# NaN, a character string, a list, and an empty vector; this check refuses
# what is left: a vector of more than one value, a logical, and an
# infinite value. Unchecked, an infinite planning value reaches the search
# and is refused as "Required n exceeded 1e6 per group", and an infinite or
# missing population value makes every simulated interval NaN. The refusal
# has class dmar_argument_invalid and names the argument as `arg`.
.equivalence_smd_check_smd <- function(x, arg) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x))
    return(invisible(x))
  .dmar_refuse_argument(c(
    "'", arg, "' must be a single finite number, a standardized mean ",
    "difference; ", .dmar_given(x), "."), arg = arg, value = x)
}

.equivalence_smd_refuse_n_max <- function(assurance) {
  .dmar_refuse_argument(
    if (assurance) "Required n exceeded 1e6 per group under assurance."
    else "Required n exceeded 1e6 per group; check inputs.", arg = NULL)
}

# .equivalence_smd_width_prob(n, true_smd, width, alpha_level)
#
# The probability that the interval of equivalence_smd(), from two
# independent normal samples of n each with a common variance and
# population standardized mean difference true_smd, is no wider than
# `width`. At a given n the width depends on the sample only through d and
# grows with |d|, so the interval is no wider than `width` exactly when
#   d^2 <= c^2 = n width^2 / t^2 - 8,   t = t_{1 - alpha, 2n - 2}.
# When c^2 <= 0 the interval at d = 0, the narrowest it can be, is already
# wider, and the probability is 0. Otherwise it is P(|T| <= c sqrt(n / 2)),
# where T = d sqrt(n / 2) is the two-group t statistic, noncentral t on
# 2n - 2 degrees of freedom with noncentrality true_smd sqrt(n / 2): the
# form of Kelley and Rausch (2006, Equation 20) for this interval. The
# probability is the integral over the chi square of
# .nc_t_abs_probability() (R/nc_t_internals.R), which does not use R's
# pt(), a normal approximation above a noncentrality of 37.62 that plans
# for narrow intervals at a large standardized mean difference reach.
.equivalence_smd_width_prob <- function(n, true_smd, width, alpha_level) {
  df <- 2 * n - 2
  t_crit <- stats::qt(1 - alpha_level, df)
  c2 <- n * width^2 / t_crit^2 - 8
  if (c2 <= 0) return(0)
  .nc_t_abs_probability(sqrt(c2) * sqrt(n / 2), df, true_smd * sqrt(n / 2))
}

# .equivalence_smd_assurance_n(n0, true_smd, width, alpha_level, assurance,
#                              n_max)
#
# The smallest n at or above n0, the plan without assurance, at which
# .equivalence_smd_width_prob() is at least `assurance`, with that
# probability. The step from n0 doubles until the probability is reached,
# and the last two sizes are bisected. Above n0 the probability increases
# with n: the bound c sqrt(n / 2) on |T| grows like n while the noncentrality
# grows like sqrt(n) and the spread of T does not grow. A plan above n_max
# per group is refused.
.equivalence_smd_assurance_n <- function(n0, true_smd, width, alpha_level,
                                         assurance, n_max = 1e6) {
  prob <- function(n) .equivalence_smd_width_prob(n, true_smd, width,
                                                  alpha_level)
  p_hi <- prob(n0)
  if (p_hi >= assurance) return(list(n = n0, attained = p_hi))
  lo <- n0
  step <- 1
  repeat {
    hi <- min(n0 + step, n_max)
    p_hi <- prob(hi)
    if (p_hi >= assurance) break
    if (hi >= n_max) .equivalence_smd_refuse_n_max(TRUE)
    lo <- hi
    step <- 2 * step
  }
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    p_mid <- prob(mid)
    if (p_mid >= assurance) {
      hi <- mid
      p_hi <- p_mid
    } else {
      lo <- mid
    }
  }
  list(n = hi, attained = p_hi)
}
