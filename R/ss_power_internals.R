# Internal helpers shared by the single-degree-of-freedom ss_power_* planners
# whose test statistic is a noncentral t (a regression coefficient, a
# polynomial-change slope). Not exported.

# Power of a t test with noncentrality 'ncp' on 'df' degrees of freedom. For a
# nondirectional test the rejection region has two tails, at +/- the two-sided
# critical value, so both must be counted: summing only the upper tail
# understates the power and, at a null effect, returns alpha/2 instead of the
# nominal alpha. For a directional test the effect is taken in the direction
# of the test, so only the upper tail of |ncp| is counted. Callers pass ncp
# with the sign of the effect or its absolute value; the power is the same.
#
# The tails are those of .nc_t_power() (R/nc_t_internals.R): R's pt() inside
# the range in which it is accurate and the exact integral outside it, where
# pt() is a normal approximation (a noncentrality above 37.62, which a
# specified sample size or a large effect reaches). Through 3338679 each
# planner called pt() directly.
#
# ncp          noncentrality parameter (its sign is not used)
# df           error degrees of freedom
# alpha_level  per-test Type I error rate
# directional  TRUE for a one-sided test, FALSE (default) for two-sided
.power_noncentral_t <- function(ncp, df, alpha_level, directional = FALSE) {
  crit <- if (directional) stats::qt(1 - alpha_level, df = df) else
    stats::qt(1 - alpha_level / 2, df = df)
  .nc_t_power(crit, df, ncp, two_sided = !directional)
}

# Validate a fixed sample-size argument: a single finite whole number at or
# above 'min'. Returns it as a double, so that a product of two sizes
# cannot overflow. Returned as an integer, a size made such a product an
# integer, which overflowed to NA once it passed 2^31 - 1 (50,000 in each
# of two groups, or 40,000 and 60,000), and a size above 2^31 - 1 became
# NA itself. The specified-size paths of the ss_power_* planners otherwise
# returned NaN, zero, or a degenerate power table for a fractional,
# negative, or too-small size (for example one that leaves no residual
# degrees of freedom).
.check_whole_n <- function(x, arg, min) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < min)
    .dmar_refuse_argument(c(
      "'", arg, "' must be a single whole number of at least ", min,
      " (so the model has residual degrees of freedom)."), arg = NULL)
  as.numeric(x)
}
