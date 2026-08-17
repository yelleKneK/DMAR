# Internal helpers shared by the single-degree-of-freedom ss_power_* planners
# whose test statistic is a noncentral t (a regression coefficient, a
# polynomial-change slope). Not exported.

# Power of a t test with noncentrality 'ncp' on 'df' degrees of freedom. For a
# nondirectional test the rejection region has two tails, at +/- the two-sided
# critical value, so both must be counted: summing only the upper tail
# understates the power and, at a null effect, returns alpha/2 instead of the
# nominal alpha. Callers pass a nonnegative ncp (for example sqrt(N) * abs(f) or
# sqrt(lambda)); for a directional test the effect is taken in the direction of
# that ncp, so only the upper tail is counted.
#
# ncp          nonnegative noncentrality parameter
# df           error degrees of freedom
# alpha_level  per-test Type I error rate
# directional  TRUE for a one-sided test, FALSE (default) for two-sided
.power_noncentral_t <- function(ncp, df, alpha_level, directional = FALSE) {
  if (directional) {
    crit <- stats::qt(1 - alpha_level, df = df)
    return(stats::pt(crit, df = df, ncp = ncp, lower.tail = FALSE))
  }
  crit <- stats::qt(1 - alpha_level / 2, df = df)
  stats::pt(-crit, df = df, ncp = ncp) +
    stats::pt(crit, df = df, ncp = ncp, lower.tail = FALSE)
}

# Validate a fixed sample-size argument: a single finite whole number at or
# above 'min'. Returns it as an integer. The specified-size paths of the
# ss_power_* planners otherwise returned NaN, zero, or a degenerate power table
# for a fractional, negative, or too-small size (for example one that leaves no
# residual degrees of freedom).
.check_whole_n <- function(x, arg, min) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < min)
    stop("'", arg, "' must be a single whole number of at least ", min,
         " (so the model has residual degrees of freedom).", call. = FALSE)
  as.integer(x)
}
