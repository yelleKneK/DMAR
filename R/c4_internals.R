# ---------------------------------------------------------------------------
# The mean and variance of a standard deviation under normality, shared by
# ss_aipe_mean(), ss_aipe_two_phase(), design_consequences(), sd_unbiased(),
# var_ete(), and expected_r() (and through it expected_partial_r()), which
# uses the square of the mean factor and, for its bias, its logarithm.
# The correction factor J of Hedges (1981), used by smd(), smd_c(),
# meta_smd(), expected_smd(), and var_smd(), is formed from the same
# logarithm at the end of this file.
#
# When S^2 = sigma^2 W / df with W chi square on df degrees of freedom,
#   E(S) / sigma = E(sqrt(W / df)) = sqrt(2 / df) Gamma((df + 1) / 2) /
#                                    Gamma(df / 2),
# the mean of a chi variable over sqrt(df), and since E(W / df) = 1,
#   Var(S) / sigma^2 = 1 - {E(S) / sigma}^2.
# For a normal sample of size n, df = n - 1 and the mean factor is the
# constant c4(n) of quality control. It is below 1 at every df > 0
# (Jensen's inequality, since the square root is concave), it increases
# with df, and for a large df
#   log{E(S) / sigma} = -1/(4 df) + 1/(24 df^3) - 1/(20 df^5)
#                       + 17/(112 df^7) - 31/(36 df^9) + 691/(88 df^11) - ...,
# the asymptotic expansion of log Gamma(x + 1/2) - log Gamma(x) - log(x)/2
# at x = df / 2, whose terms are (2^-k - 2) B_(k+1) / {k (k + 1) x^k} for
# odd k, with B_j the Bernoulli numbers. So E(S) / sigma =
# 1 - 1/(4 df) + 1/(32 df^2) + O(df^-3) and Var(S) / sigma^2 =
# 1/(2 df) - 1/(8 df^2) + O(df^-3).
#
# Both factors are formed from L(df) = log{E(S) / sigma}, the mean factor
# as exp(L) and the variance factor as -expm1(2 L), so that neither takes
# a difference of nearly equal numbers. Two simpler forms do. A difference
# of lgamma() values, used through 1.0.x, keeps only the digits left after
# they cancel: for a large df each is about (df/2) log(df/2), and their
# difference about (1/2) log(df/2), so the mean factor so computed is off
# by about 4e-11 at df = 1e5, first reaches 1 at the integer df = 10810066,
# and is 1 or more at about half of all df above 1e8. And 1 minus the
# square of even a correctly rounded mean factor loses the leading digits
# the two share, about log10(2 df) of them. The ratio formed through
# lbeta(), as sqrt(pi) / B(df/2, 1/2), cancels nothing but is a rounded
# product of rounded factors: its relative error reaches 4e-15 (at
# df = 19.95), it is not monotone in df on a grid 1e-13 apart, and from
# about df = 1.7e14 on, where the factor is within 1.5e-15 of 1, it can
# come out a few units in the last place above 1.
#
# From df = 100 on, L is the expansion above through its df^-9 term; the
# first omitted term, 691/(88 df^11), is below 1e-21 there. Below df = 100,
# L is carried down from the expansion at df + 2m, the first such point at
# or above 100, by the exact recurrence
#   L(df) = L(df + 2) - log1p(1 / {df (df + 2)}) / 2,
# which follows from Gamma(x + 3/2) = (x + 1/2) Gamma(x + 1/2) and
# Gamma(x + 1) = x Gamma(x) at x = df / 2, since (x + 1/2)^2 =
# x (x + 1) + 1/4. Every term it adds is negative, as L is, so nothing
# cancels; they are summed from the smallest. L is negative at every
# df > 0, so the mean factor is never above 1 as computed, and it is
# exactly 1 only from about df = 4.5e15 on, where the true factor rounds
# to 1; the variance factor is above 0 at every df below about 4.5e307.
#
# Accuracy. Against a 320-bit evaluation, each is right to a few units in
# the last place. Over some 490,000 df from 0.05 to 1e15, most drawn at
# random (uniformly to 120, and on the log scale) and the rest on
# fixed grids (0.005 apart to 110, every integer to 1000, a
# two-thousandth of a decade apart from 1 to 1e15, and 1e-10 apart around
# every even df to 100), the largest relative errors found are 3.7e-16 for
# the mean factor, 5.5e-16 for the variance factor, and 4.7e-16 for L,
# about 1.7, 2.5, and 2.1 times the machine epsilon of 2.2e-16.
#
# Monotonicity. On grids 1e-13 apart across every even df from 2 to 100,
# where m changes or the expansion takes over, the mean factor and L never
# decrease and the variance factor never increases. From one double to the
# next they are not exactly monotone, as no value computed with a rounding
# error can be where the true value changes by far less than a unit in the
# last place: near df = 42 the mean factor changes by about 1e-18 from one
# double to the next, and it steps back by one unit at df = 42 + 1.13e-11;
# L and the variance factor step back by one unit from the double below 62
# to 62, where m changes. Over 500 consecutive doubles from each of 400
# points drawn at random from 0.5 to 99.5, and from each of 400 drawn on
# the log scale from 100 to 1e15, none of the three steps back; the
# lbeta() ratio steps back within those 500 doubles at every one of the
# 400 points below 100.
#
# All of these are vectorized, and give NaN at df <= 0.

# log{E(S) / sigma} from its expansion through the df^-9 term, in nested
# form: -1/(4 df) {1 - v/6 + v^2/5 - 17 v^3/28 + 31 v^4/9} with v = 1/df^2.
# Right to rounding for df >= 100.
.chi_log_mean_factor_series <- function(df) {
  v <- 1 / df^2
  -(1 / (4 * df)) *
    (1 - v * (1 / 6 - v * (1 / 5 - v * (17 / 28 - v * 31 / 9))))
}

# L(df) = log{E(S) / sigma} on df degrees of freedom. expected_r() uses it
# directly for its bias, 1 - (1 + t) exp(2 L), which at a large df is far
# smaller than either factor and keeps its digits only when formed from
# the logarithms.
.chi_log_mean_factor <- function(df) {
  out <- .chi_log_mean_factor_series(df)
  small <- which(df > 0 & df < 100)
  if (length(small) > 0L) {
    d <- df[small]
    m <- ceiling((100 - d) / 2)
    # Column j holds the recurrence step from d + 2 j + 2 down to d + 2 j,
    # for j < m; the columns run from the largest j, the smallest step.
    j <- seq.int(max(m) - 1, 0)
    dj <- outer(d, 2 * j, "+")
    step <- log1p(1 / (dj * (dj + 2)))
    step[outer(m, j, "<=")] <- 0
    out[small] <- .chi_log_mean_factor_series(d + 2 * m) - rowSums(step) / 2
  }
  out[which(df <= 0)] <- NaN
  out
}

# E(S) / sigma on df degrees of freedom.
.chi_mean_factor <- function(df) exp(.chi_log_mean_factor(df))

# c4(n) = E(S) / sigma for a normal sample of size n.
.c4 <- function(n) .chi_mean_factor(n - 1)

# Var(S) / sigma^2 = 1 - {E(S) / sigma}^2 on df degrees of freedom.
.chi_var_factor <- function(df) -expm1(2 * .chi_log_mean_factor(df))

# The correction factor of Hedges (1981) on df degrees of freedom,
#   J(df) = Gamma(df / 2) / {sqrt(df / 2) Gamma((df - 1) / 2)},
# the constant for which J(df) d is an unbiased estimator of the
# standardized mean difference. E(S) / sigma on df - 1 degrees of freedom
# is sqrt(2 / (df - 1)) Gamma(df / 2) / Gamma((df - 1) / 2), and
# multiplying it by sqrt((df - 1) / df) gives J(df), so
#   log J(df) = L(df - 1) + log1p(-1 / df) / 2,
# with L the logarithm of the mean factor above. Both terms are negative,
# so nothing cancels, and J is below 1 at every df as computed. The
# difference of lgamma() values used through 1.0.x keeps only the digits
# left after two numbers of about (df/2) log(df/2) cancel: its relative
# error in J is 5e-11 at df = 1e5 and 8e-9 at df = 1e7, and at df = 1e8
# it returned a J above 1. Callers that need 1/J - 1 or 1 - 1/J^2 take them
# from the logarithm with expm1(), as expected_smd() and var_smd() do.
# Against a 60-digit evaluation the relative error of J is at most
# 6.2e-17 and that of 1/J - 1 at most 1.2e-16 for df from 2 to 1e14.
#
# At df = 1 the factor is 0, the limit of the gamma ratio (the mean of a
# standard deviation on 0 degrees of freedom does not exist, and neither
# does the mean of d), which is the value the 1.0.x form returned; below
# df = 1 it is NaN. Vectorized.
.hedges_log_J <- function(df) {
  out <- suppressWarnings(.chi_log_mean_factor(df - 1) + log1p(-1 / df) / 2)
  out[which(df == 1)] <- -Inf
  out[which(df < 1)] <- NaN
  out
}

# J(df) itself.
.hedges_J <- function(df) exp(.hedges_log_J(df))

# .hedges_require_df(df, fun, sizes)
#
# The refusal smd() and smd_c() give when unbiased = TRUE is asked for on
# 1 degree of freedom (or fewer). J(1) is 0, so the corrected estimate
# would be 0 whatever the data, which is not an estimate: d is a multiple
# of a noncentral t variable, whose mean does not exist on 1 degree of
# freedom, so no multiple of d is unbiased there. 1.0.x returned the 0.
# On 0 degrees of freedom (a control group of 1 given to smd_c() with
# summary statistics) the reason is another one: J(df) is not defined at
# df = 0, and with no degrees of freedom there is no standard deviation to
# standardize by. 1.0.x returned Inf there (NaN when the means were equal), from
# exp(lgamma(0) - lgamma(-1 / 2)) / sqrt(0).
# `fun` is the calling function's name and `sizes` the clause that says
# which sample size gives the degrees of freedom. The condition has class
# dmar_hedges_correction_undefined and carries `df`.
.hedges_require_df <- function(df, fun, sizes) {
  if (!isTRUE(df <= 1)) return(invisible(TRUE))
  opening <- c(
    fun, "(unbiased = TRUE) is not available on ", format(df),
    if (isTRUE(df == 1)) " degree" else " degrees", " of freedom (", sizes,
    "). ")
  reason <- if (isTRUE(df == 1)) {
    c("The correction factor of Hedges (1981), J(df) = Gamma(df / 2) / ",
      "{sqrt(df / 2) Gamma((df - 1) / 2)}, is 0 at df = 1, so the corrected ",
      "estimate would be 0 whatever the means and the standard deviation ",
      "are. That is not an estimate of the standardized mean difference: ",
      "the mean of d does not exist on 1 degree of freedom, so no multiple ",
      "of d is unbiased there. ")
  } else {
    c("The correction factor of Hedges (1981), J(df) = Gamma(df / 2) / ",
      "{sqrt(df / 2) Gamma((df - 1) / 2)}, is not defined at df = 0: with ",
      "no degrees of freedom the observations give no standard deviation, ",
      "so there is nothing to standardize the mean difference by and ",
      "nothing to correct. ")
  }
  .dmar_refuse(c(opening, reason,
    "Use unbiased = FALSE for the uncorrected estimate, or supply data ",
    "with at least 2 degrees of freedom."),
    "hedges_correction_undefined", df = df)
}
