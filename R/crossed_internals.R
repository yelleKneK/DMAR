# Shared numerics of the fully crossed units-by-raters design: the
# planner ss_aipe_crossed(), its Monte Carlo sibling
# ss_aipe_crossed_sensitivity(), and the analysis ci_crossed().
#
# The model is y_ur = mu + a_u + b_r + e_ur for unit u = 1, ..., n_u and
# rater r = 1, ..., n_r, one score per cell, with independent zero-mean
# unit effects (variance var_unit), rater effects (var_rater), and
# residuals (var_residual, which holds the unit-by-rater interaction a
# design with one score per cell cannot separate from error). Units and
# raters are both random, so the estimand is mu, the mean over the
# population of units and the population of raters, and the variance of
# the grand mean of a complete n_u x n_r array is
#
#   V = var_unit / n_u + var_rater / n_r + var_residual / (n_u n_r)
#     = {E(MS_u) + E(MS_r) - E(MS_e)} / (n_u n_r),
#
# with E(MS_u) = var_residual + n_r var_unit, E(MS_r) = var_residual +
# n_u var_rater, and E(MS_e) = var_residual. The analysis and the
# planner both rest on this identity: the analysis estimates V from the
# observed mean squares, and the planner evaluates it, and the
# Satterthwaite (1946) degrees of freedom of its estimate, from the
# planning values. None of these helpers is exported.


# Argument checks shared by the three exported functions. A number must
# be single and finite and lie in [lower, upper], or in (lower, upper)
# when `open`; the message says which, in words.
.crossed_check_number <- function(x, arg, lower = -Inf, upper = Inf,
                                  open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    (if (open) x > lower && x < upper else x >= lower && x <= upper)
  if (ok) return(invisible(x))
  what <- if (is.infinite(lower) && is.infinite(upper)) {
    "a single finite number"
  } else if (lower == 0 && is.infinite(upper)) {
    if (open) "a single positive number" else "a single nonnegative number"
  } else if (is.infinite(upper)) {
    sprintf(if (open) "a single number greater than %s" else
              "a single number of at least %s", format(lower))
  } else if (is.infinite(lower)) {
    sprintf(if (open) "a single number less than %s" else
              "a single number of at most %s", format(upper))
  } else if (open) {
    sprintf("a single number strictly between %s and %s",
            format(lower), format(upper))
  } else {
    sprintf("a single number between %s and %s, inclusive",
            format(lower), format(upper))
  }
  .dmar_refuse_argument(sprintf("'%s' must be %s.", arg, what), arg = arg)
}

# A count of units or raters: a single whole number of at least 2, the
# fewest that leave degrees of freedom for the mean squares.
.crossed_check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < 2)
    .dmar_refuse_argument(
      sprintf("'%s' must be a single whole number of at least 2.", arg),
      arg = arg)
  invisible(x)
}


# A count for a message: 100000 reads as "100,000", not "1e+05".
.crossed_format_count <- function(x)
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)


# The two-way ANOVA of a complete units x raters array (one score per
# cell). `x` is a numeric matrix with units in rows and raters in
# columns, already checked by the caller (complete, finite, at least two
# of each). The residual sum of squares is summed from the residuals
# themselves rather than formed as SS_total - SS_units - SS_raters, so it
# cannot come out negative through cancellation when the array is nearly
# additive.
#
# Returns a named numeric vector: grand_mean, ms_unit, ms_rater,
# ms_residual, n_units, n_raters.
.crossed_ms <- function(x) {
  n_u <- nrow(x)
  n_r <- ncol(x)
  g   <- mean(x)
  m_u <- rowMeans(x)
  m_r <- colMeans(x)
  ss_u <- n_r * sum((m_u - g)^2)
  ss_r <- n_u * sum((m_r - g)^2)
  ss_e <- sum((x - outer(m_u, m_r, "+") + g)^2)
  c(grand_mean  = g,
    ms_unit     = ss_u / (n_u - 1),
    ms_rater    = ss_r / (n_r - 1),
    ms_residual = ss_e / ((n_u - 1) * (n_r - 1)),
    n_units     = n_u,
    n_raters    = n_r)
}


# The interval for the grand mean of a complete crossed array, from its
# mean squares. Vectorized over its arguments, so the sensitivity sibling
# forms all G intervals in one call.
#
# The variance components are estimated by the ANOVA (method of moments)
# estimators truncated at zero, var_unit = max(0, (MS_u - MS_e) / n_r)
# and var_rater = max(0, (MS_r - MS_e) / n_u), with var_residual = MS_e,
# and the variance of the grand mean is estimated by plugging them into
# V. Written as a combination of mean squares, n_u n_r V-hat is
#
#   k_u MS_u + k_r MS_r + (1 - k_u - k_r) MS_e,
#
# where k_u (k_r) is 1 when MS_u (MS_r) exceeds MS_e and 0 otherwise. The
# combination is never below MS_e, so the estimate is positive whenever
# the array has any residual variation. Its degrees of freedom are
# Satterthwaite's (1946) for a linear combination of independent mean
# squares, sum(a_k MS_k)^2 / sum((a_k MS_k)^2 / df_k). The formula can be
# evaluated whatever the signs of the coefficients a_k, but Satterthwaite
# advised caution in using the approximation when some of them are
# negative, and here MS_e's coefficient is -1 whenever both components
# are kept. How well the approximation serves at a given design is what
# ss_aipe_crossed_sensitivity() lets a user check.
#
# Returns a list of vectors: grand_mean, se, df, ci_lower,
# ci_upper, var_unit, var_rater, var_residual.
.crossed_full_ci <- function(grand_mean, ms_unit, ms_rater, ms_residual,
                             n_units, n_raters, conf_level) {
  k_u <- as.numeric(ms_unit  > ms_residual)
  k_r <- as.numeric(ms_rater > ms_residual)
  k_e <- 1 - k_u - k_r
  comb <- k_u * ms_unit + k_r * ms_rater + k_e * ms_residual
  df <- comb^2 / ((k_u * ms_unit)^2 / (n_units - 1) +
                  (k_r * ms_rater)^2 / (n_raters - 1) +
                  (k_e * ms_residual)^2 / ((n_units - 1) * (n_raters - 1)))
  se <- sqrt(comb / (n_units * n_raters))
  q  <- stats::qt(1 - (1 - conf_level) / 2, df)
  list(grand_mean   = grand_mean,
       se           = se,
       df           = df,
       ci_lower  = grand_mean - q * se,
       ci_upper  = grand_mean + q * se,
       var_unit     = k_u * (ms_unit  - ms_residual) / n_raters,
       var_rater    = k_r * (ms_rater - ms_residual) / n_units,
       var_residual = ms_residual)
}


# The design-stage counterpart of .crossed_full_ci(): at n_u units and
# n_r raters, the variance of the grand mean implied by the planning
# values, and the Satterthwaite degrees of freedom its plug-in estimate
# (MS_u + MS_r - MS_e) / (n_u n_r) would have, computed from the
# expected mean squares. The planning degrees of freedom use the
# untruncated combination, the one whose expectation is n_u n_r V; when a
# planning component is near zero the analysis truncates its estimate
# often and its realized degrees of freedom differ from these, which is
# what ss_aipe_crossed_sensitivity() lets a user see. Vectorized over
# n_u and n_r.
.crossed_design <- function(n_u, n_r, var_unit, var_rater, var_residual) {
  e_u <- var_residual + n_r * var_unit
  e_r <- var_residual + n_u * var_rater
  e_e <- var_residual
  comb <- e_u + e_r - e_e
  list(variance = comb / (n_u * n_r),
       df = comb^2 / (e_u^2 / (n_u - 1) + e_r^2 / (n_r - 1) +
                        e_e^2 / ((n_u - 1) * (n_r - 1))))
}


# The half-width the plan attributes to the interval at (n_u, n_r).
# With the normal quantile it is z sqrt(V), the half-width when the
# variance of the mean is known. With the t quantile it is t_df sqrt(V),
# with the Satterthwaite df of .crossed_design(). With an assurance gamma
# it is the gamma quantile of the half-width under the approximation that
# df V-hat / V is chi square on df, t_df sqrt(V q_gamma / df), with
# q_gamma the gamma quantile of that chi square: the realized half-width
# is at most this value with probability gamma. That quantile is below
# the expected-width half-width t_df sqrt(V) whenever gamma is below
# P(chi square_df <= df), which exceeds one half (.58 at 5 df, .53 at 30)
# because the chi square's median lies below its mean. The plan without
# an assurance already delivers that probability, so the factor is
# floored at 1: an assurance plan is never smaller than the plan without
# one, and an assurance at or below one half returns that plan. The
# floor only raises the half-width, so the feasible designs with an
# assurance are a subset of those without. Vectorized over n_u and n_r.
.crossed_half_width <- function(n_u, n_r, var_unit, var_rater, var_residual,
                                conf_level, quantile, assurance) {
  d <- .crossed_design(n_u, n_r, var_unit, var_rater, var_residual)
  p <- 1 - (1 - conf_level) / 2
  if (quantile == "normal") return(stats::qnorm(p) * sqrt(d$variance))
  inflate <- if (is.null(assurance)) 1 else
    pmax(1, stats::qchisq(assurance, d$df) / d$df)
  stats::qt(p, d$df) * sqrt(d$variance * inflate)
}


# The smallest number of units, in [2, max_units], at which the planned
# half-width with n_r raters is at most h; NA where none is. Vectorized
# over n_r, with the searches for every n_r run side by side.
#
# With the normal quantile the planned half-width falls steadily as units
# are added. With the t quantile it does not: as n_u grows, V approaches
# its rater floor var_rater / n_r while the Satterthwaite degrees of
# freedom eventually fall toward n_r - 1, so past some number of units
# the larger quantile outweighs the smaller variance and the planned
# interval lengthens again. Over a wide grid of planning values, numbers of
# raters, confidence levels, and assurances, the planned half-width
# first falls and then (if at all) rises, never the other way, and the
# search relies on that shape. It first locates the turning point, the
# smallest n_u at which one more unit no longer shortens the interval
# (max_units when that never happens), by bisection on the sign of the
# forward difference. The number of raters is infeasible when even the
# turning point misses h; otherwise a second bisection finds the smallest
# feasible n_u on the falling stretch [2, turning point]. Each search
# takes about log2(max_units) evaluations, so the planner can afford to
# consider every number of raters up to max_raters.
.crossed_n_units <- function(n_r, h, max_units, half_width) {
  m  <- length(n_r)
  lo <- rep(2, m)
  hi <- rep(max_units, m)
  while (any(lo < hi)) {
    act  <- lo < hi
    mid  <- floor((lo + hi) / 2)
    rise <- half_width(mid + 1, n_r) >= half_width(mid, n_r)
    hi   <- ifelse(act & rise, mid, hi)
    lo   <- ifelse(act & !rise, mid + 1, lo)
  }
  turn <- lo
  ok_at_turn <- half_width(turn, n_r) <= h
  lo <- rep(2, m)
  hi <- turn
  while (any(ok_at_turn & lo < hi)) {
    act <- ok_at_turn & lo < hi
    mid <- floor((lo + hi) / 2)
    ok  <- half_width(mid, n_r) <= h
    hi  <- ifelse(act & ok, mid, hi)
    lo  <- ifelse(act & !ok, mid + 1, lo)
  }
  ifelse(ok_at_turn, lo, NA_real_)
}


# The planner's search: the least costly complete design whose planned
# half-width is at most width / 2. Every number of raters from 2 to
# max_raters (or the one fixed number n_raters) is paired with its
# smallest feasible number of units, and the cheapest pair is kept; cost
# rises with the number of units at a fixed number of raters, so nothing
# cheaper is missed. Among designs of equal cost the one with fewer
# raters is kept. The arguments are checked by the caller.
#
# Returns a list: n_units, n_raters, cost, df (all at the chosen design;
# NA throughout when no design within the limits meets the target), and
# the logical at_limit, TRUE when the cheapest design uses max_raters
# raters, so a larger limit might find a cheaper one. The planned
# half-width is not returned: the caller reports the expected width,
# which leaves out the assurance factor, and computes it itself.
.crossed_plan <- function(var_unit, var_rater, var_residual, width,
                          conf_level, assurance, n_raters, unit_cost,
                          rater_cost, judgment_cost, quantile, max_units,
                          max_raters) {
  hw <- function(n_u, n_r)
    .crossed_half_width(n_u, n_r, var_unit, var_rater, var_residual,
                        conf_level, quantile, assurance)
  n_r <- if (is.null(n_raters)) seq.int(2, max_raters) else n_raters
  n_u <- .crossed_n_units(n_r, width / 2, max_units, hw)
  cost <- unit_cost * n_u + rater_cost * n_r + judgment_cost * n_u * n_r
  if (all(is.na(cost)))
    return(list(n_units = NA_real_, n_raters = NA_real_, cost = NA_real_,
                df = NA_real_, at_limit = FALSE))
  i <- which.min(cost)
  list(n_units  = n_u[i],
       n_raters = n_r[i],
       cost     = cost[i],
       df       = .crossed_design(n_u[i], n_r[i], var_unit, var_rater,
                                  var_residual)$df,
       at_limit = is.null(n_raters) && n_r[i] == max_raters)
}


# The message for a target no design within the limits meets, naming the
# limit that binds. The search is repeated with room for far more units
# (a billion, or a thousand times max_units, whichever is larger, capped
# at 2^52 so that the bisection's midpoints stay whole numbers in double
# precision). If a design is then found, more units are the remedy,
# and the message gives the design. If none is, units cannot help with
# the raters allowed: a fixed panel needs more raters, and so does a
# search that already considered every number up to max_raters. The rater
# floor is named when it explains the shortfall, which is when the panel
# is at or below it; otherwise the t quantile's loss of degrees of
# freedom (or the assurance) does.
.crossed_infeasible_message <- function(var_unit, var_rater, var_residual,
                                        width, conf_level, assurance,
                                        n_raters, unit_cost, rater_cost,
                                        judgment_cost, quantile, max_units,
                                        max_raters, rater_floor) {
  fmt <- .crossed_format_count
  head <- if (is.null(n_raters)) {
    sprintf("No design with at most %s raters ('max_raters') and %s units ('max_units') has an interval as narrow as 'width'.",
            fmt(max_raters), fmt(max_units))
  } else {
    sprintf("No number of units up to %s ('max_units') gives an interval as narrow as 'width' with %s raters.",
            fmt(max_units), fmt(n_raters))
  }
  big <- min(max(1e9, 1000 * max_units), 2^52)
  wide <- .crossed_plan(var_unit, var_rater, var_residual, width,
                        conf_level, assurance, n_raters, unit_cost,
                        rater_cost, judgment_cost, quantile, big,
                        max_raters)
  if (!is.na(wide$n_units)) {
    tail <- if (is.null(n_raters)) {
      sprintf(" The least costly design without that limit on units has %s units and %s raters; a larger 'max_units' allows it.",
              fmt(wide$n_units), fmt(wide$n_raters))
    } else {
      sprintf(" With %s raters the target needs %s units; a larger 'max_units' allows them.",
              fmt(n_raters), fmt(wide$n_units))
    }
    return(paste0(head, tail))
  }
  panel <- if (is.null(n_raters)) max_raters else n_raters
  floor_text <- sprintf("However many units are added, the half-width with n raters stays above z * sqrt(var_rater / n), so the target needs more than %.4g raters (the rater floor) with the normal quantile, and more with the t quantile or an assurance.",
                        rater_floor)
  why <- if (panel <= rater_floor) {
    floor_text
  } else if (quantile == "normal") {
    sprintf("The panel is barely above the rater floor of %.4g raters, so not even %s units suffice.",
            rater_floor, fmt(big))
  } else {
    sprintf("The panel is above the rater floor of %.4g raters, but with the t quantile%s the planned interval stops narrowing as units are added: the Satterthwaite degrees of freedom fall toward the number of raters less one, and the larger quantile outweighs the smaller variance.",
            rater_floor, if (is.null(assurance)) "" else " and the assurance")
  }
  remedy <- if (is.null(n_raters)) {
    " More raters than 'max_raters' are needed, whatever the number of units."
  } else {
    " More raters are needed, whatever the number of units."
  }
  paste0(head, " ", why, remedy)
}
