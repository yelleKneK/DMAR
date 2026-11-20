#' @rdname ss_aipe_smd
#' @name ss_aipe_smd
#'
#' @title Sample Size Planning for the Standardized Mean Difference (AIPE)
#'
#' @description
#' Determines the per-group sample size needed for a two-independent-groups
#' design so that the confidence interval for the population standardized
#' mean difference, denoted \eqn{\delta}, is no wider than a user-specified
#' value when the sample standardized mean difference (Cohen's \emph{d})
#' equals \eqn{\delta}. This is the accuracy in parameter estimation (AIPE)
#' approach of Kelley and Rausch (2006), the standardized mean difference
#' companion to power-based planning via \code{\link{ss_power_smd}}.
#' Optionally, supplying \code{assurance} returns the smallest sample size,
#' not below the plan without it, at which the realized interval is no wider
#' than the target with at least that probability, computed exactly.
#'
#' @param true_smd The supposed value of the population standardized mean
#'   difference \eqn{\delta} the sample size is planned against: a value the
#'   researcher posits, either a minimally important effect or a value believed
#'   to be true in the population, never a sample estimate. Echoed in the
#'   returned table as the \code{true_smd} row.
#' @param conf_level Desired confidence level (i.e., \eqn{1-\alpha}, where
#'   \eqn{\alpha} is the Type I error rate). Default \code{0.95}. A level
#'   above .999999998, whose tail areas are below 1e-9, is refused (see
#'   "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param width Desired (full) width of the two-sided confidence interval
#'   on \eqn{\delta}.
#' @param assurance Optional probability with which the realized
#'   confidence interval is to be no wider than \code{width}. When
#'   \code{NULL} (the default), the plan is for the width of the interval
#'   at \eqn{d = \delta}; when supplied (e.g., 0.80, 0.90, 0.99), the
#'   function returns the smallest per-group sample size, not below that
#'   plan, at which the probability that the realized interval is no wider
#'   than \code{width} is at least \code{assurance} (see Details). Must be
#'   \code{NULL} or strictly between 0.50 and 1. An assurance above
#'   .999999998 is refused with a condition of class
#'   \code{dmar_argument_invalid}: the floor of 1e-9 on a tail area (see
#'   "Small tail areas" in \code{\link{ci_nc_t}}) applies to
#'   \code{(1 - assurance) / 2}.
#' @param assurance_method How the probability that the realized interval
#'   is no wider than \code{width} is found, used only with
#'   \code{assurance}. Only \code{"analytic"} (the default) is offered here:
#'   the probability is computed exactly, and the plan is the smallest
#'   sample size, not below the plan without \code{assurance}, at which it
#'   reaches \code{assurance}. \code{"simulation"}, which
#'   \code{\link{ss_aipe_sc_ancova}}, \code{\link{ss_aipe_reg_coef}} with
#'   \code{noncentral = TRUE}, and \code{\link{ss_aipe_src}} offer, where
#'   the probability has no exact form, is refused with a condition of class
#'   \code{dmar_argument_invalid}, whether or not \code{assurance} is given,
#'   since a simulated estimate of an exact probability would only add Monte
#'   Carlo error.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}:
#' \describe{
#'   \item{\code{necessary_n_per_group}}{the necessary per-group sample
#'     size.}
#'   \item{\code{true_smd}, \code{width}}{the planning inputs, so the
#'     assumptions the sample size was planned under travel with the
#'     result. The \code{true_smd} row is the supposed effect the plan is
#'     built on, as given (with its sign): a value the researcher posits,
#'     either a minimally important effect or a value believed to be true in
#'     the population, never a sample estimate.}
#'   \item{\code{assurance}}{the requested probability; present only when
#'     \code{assurance} was given.}
#'   \item{\code{assurance_attained}}{the exact probability, at
#'     \code{necessary_n_per_group} and \code{true_smd}, that the realized
#'     interval is no wider than \code{width}, at least \code{assurance};
#'     present only when \code{assurance} was given.}
#' }
#' The confidence level is reported in the printed footer, and with
#' \code{assurance} the footer names the assurance method; the attribute
#' \code{assurance_method} records it, and \code{assurance_probability}
#' records that \code{assurance_attained} is \code{"exact"}.
#'
#' @details
#' \strong{The plan without assurance.} The interval is the one
#' \code{\link{ci_smd}} reports: the confidence interval for the
#' noncentrality parameter of the noncentral \emph{t} distribution, with
#' \eqn{2n - 2} degrees of freedom, at the observed \emph{t} statistic,
#' rescaled to \eqn{\delta} by \eqn{\sqrt{2/n}}, with \emph{n} the
#' per-group sample size. The function uses \eqn{\delta} as if it were
#' \emph{d} and returns the smallest \emph{n} at which that interval is no
#' wider than \code{width} (Kelley & Rausch, 2006, pp. 371--372). The width
#' depends on \eqn{\delta} only through \eqn{|\delta|}.
#'
#' \strong{The assurance step.} At a given \emph{n} the width of the
#' interval depends on the sample only through \emph{d} and grows as
#' \eqn{|d|} grows, so the realized interval is no wider than the target
#' exactly when \eqn{|d| \le d^*(n)}, the value of \eqn{|d|} at which the
#' width equals the target. With \eqn{\lambda = |\delta| \sqrt{n/2}},
#' \eqn{\nu = 2n - 2}, and \eqn{t^* = d^*(n) \sqrt{n/2}}, the probability
#' of that event under normality with equal variances is
#' \deqn{P_n \;=\; \int_0^\infty \Big\{
#'   \Phi\big(t^* \sqrt{v/\nu} - \lambda\big) -
#'   \Phi\big(-t^* \sqrt{v/\nu} - \lambda\big) \Big\}\,
#'   f_\nu(v)\, dv,}
#' where \eqn{\Phi} is the standard normal distribution function and
#' \eqn{f_\nu} the density of a chi square on \eqn{\nu} degrees of
#' freedom: the \emph{t} statistic is \eqn{(Z + \lambda) / \sqrt{V/\nu}}
#' with \eqn{Z} standard normal and \eqn{V} chi square, independent. The
#' integral involves only the normal distribution and the chi square, so
#' it does not depend on R's noncentral \emph{t}, which is a normal
#' approximation above a noncentrality of 37.62. It depends on
#' \eqn{\delta} only through \eqn{|\delta|}, and the plan is the same for
#' \code{true_smd} and \code{-true_smd}.
#'
#' The returned sample size is the smallest \emph{n} at or above the plan
#' without assurance, \eqn{n_0}, at which \eqn{P_n} is at least
#' \code{assurance} (\code{assurance_method = "analytic"}), found by
#' doubling a step upward from \eqn{n_0} and then bisecting, and
#' \code{assurance_attained} is \eqn{P_n} there. The search assumes that
#' \eqn{P_n} does not decrease as \emph{n} grows from \eqn{n_0}. In 45
#' planning sets (\eqn{|\delta|} of 0, .2, .5, 1, and 2, widths of .25, .5,
#' and 1, and assurances of .80, .90, and .99) it did not decrease at any
#' \emph{n} from \eqn{n_0} to ten past the plan, and the plan was the
#' smallest \emph{n} found by raising \emph{n} one at a time from
#' \eqn{n_0}. Because the width changes little with \emph{d}, \eqn{P_n}
#' rises steeply with \emph{n}: at \code{true_smd = .5} and a width of .5 it
#' is .838 at 129 and .916 at 130.
#'
#' The probability was checked against the noncentral \emph{t}
#' distribution computed by R's \code{pt()} where it is exact, against a
#' second integral that conditions on the normal variate instead of the chi
#' square, and against the proportion of intervals no wider than the target
#' among samples of two groups of normal data.
#'
#' \strong{Large noncentrality.} R's \code{pt()} is a normal approximation
#' to the noncentral \emph{t} distribution above a noncentrality of 37.62,
#' and its series loses accuracy above 33 with more than 2,000 degrees of
#' freedom, ranges that narrow targets with a large \eqn{|\delta|} reach
#' (\code{true_smd = 2} with a width of .25, or \code{true_smd = 1} with a
#' width of .1). The interval \code{\link{ci_smd}} reports is computed
#' exactly there (see \code{\link{ci_nc_t}}), so the width of the interval
#' increases with \eqn{|d|} at every sample size, the event above is
#' \eqn{|d| \le d^*(n)}, and \eqn{P_n} is exact.
#'
#' \strong{The published procedure.} Kelley and Rausch (2006, pp.
#' 372--374) plan with \code{assurance} as follows: at \eqn{n_0}, find the
#' value \eqn{\delta_\gamma} that \eqn{|d|} exceeds with probability
#' \eqn{1 - \gamma}, where \eqn{\gamma} is \code{assurance}, by searching
#' between the upper limit of a one-sided interval and the upper limit of a
#' two-sided interval computed at \eqn{d = \delta} (their Equations 25 to
#' 31), and plan for \eqn{\delta_\gamma} as for \eqn{\delta}. The procedure
#' does not compute \eqn{P_n}. Since \eqn{\delta_\gamma} is found at
#' \eqn{n_0} and the plan for it is larger, that plan can be larger than
#' the smallest sample size that attains the request: in the 45 sets above
#' it attained the request in every set and was one case larger in 7, all
#' but one at an assurance of .99; at \code{true_smd = 1}, a width of .5,
#' and an assurance of .99 it is 150, where \eqn{P_n} is .9943, and this
#' function's plan is 149, where it is .9903. It can also fall short: at
#' \code{true_smd = 0}, a width of .8, and an assurance of .95 it is 49 per
#' group, where \eqn{P_n} is .9499, and this function's plan is 50, where
#' it is .9946.
#'
#' @references
#' Anderson, S. F., & Kelley, K. (2024). Sample size planning for
#'   replication studies: The devil is in the design. \emph{Psychological
#'   Methods, 29}(5), 844--867. \doi{10.1037/met0000520}
#'
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
#'   planning for more accurate statistical power: A method adjusting
#'   sample effect sizes for publication bias and uncertainty.
#'   \emph{Psychological Science, 28}(11), 1547--1562.
#'   \doi{10.1177/0956797617723724}
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
#' difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample-size planning.
#'   \emph{Evaluation & the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @section Warning: The returned value is the sample size \emph{per group}.
#'
#' @seealso
#' \code{\link{smd}}, \code{\link{smd_c}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, 
#' \code{\link{ci_nc_t}}, \code{\link[stats:power.t.test]{stats::power.t.test()}}
#'
#' @examples
#' ss_aipe_smd(true_smd = .5, conf_level = .95, width = .30)
#'
#' # With assurance the last row is the exact probability, at the planned
#' # sample size, that the interval is no wider than the target.
#' ss_aipe_smd(true_smd = .5, conf_level = .95, width = .30, assurance = .8)
#' ss_aipe_smd(true_smd = .5, conf_level = .95, width = .30, assurance = .95)
#'
#' # Reversing the groups reverses the sign of the effect and leaves the
#' # plan unchanged.
#' ss_aipe_smd(true_smd = -.5, conf_level = .95, width = .30, assurance = .95)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_smd <- function(true_smd, conf_level = 0.95, width, assurance = NULL,
                        assurance_method = "analytic") {
  .dmar_check_assurance_single(assurance)
  # The probability is exact here, so "simulation" is not offered
  # (R/ss_aipe_assurance_internals.R).
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_smd", simulation = FALSE, both = FALSE)
  # Validate the planning targets at entry so a boundary input errors clearly
  # rather than sending qnorm()/qt() out of domain and crashing the fixed-point
  # search with a missing value.
  if (missing(true_smd) || !is.numeric(true_smd) || length(true_smd) != 1L || !is.finite(true_smd))
    .dmar_refuse_argument(
      "'true_smd' must be a single finite number.", arg = "true_smd")
  if (missing(width) || !is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single finite positive number.", arg = "width")
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument(
      "'conf_level' must be a single number strictly between 0 and 1.",
      arg = "conf_level")
  if (!is.null(assurance)) {
    if (!is.finite(assurance) || (assurance <= 0) || (assurance >= 1)) .dmar_refuse_argument(
      "The 'assurance' must either be NULL or some value greater than zero and less than unity.",
      arg = "assurance")
    if (assurance <= .50) .dmar_refuse_argument(
      "The 'assurance' should be > 0.50 (but less than 1).", arg = "assurance")
    # The floor on a tail area (R/dmar_conditions.R) applies to
    # (1 - assurance) / 2; an assurance that would need a smaller tail is
    # refused in its own name.
    .dmar_check_assurance_tail_floor(assurance)
  }

  # The width of the interval depends on the estimate only through its
  # magnitude, and the distribution of |d| only on |true_smd|, so the plan
  # is made for |true_smd|; the signed value is echoed in the result. (DMAR
  # 1.0.x worked on the signed value in the assurance step, and reversing
  # the groups changed the plan.)
  delta <- abs(true_smd)

  n_0 <- .smd_expected_width_n(delta, conf_level, width)
  n <- n_0
  attained <- NULL

  # The assurance step: the smallest sample size at or above n_0, the plan
  # without assurance, at which the exact probability that the realized
  # interval is no wider than `width` (.smd_width_probability()) is at least
  # `assurance` (.smd_assurance_n()), as ss_aipe_partial_r() and ss_aipe_r()
  # plan.
  if (!is.null(assurance)) {
    step <- .smd_assurance_n(n_0, delta, width, conf_level, assurance)
    n <- step$n
    attained <- step$attained
  }

  # The necessary per-group sample size comes first, followed by rows that echo
  # the user-supplied planning inputs, so the assumptions the sample size was
  # planned under travel with the result. true_smd is the supposed effect
  # the plan is built on, a value the researcher posits (a minimally important
  # effect or a value believed to be true in the population), never a sample
  # estimate. The confidence level is reported in the printed footer; the value
  # column stays numeric.
  terms  <- c("necessary_n_per_group", "true_smd", "width")
  values <- c(n, true_smd, width)
  if (!is.null(assurance)) {
    terms  <- c(terms, "assurance", "assurance_attained")
    values <- c(values, assurance, attained)
  }
  out <- .as_dmar_tbl(data.frame(term = terms, value = values),
                      conf_level = conf_level, subclass = "dmar_ss_aipe")
  if (!is.null(assurance))
    out <- .ss_aipe_assurance_record(out, assurance_method, "exact")
  out
}

# ---- Internal helpers of ss_aipe_smd() -----------------------------------
#
# Notation: n per group, df = 2n - 2, the observed noncentrality
# t = d sqrt(n / 2), and the population noncentrality
# lambda = |delta| sqrt(n / 2). The interval is the one ci_smd() reports.

# The interval, and every noncentral t probability below, is exact: ci_nc_t()
# and .nc_t_prob() use R's pt() only where it is accurate and the integral
# of R/nc_t_internals.R elsewhere.

# .smd_interval(t, n, conf_level)
#
# The width, on the scale of d, of the interval ci_smd() reports at the
# observed noncentrality t with n per group.
.smd_interval <- function(t, n, conf_level) {
  ci <- ci_smd(ncp = t, n_1 = n, n_2 = n, conf_level = conf_level)
  ci$value[ci$term == "ci_upper"] - ci$value[ci$term == "ci_lower"]
}

# .smd_expected_width_n(delta, conf_level, width)
#
# The plan without assurance: the smallest n per group, at least 4, at
# which the interval computed at d = delta is no wider than `width`
# (Kelley & Rausch, 2006). The fixed point of the central t gives a start a
# few below the answer, and n then steps upward.
.smd_expected_width_n <- function(delta, conf_level, width) {
  alpha <- 1 - conf_level
  # An explicit cap on each fixed-point search so a pathological input cannot
  # loop without terminating.
  max_iter <- 100000L
  crit <- function(n) stats::qt(1 - alpha / 2, max(2 * n - 2, 1))

  # Initial starting value for n using the z distribution.
  n_0 <- 2 * (stats::qnorm(1 - alpha / 2) / (width / 2))^2
  # Second starting value for n using the central t distribution. Floor the
  # error df at 1 so that a very wide target (whose z-based start can be below
  # 1 per group) never passes a nonpositive df to qt() and yields NaN; the
  # final n is floored at 4 below regardless.
  n <- 2 * (crit(n_0) / (width / 2))^2
  difference <- abs(n - n_0)
  iter <- 0L
  while (difference > .000001) {
    iter <- iter + 1L
    if (iter > max_iter)
      .dmar_refuse(c(
        "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
        "computation_failed")
    n_p <- n
    n <- 2 * (crit(n) / (width / 2))^2
    difference <- abs(n - n_p)
  }
  # Four per group is the admissible minimum this planner returns, so an
  # impossibly wide target settles here rather than at a degenerate value.
  n <- max(4, ceiling(n) - 5)

  iter <- 0L
  while (.smd_interval(delta * sqrt(n / 2), n, conf_level) > width) {
    iter <- iter + 1L
    if (iter > max_iter)
      .dmar_refuse(c(
        "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
        "computation_failed")
    n <- n + 1
  }
  n
}

# The width probability. For the observed noncentrality t of two groups of
# n normal observations each with population standardized mean difference
# delta, t = (Z + lambda) / S with S = sqrt(V / df), df = 2n - 2, and
# lambda = |delta| sqrt(n / 2), so P(|d| <= k / sqrt(n / 2)) = P(|t| <= k)
# is .nc_t_abs_probability(k, 2n - 2, lambda), the integral over the chi
# square of {Phi(k s - lambda) - Phi(-k s - lambda)}.

# .smd_width_probability(n, delta, width, conf_level)
#
# The probability, at n per group and population value delta, that the
# interval ci_smd() reports is no wider than `width`. The width increases
# with |t|, so the event is |t| <= t*, with t* the observed noncentrality
# at which the width equals the target, and the probability is
# P(|t| <= t*) at the noncentrality |delta| sqrt(n / 2). When the interval
# at t = 0, the narrowest it can be, is already wider than the target, it
# is 0. (Before the interval was exact everywhere, the width it reported
# dropped by about 1e-4 where its upper limit passed 37.62, the event
# could differ from |t| <= t*, and the probability was reported as NA when
# that mattered.)
.smd_width_probability <- function(n, delta, width, conf_level) {
  at <- function(t) .smd_interval(t, n, conf_level)
  t_star <- 0
  if (at(0) <= width) {
    hi <- max(1, abs(delta) * sqrt(n / 2))
    doublings <- 0L
    while (at(hi) <= width) {
      hi <- 2 * hi
      doublings <- doublings + 1L
      if (doublings > 60L)
        .dmar_refuse(c(
          "The search for the estimate at which the interval reaches ",
          "'width' did not end; the target may be infeasible."),
          "computation_failed")
    }
    t_star <- stats::uniroot(function(t) at(t) - width,
                             lower = 0, upper = hi, tol = 1e-10)$root
  }
  .nc_t_abs_probability(t_star, 2 * n - 2, abs(delta) * sqrt(n / 2))
}

# .smd_assurance_n(n_start, delta, width, conf_level, assurance)
#
# The smallest n per group at or above n_start at which
# .smd_width_probability() is at least `assurance`. The step upward doubles
# until the request is met, and a bisection then finds the first n that
# meets it. Returns list(n, attained).
.smd_assurance_n <- function(n_start, delta, width, conf_level,
                             assurance) {
  prob <- function(n) .smd_width_probability(n, delta, width, conf_level)
  lo <- max(4, n_start)
  p_lo <- prob(lo)
  if (p_lo >= assurance) return(list(n = lo, attained = p_lo))
  step <- 1
  repeat {
    hi <- lo + step
    p_hi <- prob(hi)
    if (p_hi >= assurance) break
    lo <- hi
    step <- 2 * step
    if (hi > 1e9)
      .dmar_refuse_argument(c(
        "The sample size search for 'assurance' passed 1e9 per group ",
        "without reaching the requested probability."), arg = "assurance")
  }
  while (hi - lo > 1) {
    mid <- lo + (hi - lo) %/% 2
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
