#' Sequential Sample Size for a Fixed-Width Interval for a Mean Under Two-Phase Sampling
#'
#' @description
#' Evaluates, at the data in hand, a sequential stopping rule for a
#' two-phase design: the proxy measurement is available on every unit of a
#' pool, the costly target measurement is added to units drawn at random
#' from the pool a few at a time, and sampling stops at the first number of
#' labeled units, at or after a minimum of \code{m0}, at which the
#' confidence interval from \code{\link{ci_two_phase}} has a half-width no
#' larger than \code{half_width}. Unlike \code{\link{ss_aipe_two_phase}},
#' the rule needs no planning values for the target's standard deviation or
#' its correlation with the proxy: the data supply both as they accrue. It
#' follows \code{\link{ss_seq_c}} in its stopping rule, in the
#' \code{half_width} and \code{alpha_level} arguments, and in the rows it
#' returns (with the sizes named for the labeled units), and like it keeps
#' no state: call it again after each batch. It differs from
#' \code{ss_seq_c()} in three ways: it takes the data rather than summary
#' statistics; it has no \code{quantile} argument, since its criterion is
#' the half-width of the interval \code{\link{ci_two_phase}} reports, which
#' uses the \emph{t} quantile (the default of \code{ss_seq_c()}); and its
#' minimum \code{m0} counts labeled units and defaults to 30 rather than
#' \code{ss_seq_c()}'s 10 per group, because the slope as well as the
#' variance must be estimated before the rule can be trusted.
#'
#' @param y The target measurement on the units labeled so far.
#' @param x The proxy measurement on the same units, aligned with \code{y}.
#' @param x_unlabeled The proxy measurement on the units of the pool not
#'   yet labeled.
#' @param half_width The target half-width \eqn{h} of the
#'   100(1 - 2\eqn{\alpha})\% confidence interval, in the units of the
#'   target. For a noninferiority or equivalence decision at bounds
#'   \eqn{\pm\delta}, \code{\link{ss_seq_c}} explains why \eqn{\delta/2} is
#'   the usual choice.
#' @param alpha_level One-sided rate per bound; the interval is at
#'   confidence level 1 - 2\eqn{\alpha}. Default \code{0.05} (a 90\%
#'   interval), as in \code{\link{ss_seq_c}}.
#' @param m0 The minimum number of labeled units before the rule may stop,
#'   default \code{30}. The early variance and slope estimates are unstable,
#'   and a rule allowed to stop on them stops too early too often.
#'
#' @details
#' \strong{The rule.} With \eqn{n} labeled units and a pool of \eqn{M},
#' sampling stops when \eqn{n \ge m_0} and
#' \deqn{t_{1-\alpha,\, n-2} \left\{\frac{S_e^2}{n} +
#'   \frac{\hat b^2 S_X^2}{M}\right\}^{1/2} \le h,}
#' the left side being the half-width of the regression interval of
#' \code{\link{ci_two_phase}} at confidence level 1 - 2\eqn{\alpha}
#' (\eqn{S_e^2} the residual mean square of the target on the proxy,
#' \eqn{\hat b} the slope, \eqn{S_X^2} the proxy's variance over the pool).
#' This carries the fixed-width sequential approach of Chow and Robbins
#' (1965), which \code{\link{ss_seq_c}} applies to a contrast, over to the
#' two-phase interval. Labeled units must be added at random from the units
#' not yet labeled, so that at every step they remain a simple random
#' subsample of the pool.
#' In an evaluation of a system whose outputs an automated judge labels, for
#' example, a human assessor labels outputs drawn at random from those the
#' judge has labeled until the rule says stop.
#'
#' \strong{The projection.} \code{n_labeled_projected} is the number of
#' labeled units at which the rule would be met if the current estimates
#' held, from the normal approximation
#' \eqn{n = S_e^2 / \{(h / z_{1-\alpha})^2 - \hat b^2 S_X^2 / M\}}, the same
#' device \code{\link{ss_seq_c}} uses: a guide to how much further to go,
#' not a commitment. When the current estimates say that labeling every
#' unit of the pool but one (the interval needs at least one unit with the
#' proxy alone) would not reach the target, the projection is \code{NA}
#' with a warning: the pool, not the labeling, is then the constraint.
#'
#' \strong{A constant target, or one the proxy reproduces exactly.} While
#' every labeled unit has the same target value (for example, a binary
#' target that has been 1 on every unit so far), the residual variation and
#' the slope are both zero, and the interval would have zero width. While
#' the target is an exact linear function of the proxy on the labeled units
#' (for example, a binary judge and a binary human label that have agreed
#' on every unit so far), the residual variation is zero and the width of
#' the interval would come from the size of the pool alone. In both cases
#' the rule could stop on no evidence, so the function stops with an error
#' asking for more labeled units, whether or not the rule would have been
#' met, as it does while the proxy is constant on them. The residual
#' variance is taken to be zero when it is no larger than
#' \code{.Machine$double.eps} times the variance of the target on the
#' labeled units, which allows for rounding error.
#'
#' \strong{Coverage.} The stopping size is random and the interval at
#' stopping is computed as though it had been fixed, so its coverage is not
#' guaranteed in finite samples; it can fall short of nominal when \eqn{n}
#' at stopping is small, or for a bounded target (a proportion, for
#' example) whose variance depends on its mean. The minimum \code{m0} is
#' the guard against the worst of this, and a simulation of the whole
#' procedure in a plausible population is the check for a given design.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{stop} (1 when the rule is met, so labeling stops; 0 to continue),
#'   \code{half_width_current} (the half-width of the interval now),
#'   \code{half_width_target}, \code{n_labeled_current} (the labeled units
#'   used, after any with a missing value are set aside, as in
#'   \code{\link{ci_two_phase}}), \code{n_labeled_projected}, and
#'   \code{n_pool} (\eqn{M}).
#'
#' @references
#' Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of
#'   fixed-width sequential confidence intervals for the mean.
#'   \emph{The Annals of Mathematical Statistics, 36}(2), 457--462.
#'   \doi{10.1214/aoms/1177700156}
#'
#' Neyman, J. (1938). Contribution to the theory of sampling human
#'   populations. \emph{Journal of the American Statistical Association,
#'   33}(201), 101--116. \doi{10.1080/01621459.1938.10503378}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_seq_c}}, the sequential rule for a contrast;
#'   \code{\link{ci_two_phase}}; \code{\link{ss_aipe_two_phase}} for a
#'   design fixed in advance from planning values.
#'
#' @family sequential estimation
#'
#' @keywords design
#'
#' @examples
#' # A proxy is available on a pool of 5,000 units; the target has been
#' # measured on a random 40 of them so far. Is a 90% interval with
#' # half-width 1 reached?
#' set.seed(113)
#' proxy  <- rnorm(5000, mean = 50, sd = 10)
#' target <- 20 + 0.6 * proxy + rnorm(5000, sd = 5)
#' order  <- sample(5000)
#' first  <- order[1:40]
#' ss_seq_two_phase(y = target[first], x = proxy[first],
#'                  x_unlabeled = proxy[-first], half_width = 1)
#'
#' # Not yet; the projection suggests about 50. After labeling the next 30
#' # units at random:
#' first <- order[1:70]
#' ss_seq_two_phase(y = target[first], x = proxy[first],
#'                  x_unlabeled = proxy[-first], half_width = 1)
#'
#' @export
ss_seq_two_phase <- function(y, x, x_unlabeled, half_width,
                             alpha_level = 0.05, m0 = 30) {
  if (!is.numeric(half_width) || length(half_width) != 1L ||
      !is.finite(half_width) || half_width <= 0)
    .dmar_refuse_argument(
      "'half_width' must be a single positive number.", arg = "half_width")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 0.5)
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 0.5).", arg = "alpha_level")
  if (!is.numeric(m0) || length(m0) != 1L || !is.finite(m0) || m0 < 3 ||
      m0 != round(m0))
    .dmar_refuse_argument(
      "'m0' must be a single whole number of at least 3.", arg = "m0")

  conf_level <- 1 - 2 * alpha_level
  d <- .two_phase_prepare(y, x, x_unlabeled)
  xu <- d$x_unlabeled
  N <- length(xu)
  mean_u <- mean(xu)
  ss_u <- sum((xu - mean_u)^2)
  r <- .two_phase_regression(d$y, d$x, N, mean_u, ss_u)
  if (is.na(r[["slope"]]))
    .dmar_refuse_argument(c(
      "The proxy takes a single value on the labeled units, so the ",
      "two-phase interval is not yet defined; label more units."), arg = NULL)
  if (all(d$y == d$y[1L]))
    .dmar_refuse_argument(c(
      "The target takes a single value on the labeled units, so the ",
      "two-phase interval would have zero width and the rule could stop ",
      "on no evidence; label more units."), arg = NULL)
  # A target that is an exact linear function of the proxy on the labeled
  # units (a binary judge and a binary human label that agree on every unit
  # so far) has a residual variance of zero, or of rounding error. The
  # half-width would then be t * |b| * S_X / sqrt(M), the pool's part alone,
  # and the rule could stop on no evidence of how far the target departs
  # from the proxy. The refusal is made whether or not the rule would have
  # been met (below m0, or with a pool whose part alone exceeds the target),
  # as for a constant target. A residual standard deviation below about
  # 1.5e-8 of the target's is rounding error, not data.
  if (r[["var_residual"]] <= .Machine$double.eps * stats::var(d$y))
    .dmar_refuse_argument(c(
      "The target is an exact linear function of the proxy on the labeled ",
      "units, so the residual variance is zero and the two-phase interval ",
      "would reflect the size of the pool alone, with no evidence of how ",
      "far the target departs from the proxy; label more units."),
      arg = c("y", "x"))
  ci <- .ci_two_phase_fast(d$y, d$x, N, mean_u, ss_u, conf_level)
  n <- r[["n_labeled"]]
  M <- r[["n_pool"]]
  half_width_current <- (ci[["ci_upper"]] - ci[["ci_lower"]]) / 2
  stop_flag <- as.integer(n >= m0 && half_width_current <= half_width)

  # Projection at the current estimates with the normal quantile, as
  # ss_seq_c() projects: the df at the projected size is not yet known,
  # and the rule itself is re-evaluated as units accrue.
  z <- stats::qnorm(1 - alpha_level)
  room <- (half_width / z)^2 - r[["slope"]]^2 * r[["var_pool"]] / M
  n_projected <- if (room > 0) ceiling(r[["var_residual"]] / room) else Inf
  # The interval needs at least one unit with the proxy alone, so at most
  # M - 1 units can be labeled.
  if (n_projected > M - 1) {
    .dmar_warn(c(
      "At the current estimates, labeling every unit of the pool ",
      "but one (", M - 1, " of ", M, ") would not reach a half-width ",
      "of ", signif(half_width, 4), "; the projection is NA. A larger ",
      "pool or a wider target is needed unless the estimates change as ",
      "units accrue."), "result_qualified")
    n_projected <- NA_real_
  } else {
    n_projected <- max(n, n_projected)
  }

  out <- data.frame(
    term  = c("stop", "half_width_current", "half_width_target",
              "n_labeled_current", "n_labeled_projected", "n_pool"),
    value = c(stop_flag, half_width_current, half_width, n, n_projected, M),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
