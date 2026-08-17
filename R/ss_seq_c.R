# Sequential fixed-width confidence interval for a linear contrast.
#' Sequential Sample Size for a Fixed-Width Contrast Interval
#'
#' Implements the purely sequential fixed-width confidence interval
#' procedure for a linear contrast \eqn{\psi = \sum_j c_j \mu_j},
#' following Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal
#' (2025). The goal is a 100(1 - 2\eqn{\alpha})\% confidence interval
#' \eqn{\hat\psi \pm h} whose half-width \eqn{h} is fixed in advance,
#' which is what makes a noninferiority or equivalence verdict
#' reachable by design: with bounds \eqn{\pm\delta}, equivalence can
#' never be declared unless \eqn{h < \delta}, and targeting
#' \eqn{h = \delta/2} gives a truly equivalent contrast about a 90\%
#' chance of being declared equivalent at \eqn{\alpha = .05}. Because
#' no fixed sample size can guarantee a bounded-width interval when
#' the error variance is unknown (Dantzig, 1940), the procedure is
#' sequential: begin with a pilot, then keep sampling, re-estimating
#' the variance, until the stopping criterion is met. Used with
#' \code{pilot = TRUE} the function plans the pilot and the
#' cost-optimal allocation; with \code{pilot = FALSE} it evaluates the
#' stopping criterion at the data in hand, in the style of
#' \code{\link{mr_smd}}.
#'
#' @param c_weights The contrast weights. The weights must sum to zero
#'   with the positive weights summing to 1 and the negative weights
#'   to -1, so that \code{half_width} is on the raw scale of the
#'   response.
#' @param half_width The target half-width \eqn{h} of the
#'   100(1 - 2\eqn{\alpha})\% confidence interval, in raw units of the
#'   response. For a noninferiority or equivalence decision at bounds
#'   \eqn{\pm\delta}, the recommended target is \eqn{\delta/2}.
#' @param s The current estimate(s) of the error standard deviation:
#'   either a single pooled value or one value per group (aligned with
#'   \code{c_weights}). Required when \code{pilot = FALSE}; optional
#'   planning values when \code{pilot = TRUE} (used only to shape the
#'   allocation).
#' @param n The current per-group sample sizes, aligned with
#'   \code{c_weights}. Required when \code{pilot = FALSE}.
#' @param cost Optional per-observation sampling costs, one per group,
#'   aligned with \code{c_weights}. When supplied, the reported
#'   allocation is the cost-optimal
#'   \eqn{n_j \propto |c_j|\, \sigma_j / \sqrt{\mathrm{cost}_j}}
#'   (Chattopadhyay et al., 2025); with equal costs and standard
#'   deviations this reduces to allocation proportional to
#'   \eqn{|c_j|}.
#' @param alpha_level One-sided rate per bound; the interval is at
#'   confidence level 1 - 2\eqn{\alpha}. Default \code{0.05} (a 90\%
#'   interval).
#' @param quantile \code{"t"} (default) uses the \emph{t} quantile on
#'   the current error degrees of freedom in the stopping criterion, a
#'   finite-sample refinement; \code{"normal"} uses the normal
#'   quantile of the classical statement of the rule (Chow & Robbins,
#'   1965). The normal rule stops slightly early at wide targets, the
#'   finite-sample undershoot anticipated by Woodroofe (1977); the
#'   \emph{t} rule corrects it at a small cost in sample size.
#' @param pilot \code{TRUE} to plan the pilot stage; \code{FALSE}
#'   (default) to evaluate the stopping criterion at the current data.
#' @param m0 The minimum pilot sample size per group. Default
#'   \code{10}. A pilot that is too small makes the variance estimate
#'   that drives the early steps unstable.
#'
#' @return With \code{pilot = TRUE}, a \code{data.frame} with row
#'   \code{pilot_n_per_group} (the pilot size for each group with a
#'   nonzero weight) followed by one \code{allocation_j} row per group
#'   giving the recommended sampling proportions for the accrual
#'   stage. With \code{pilot = FALSE}, rows \code{stop} (1 = the
#'   criterion is met, stop sampling; 0 = continue),
#'   \code{half_width_current} (the half-width the interval would have
#'   now), \code{half_width_target}, \code{N_current},
#'   \code{N_projected} (the approximate total at which the criterion
#'   would be met under the current allocation, from the normal
#'   approximation), and the \code{allocation_j} rows for the next
#'   round of sampling.
#'
#' @details
#' \strong{The stopping criterion.} Sampling stops at the first
#' \eqn{N} for which
#' \deqn{q^2 \sum_j c_j^2 s_j^2 / n_j \;\le\; h^2,}
#' where \eqn{q} is the \emph{t} or normal quantile at \eqn{1-\alpha}.
#' With a pooled \eqn{s} and equal allocation this is the Chow and
#' Robbins (1965) rule specialized to a contrast; the procedure is
#' asymptotically consistent (coverage approaches 1 - 2\eqn{\alpha})
#' and first-order efficient (the mean stopping size approaches the
#' oracle \eqn{n^*} an investigator with known variance would use).
#' See \code{\link{ss_seq_c_sensitivity}} for a Monte Carlo evaluation
#' of both properties.
#'
#' \strong{Degrees of freedom.} With a pooled \code{s} the criterion
#' uses \eqn{\nu = N - J}. With per-group \code{s} it uses the
#' Satterthwaite approximation, which is the appropriate error law
#' when the variances are not assumed homogeneous.
#'
#' \strong{Batches.} Observations may be added in batches rather than
#' one at a time; the asymptotic properties survive batching. Re-call
#' the function after each batch.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of
#'   fixed-width sequential confidence intervals for the mean.
#'   \emph{The Annals of Mathematical Statistics, 36}(2), 457--462.
#'
#' Dantzig, G. B. (1940). On the non-existence of tests of "Student's"
#'   hypothesis having power functions independent of \eqn{\sigma}.
#'   \emph{The Annals of Mathematical Statistics, 11}(2), 186--192.
#'
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential
#'   methods and their applications}. CRC Press.
#'
#' Woodroofe, M. (1977). Second order approximations for sequential
#'   point and interval estimation. \emph{The Annals of Statistics,
#'   5}(5), 984--995.
#'
#' @seealso \code{\link{ss_seq_c_sensitivity}}, \code{\link{ss_aipe_c}},
#'   \code{\link{ss_power_equivalence_c}}, \code{\link{equivalence_c}},
#'   \code{\link{mr_smd}}
#'
#' @examples
#' # 1. Plan the pilot for a two-group contrast, target half-width 2.5
#' #    (bounds of 5 with the h = delta/2 rule):
#' ss_seq_c(c_weights = c(1, -1), half_width = 2.5, pilot = TRUE)
#'
#' # 2. Evaluate the stopping criterion mid-study: pooled s = 15.4 at
#' #    n = 60 per group. Too imprecise to stop; the projection says
#' #    roughly how much further to go.
#' ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
#'          s = 15.4, n = c(60, 60))
#'
#' # 3. Costs differ: sampling the second group costs four times as
#' #    much per observation, so its share of new observations drops.
#' ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
#'          s = c(15.4, 15.4), n = c(60, 60), cost = c(1, 4))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family sequential estimation
#'
#' @export

ss_seq_c <- function(c_weights, half_width, s = NULL, n = NULL, cost = NULL,
                     alpha_level = 0.05, quantile = c("t", "normal"),
                     pilot = FALSE, m0 = 10) {
  quantile <- match.arg(quantile)

  if (!is.numeric(half_width) || length(half_width) != 1L || half_width <= 0)
    stop("'half_width' must be a single positive number.")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be a single number in (0, 0.5).")
  if (!identical(round(sum(c_weights), 5), 0))
    stop("The sum of the contrast weights ('c_weights') should equal zero.")
  pos <- sum(c_weights[c_weights > 0])
  neg <- sum(c_weights[c_weights < 0])
  if (!isTRUE(all.equal(pos, 1)) || !isTRUE(all.equal(neg, -1)))
    stop("The positive weights must sum to 1 and the negative weights to ",
         "-1, so that 'half_width' is on the raw scale of the response.")

  J        <- length(c_weights)
  involved <- which(c_weights != 0)

  if (!is.null(cost)) {
    if (length(cost) != J || any(cost <= 0))
      stop("'cost' must give one positive per-observation cost per group, ",
           "aligned with 'c_weights'.")
  } else {
    cost <- rep(1, J)
  }

  # Cost-optimal allocation over the involved groups:
  # n_j proportional to |c_j| * sigma_j / sqrt(cost_j). Planning s defaults
  # to equal across groups, in which case the shape comes from the weights
  # and costs alone.
  s_alloc <- if (is.null(s)) rep(1, J) else if (length(s) == 1L) rep(s, J) else s
  if (length(s_alloc) != J)
    stop("'s' must be a single pooled value or one value per group.")
  if (any(s_alloc[involved] <= 0))
    stop("'s' must be positive for every group with a nonzero weight.")
  raw_alloc <- abs(c_weights) * s_alloc / sqrt(cost)
  allocation <- raw_alloc / sum(raw_alloc)

  alloc_terms <- paste0("allocation_", seq_len(J))

  if (pilot) {
    if (!is.numeric(m0) || length(m0) != 1L || m0 < 2)
      stop("'m0' must be a single number of at least 2.")
    out <- data.frame(
      term  = c("pilot_n_per_group", alloc_terms),
      value = c(ceiling(m0), allocation),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    return(.as_dmar_tbl(out, conf_level = 1 - 2 * alpha_level))
  }

  if (is.null(s) || is.null(n))
    stop("With 'pilot = FALSE', supply the current 's' and 'n'.")
  if (length(n) != J)
    stop("'n' must give one current sample size per group, aligned with ",
         "'c_weights'.")
  if (any(n[involved] < 2))
    stop("Each group with a nonzero weight needs at least 2 observations.")

  s_j <- if (length(s) == 1L) rep(s, J) else s
  var_hat <- sum(c_weights^2 * s_j^2 / n)

  if (length(s) == 1L) {
    nu <- sum(n) - J
  } else {
    # Satterthwaite approximation for heterogeneous variances.
    terms <- c_weights^2 * s_j^2 / n
    nu <- var_hat^2 / sum(ifelse(c_weights != 0,
                                 terms^2 / (n - 1), 0))
  }
  if (nu <= 0) stop("The error degrees of freedom must be positive.")

  q <- if (quantile == "t") stats::qt(1 - alpha_level, df = nu) else
    stats::qnorm(1 - alpha_level)

  half_width_current <- q * sqrt(var_hat)
  stop_flag <- as.integer(half_width_current <= half_width)

  # Projection under the current allocation shape: scaling every n_j by f
  # scales the variance by 1/f, so the criterion is met at approximately
  # f = (half_width_current / half_width)^2 times the current total. The
  # normal quantile is used for the projection because the df at the
  # projected N is not yet known; the rule itself is re-evaluated as data
  # accrue, so the projection is a guide, not a commitment.
  z <- stats::qnorm(1 - alpha_level)
  f <- (z * sqrt(var_hat) / half_width)^2
  N_projected <- max(sum(n), ceiling(f * sum(n)))

  out <- data.frame(
    term  = c("stop", "half_width_current", "half_width_target",
              "N_current", "N_projected", alloc_terms),
    value = c(stop_flag, half_width_current, half_width,
              sum(n), N_projected, allocation),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = 1 - 2 * alpha_level)
}
