# Sensitivity analysis for AIPE on Cliff's delta.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for Cliff's Delta
#'
#' @description
#' Quantifies how much misspecification of the population Cliff's delta
#' (\eqn{\delta = \Pr(X > Y) - \Pr(X < Y)}) distorts an AIPE-based
#' sample size plan. On each replication the function simulates two
#' independent samples whose population Cliff's delta equals
#' \code{true_delta}, computes the sample \code{\link{cliff_delta}}
#' and its CI, and summarizes the realized widths and coverage.
#'
#' \strong{Replications without an interval.} When the two simulated
#' groups are completely separated (every score of one above every score
#' of the other), or every score is equal, the estimated variance of
#' Cliff's delta is 0 and \code{\link{cliff_delta}} reports both limits
#' equal to the estimate, with a warning that they are not a confidence
#' interval. Such a replication has no interval here. Its estimate still
#' enters \code{mean_cliff_delta}, \code{median_cliff_delta}, and
#' \code{sd_cliff_delta}. The width, width-target, and coverage rows
#' (\code{mean_ci_width}, \code{median_ci_width}, \code{sd_ci_width},
#' \code{pct_ci_less_width}, \code{pct_ci_miss_low}, \code{pct_ci_miss_high},
#' \code{total_type_I_error}) describe the replications that have an
#' interval, and are \code{NA} when none has; \code{pct_ci_less_width} is
#' the number of intervals no wider than \code{width} divided by the
#' number of replications with an interval. The warning of each
#' replication is not passed on. One warning of class
#' \code{dmar_replications_without_interval}, with the fields
#' \code{n_without_interval} and \code{G}, states how many there were,
#' and the \code{simulation} attribute records the counts. Separation is
#' common when \code{true_delta} is large in absolute value and the
#' groups are small: with \code{true_delta = 0.9} and six per group,
#' about four replications in ten are separated.
#'
#' \strong{Data generating mechanism.} The simulator draws each sample
#' from a normal distribution and chooses the mean shift so that the
#' implied Cliff's delta equals \code{true_delta}. For normal samples
#' \eqn{\delta = 2 \Phi(\Delta/\sqrt{2}) - 1} where \eqn{\Delta} is the
#' standardized mean difference, so the simulator sets
#' \eqn{\Delta = \sqrt{2} \cdot \Phi^{-1}((1 + \delta)/2)}.
#'
#' @param true_delta Population Cliff's delta (the data generating
#'   value); in \eqn{(-1, 1)}.
#' @param supposed_delta The supposed population Cliff's delta, the
#'   planning value passed to \code{\link{ss_aipe_cliff_delta}} (not an
#'   estimate computed from data); supply this or \code{specified_N} but
#'   not both.
#' @param ratio Allocation ratio \eqn{n_1 / n_2} (default 1).
#' @param width Desired full width of the CI on Cliff's delta.
#' @param specified_N Total sample size to evaluate (split per
#'   \code{ratio}).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample Cliff's delta, the two confidence limits, the
#'   interval width, and two indicators of whether the interval missed
#'   \code{true_delta} below or above): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a
#'   header row is created otherwise, an existing file at that path is
#'   appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}. A replication without an interval
#'   is recorded with its estimate and with \code{NA} in the other
#'   columns.
#' @param seed Optional integer seed. The default \code{NULL} draws from,
#'   and advances, the current random number stream; an integer makes the
#'   call reproducible, and the caller's random number state is restored
#'   when the call returns.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with columns
#'   \code{term} and \code{value} and these rows:
#'   \code{mean_cliff_delta}, \code{median_cliff_delta}, and
#'   \code{sd_cliff_delta} (summaries of the \code{G} estimates);
#'   \code{mean_ci_width}, \code{median_ci_width}, and \code{sd_ci_width}
#'   (summaries of the widths of the replications with an interval);
#'   \code{pct_ci_less_width} (the proportion of all \code{G} replications
#'   whose interval is no wider than \code{width});
#'   \code{pct_ci_miss_low}, \code{pct_ci_miss_high}, and
#'   \code{total_type_I_error} (among the replications with an interval,
#'   the proportions whose interval lies wholly above \code{true_delta},
#'   wholly below it, and either); \code{n_1}, \code{n_2}, and
#'   \code{total_N} (the sample sizes simulated); and the inputs
#'   \code{true_delta}, \code{supposed_delta} (\code{NA} when
#'   \code{specified_N} was given), \code{ratio}, \code{width},
#'   \code{conf_level}, and \code{assurance} (present only when an
#'   assurance was supplied). The \code{simulation} attribute records the
#'   numbers of replications attempted, with an interval, and without one
#'   (see the section on the record of the simulation).
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
#'   ordinal questions. \emph{Psychological Bulletin, 114}(3), 494--509.
#'   \doi{10.1037/0033-2909.114.3.494}
#'
#' Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
#'   improved confidence interval. \emph{Journal of Modern Applied
#'   Statistical Methods, 3}(2), 322--332. \doi{10.22237/jmasm/1099267560}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_cliff_delta}}, \code{\link{cliff_delta}}
#'
#' @examples
#' # Small G keeps the Monte Carlo sweep fast; raise G for a real plan.
#' ss_aipe_cliff_delta_sensitivity(
#'   true_delta = 0.30, supposed_delta = 0.30,
#'   width = 0.30, G = 25, print_iter = FALSE, seed = 113
#' )
#'
#' # A large effect with six per group: some replications are completely
#' # separated and have no interval. The width and coverage rows describe
#' # the others, and one warning says how many were without an interval.
#' ss_aipe_cliff_delta_sensitivity(
#'   true_delta = 0.90, specified_N = 12, width = 0.40, G = 100,
#'   seed = 113
#' )
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_cliff_delta_sensitivity <- function(true_delta = NULL,
                                            supposed_delta = NULL,
                                            ratio = 1,
                                            width,
                                            specified_N = NULL,
                                            conf_level = 0.95,
                                            assurance = NULL,
                                            G = 1000, print_iter = FALSE,
                                            filename = NULL,
                                            seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  if (is.null(supposed_delta) && is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify either 'supposed_delta' or 'specified_N'.",
      arg = c("supposed_delta", "specified_N"))
  if (!is.null(supposed_delta) && !is.null(specified_N))
    .dmar_refuse_argument("You must specify 'supposed_delta' or 'specified_N', but not both.",
                          arg = c("supposed_delta", "specified_N"))
  if (is.null(true_delta) || !is.numeric(true_delta) || abs(true_delta) >= 1)
    .dmar_refuse_argument(
      "'true_delta' must be a single value in (-1, 1).", arg = "true_delta")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  if (!is.null(supposed_delta)) {
    plan <- ss_aipe_cliff_delta(true_delta = supposed_delta, width = width,
                                conf_level = conf_level, ratio = ratio,
                                assurance = assurance)
    n1 <- plan$value[plan$term == "n_1"]
    n2 <- plan$value[plan$term == "n_2"]
  } else {
    n2 <- floor(specified_N / (1 + ratio))
    n1 <- specified_N - n2
  }
  n1 <- as.integer(n1); n2 <- as.integer(n2)
  if (n1 < 2 || n2 < 2)
    .dmar_refuse_argument("Resolved per-group sample size is < 2.", arg = NULL)

  # Normal-population shift that produces the requested Cliff's delta:
  # delta = 2 * Phi(Delta / sqrt(2)) - 1.
  Delta <- sqrt(2) * stats::qnorm((1 + true_delta) / 2)

  cd_hat   <- numeric(G)
  ci_lo    <- rep(NA_real_, G)
  ci_hi    <- rep(NA_real_, G)
  ci_width <- rep(NA_real_, G)
  tI_lower <- rep(NA, G)
  tI_upper <- rep(NA, G)

  # A replication whose groups are completely separated (or whose scores
  # are all equal) has an estimated variance of 0, and cliff_delta() then
  # reports both limits equal to the estimate with a warning of class
  # dmar_interval_degenerate. Those limits are not a confidence interval,
  # so the replication has none: its estimate is kept and its limits,
  # width, and coverage indicators stay NA. The warning is muffled by its
  # class here, counted, and reported once after the loop.
  state <- new.env(parent = emptyenv())
  state$degenerate <- FALSE

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    g1 <- stats::rnorm(n1, mean = Delta, sd = 1)
    g2 <- stats::rnorm(n2, mean = 0,    sd = 1)
    state$degenerate <- FALSE
    res <- withCallingHandlers(
      cliff_delta(group_1 = g1, group_2 = g2, conf_level = conf_level),
      dmar_interval_degenerate = function(w) {
        state$degenerate <- TRUE
        invokeRestart("muffleWarning")
      })
    cd_hat[g] <- res$value[res$term == "cliff_delta"]
    if (!state$degenerate) {
      lo <- res$value[res$term == "ci_lower"]
      hi <- res$value[res$term == "ci_upper"]
      ci_lo[g]    <- lo
      ci_hi[g]    <- hi
      ci_width[g] <- hi - lo
      tI_lower[g] <- true_delta < lo
      tI_upper[g] <- true_delta > hi
    }
  }

  has_interval <- !is.na(ci_width)
  met_width    <- has_interval & ci_width <= width
  n_without_interval <- sum(!has_interval)
  if (n_without_interval > 0L) {
    .dmar_warn(sprintf(paste0(
      "%d of %d replications had completely separated groups (or no ",
      "variation), where the estimated variance of Cliff's delta is 0 ",
      "and no confidence interval exists. The width, width-target, and ",
      "coverage rows describe the other %d replications."),
      n_without_interval, G, G - n_without_interval),
      "replications_without_interval",
      n_without_interval = n_without_interval, G = G)
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(cliff_delta = cd_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  # A summary over the replications with an interval (the others hold NA)
  # is NA when there are none (mean() of no values is NaN).
  over <- function(f, v) if (any(has_interval)) f(v, na.rm = TRUE) else
    NA_real_
  out <- data.frame(
    term  = c("mean_cliff_delta", "median_cliff_delta", "sd_cliff_delta",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_1", "n_2", "total_N",
              "true_delta", "supposed_delta",
              "ratio", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(cd_hat, na.rm = TRUE),
              stats::median(cd_hat, na.rm = TRUE),
              stats::sd(cd_hat, na.rm = TRUE),
              over(mean, ci_width),
              over(stats::median, ci_width),
              over(stats::sd, ci_width),
              if (any(has_interval)) sum(met_width) / sum(has_interval) else
                NA_real_,
              over(mean, tI_lower),
              over(mean, tI_upper),
              over(mean, tI_lower | tI_upper),
              n1, n2, n1 + n2,
              true_delta,
              if (is.null(supposed_delta)) NA_real_ else supposed_delta,
              ratio, width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "simulation") <- .ss_sens_simulation(G, has_interval, met_width,
                                                 conf_level)
  out
}
