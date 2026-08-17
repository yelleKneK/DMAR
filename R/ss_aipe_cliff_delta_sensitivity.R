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
#' \strong{Data generating mechanism.} The simulator draws each sample
#' from a normal distribution and chooses the mean shift so that the
#' implied Cliff's delta equals \code{true_delta}. For normal samples
#' \eqn{\delta = 2 \Phi(\Delta/\sqrt{2}) - 1} where \eqn{\Delta} is the
#' standardized mean difference, so the simulator sets
#' \eqn{\Delta = \sqrt{2} \cdot \Phi^{-1}((1 + \delta)/2)}.
#'
#' @param true_delta Population Cliff's delta (the data generating
#'   value); in \eqn{(-1, 1)}.
#' @param estimated_delta Planning value passed to
#'   \code{\link{ss_aipe_cliff_delta}}; supply this or \code{specified_N}
#'   but not both.
#' @param ratio Allocation ratio \eqn{n_1 / n_2} (default 1).
#' @param width Desired full width of the CI on Cliff's delta.
#' @param specified_N Total sample size to evaluate (split per
#'   \code{ratio}).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized Cliff's delta and CI width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of \code{true_delta}, and the input echoes, including \code{assurance} (present only when an
#'   assurance was supplied).
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
#'   ordinal questions. \emph{Psychological Bulletin, 114}(3), 494--509.
#'   \doi{10.1037/0033-2909.114.3.494}
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
#' set.seed(113)
#' # Small G keeps the Monte Carlo sweep fast; raise G for a real plan.
#' ss_aipe_cliff_delta_sensitivity(
#'   true_delta = 0.30, estimated_delta = 0.30,
#'   width = 0.30, G = 25, print_iter = FALSE
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
#' @export
ss_aipe_cliff_delta_sensitivity <- function(true_delta = NULL,
                                            estimated_delta = NULL,
                                            ratio = 1,
                                            width,
                                            specified_N = NULL,
                                            conf_level = 0.95,
                                            assurance = NULL,
                                            G = 1000, print_iter = FALSE,
                                            save = FALSE,
                                            filename = "ss_aipe_cliff_delta_sensitivity_result.csv") {
  if (is.null(estimated_delta) && is.null(specified_N))
    stop("You must specify either 'estimated_delta' or 'specified_N'.", call. = FALSE)
  if (!is.null(estimated_delta) && !is.null(specified_N))
    stop("You must specify 'estimated_delta' or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_delta) || !is.numeric(true_delta) || abs(true_delta) >= 1)
    stop("'true_delta' must be a single value in (-1, 1).", call. = FALSE)

  if (!is.null(estimated_delta)) {
    plan <- ss_aipe_cliff_delta(delta = estimated_delta, width = width,
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
    stop("Resolved per-group sample size is < 2.", call. = FALSE)

  # Normal-population shift that produces the requested Cliff's delta:
  # delta = 2 * Phi(Delta / sqrt(2)) - 1.
  Delta <- sqrt(2) * stats::qnorm((1 + true_delta) / 2)

  cd_hat   <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_width <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    g1 <- stats::rnorm(n1, mean = Delta, sd = 1)
    g2 <- stats::rnorm(n2, mean = 0,    sd = 1)
    res <- cliff_delta(group_1 = g1, group_2 = g2, conf_level = conf_level)
    cd <- res$value[res$term == "cliff_delta"]
    lo <- res$value[res$term == "lower_limit"]
    hi <- res$value[res$term == "upper_limit"]
    cd_hat[g]   <- cd
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_delta < lo
    tI_upper[g] <- true_delta > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(cliff_delta = cd_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_cliff_delta", "median_cliff_delta", "sd_cliff_delta",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_1", "n_2", "total_N",
              "true_delta", "estimated_delta",
              "ratio", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(cd_hat, na.rm = TRUE),
              stats::median(cd_hat, na.rm = TRUE),
              stats::sd(cd_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n1, n2, n1 + n2,
              true_delta,
              if (is.null(estimated_delta)) NA_real_ else estimated_delta,
              ratio, width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
