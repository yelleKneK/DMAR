# Sensitivity analysis for AIPE on a mediation indirect effect (a*b).
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Indirect Effect
#'
#' @description
#' Quantifies how much misspecification of the population mediation path
#' coefficients \eqn{a} and \eqn{b} distorts an AIPE-based sample size
#' plan for the indirect effect \eqn{ab}. On each replication the
#' function simulates a three-variable mediation system
#' \eqn{X \to M \to Y} of size \emph{n} with population path
#' coefficients \code{true_a} and \code{true_b}, fits the two
#' regressions of \emph{M} on \emph{X} and \emph{Y} on \emph{M} and
#' \emph{X}, computes the sample indirect effect \eqn{\hat a \hat b},
#' and forms the interval the plan targeted: the symmetric Wald interval
#' from the delta method standard error
#' (\code{method = "closed_form"}) or the Monte Carlo interval
#' (\code{method = "monte_carlo"}), matching
#' \code{\link{ss_aipe_indirect_effect}}.
#'
#' @param true_a Population path coefficient \emph{a} (from \emph{X} to
#'   \emph{M}); the data generating value.
#' @param true_b Population path coefficient \emph{b} (from \emph{M} to
#'   \emph{Y} after controlling for \emph{X}); the data generating
#'   value.
#' @param estimated_a,estimated_b Path coefficients used to plan the
#'   study (passed to \code{\link{ss_aipe_indirect_effect}}). Supply
#'   both or neither (if neither, supply \code{specified_N}).
#' @param width Desired full width of the CI on \eqn{ab}.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{estimated_a} / \code{estimated_b}).
#' @param method One of \code{"closed_form"} (default) or
#'   \code{"monte_carlo"}; the interval computed on each replication,
#'   also forwarded to the planner when the sample size is planned from
#'   \code{estimated_a} and \code{estimated_b}. A planning call with
#'   \code{method = "monte_carlo"} runs the planner's a priori Monte
#'   Carlo search at its default \code{G}, so it takes a few seconds.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param B Number of Monte Carlo draws used for the indirect-effect CI
#'   when \code{method = "monte_carlo"} (default 5000).
#' @param G Number of outer simulation replications (default 1000).
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   \eqn{\hat a \hat b} and the CI width, the proportion of intervals
#'   at or below \code{width}, tail-specific and overall non-coverage
#'   of the population value \code{true_a * true_b}, and the input
#'   echoes.
#'
#' @references
#' Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
#'   models: Quantitative strategies for communicating indirect effects.
#'   \emph{Psychological Methods, 16}(2), 93--115. \doi{10.1037/a0022658}
#'
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
#'   analysis: Introducing the model-based constrained optimization
#'   procedure. \emph{Psychological Methods, 25}, 496--515.
#'   \doi{10.1037/met0000259}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_indirect_effect}}, \code{\link{var_indirect_effect}}
#'
#' @examples
#' # Reduced replications and a wide target width keep this fast.
#' set.seed(113)
#' ss_aipe_indirect_effect_sensitivity(
#'   true_a = 0.4, true_b = 0.3,
#'   estimated_a = 0.4, estimated_b = 0.3,
#'   width = 0.40, method = "closed_form",
#'   G = 50, print_iter = FALSE
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
ss_aipe_indirect_effect_sensitivity <- function(true_a = NULL, true_b = NULL,
                                                estimated_a = NULL, estimated_b = NULL,
                                                width,
                                                specified_N = NULL,
                                                method = c("closed_form",
                                                           "monte_carlo"),
                                                conf_level = 0.95,
                                                B = 5000L,
                                                G = 1000, print_iter = FALSE,
                                                save = FALSE,
                                                filename = "ss_aipe_indirect_effect_sensitivity_result.csv") {
  method <- match.arg(method)
  est_supplied <- !is.null(estimated_a) && !is.null(estimated_b)
  if (!est_supplied && is.null(specified_N))
    stop("Supply either both 'estimated_a' and 'estimated_b', or 'specified_N'.", call. = FALSE)
  if (est_supplied && !is.null(specified_N))
    stop("Supply estimated paths or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_a) || is.null(true_b) || !is.numeric(true_a) || !is.numeric(true_b))
    stop("'true_a' and 'true_b' must both be numeric.", call. = FALSE)

  if (est_supplied) {
    plan <- ss_aipe_indirect_effect(a = estimated_a, b = estimated_b,
                                    width = width, method = method,
                                    conf_level = conf_level, B = B)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < 10) stop("Resolved sample size is < 10; too small for a stable mediation fit.")

  true_ab <- true_a * true_b
  crit    <- stats::qnorm(1 - (1 - conf_level) / 2)

  ab_hat   <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_width <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    # Standardized variables; residual SDs chosen so each equation has
    # unit total variance.
    X  <- stats::rnorm(n)
    eM <- stats::rnorm(n, sd = sqrt(max(1 - true_a^2, .Machine$double.eps)))
    M  <- true_a * X + eM
    # Y given X and M; in the partial setup Y = b*M + c'*X + e_Y; let c' = 0.
    eY <- stats::rnorm(n, sd = sqrt(max(1 - true_b^2, .Machine$double.eps)))
    Y  <- true_b * M + eY

    fit_a <- stats::lm(M ~ X)
    fit_b <- stats::lm(Y ~ M + X)
    a_hat <- stats::coef(fit_a)["X"]
    b_hat <- stats::coef(fit_b)["M"]
    se_a  <- sqrt(stats::vcov(fit_a)["X", "X"])
    se_b  <- sqrt(stats::vcov(fit_b)["M", "M"])

    ab <- as.numeric(a_hat * b_hat)
    if (method == "closed_form") {
      se_ab <- sqrt(a_hat^2 * se_b^2 + b_hat^2 * se_a^2)
      lo <- ab - crit * se_ab
      hi <- ab + crit * se_ab
    } else {
      a_draws <- stats::rnorm(B, mean = a_hat, sd = se_a)
      b_draws <- stats::rnorm(B, mean = b_hat, sd = se_b)
      ab_draws <- a_draws * b_draws
      qs <- stats::quantile(ab_draws, probs = c((1 - conf_level) / 2,
                                                1 - (1 - conf_level) / 2),
                            names = FALSE, na.rm = TRUE)
      lo <- qs[1]; hi <- qs[2]
    }
    ab_hat[g]   <- ab
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_ab < lo
    tI_upper[g] <- true_ab > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(ab = ab_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_ab", "median_ab", "sd_ab",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N",
              "true_a", "true_b", "true_ab",
              "estimated_a", "estimated_b",
              "width", "conf_level"),
    value = c(mean(ab_hat, na.rm = TRUE),
              stats::median(ab_hat, na.rm = TRUE),
              stats::sd(ab_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n,
              true_a, true_b, true_ab,
              if (is.null(estimated_a)) NA_real_ else estimated_a,
              if (is.null(estimated_b)) NA_real_ else estimated_b,
              width, conf_level),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
