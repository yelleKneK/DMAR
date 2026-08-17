# Sensitivity analysis for AIPE on a polynomial change parameter.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Polynomial Change Parameter
#'
#' @description
#' Quantifies how much misspecification of the population
#' between-subject slope variance and within-subject error variance
#' distorts an AIPE-based sample size plan for the group-by-time
#' polynomial change parameter. On each replication the function
#' simulates two independent groups of \emph{n} subjects each, measured
#' at \eqn{M = f \times D + 1} timepoints, where every subject has a
#' true linear slope drawn from
#' \eqn{N(0, \mathrm{true\_variance\_trend})} and within-subject
#' observations have residual variance \code{true_error_variance}.
#' Subject-level OLS slopes are computed in each group, the
#' between-group difference in mean slopes (the change parameter
#' \eqn{\beta_{m1}} that \code{\link{ss_aipe_pcm}} plans for) is
#' estimated, and a two-group \emph{t}-confidence interval on that
#' difference (pooled standard error, \eqn{2n - 2} degrees of freedom)
#' is recorded. The function only handles \code{trend = "linear"} in
#' the simulator; for quadratic / cubic trends, the planner's
#' closed-form solution is still available via \code{\link{ss_aipe_pcm}}.
#'
#' @param true_variance_trend Population between-subject variance of
#'   the polynomial change coefficient (the data generating
#'   \eqn{\sigma^2_{\upsilon_m}} of Kelley & Rausch, 2011).
#' @param true_error_variance Population within-subject error variance
#'   (\eqn{\sigma^2_\epsilon}).
#' @param estimated_variance_trend Planning value of
#'   \code{variance_trend} passed to \code{\link{ss_aipe_pcm}}; supply
#'   this and \code{estimated_error_variance}, or supply
#'   \code{n_per_group}.
#' @param estimated_error_variance Planning value of
#'   \code{error_variance} passed to \code{\link{ss_aipe_pcm}}.
#' @param duration Study duration (in time units).
#' @param frequency Number of measurements per unit time. Total
#'   timepoints = \eqn{f \times D + 1}.
#' @param width Desired full width of the CI on the between-group
#'   difference in change parameters (\eqn{\beta_{m1}}).
#' @param n_per_group Number of subjects to evaluate (incompatible
#'   with the estimated-variance arguments).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to the
#'   planner.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized estimated slope difference and CI width, the
#'   proportion of intervals at or below \code{width}, tail-specific and
#'   overall non-coverage of the population slope difference (0 by
#'   construction in this simulator), and the input echoes, including \code{assurance} (present only when an
#'   assurance was supplied).
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2011). Sample size planning for
#'   longitudinal models: Accuracy in parameter estimation for
#'   polynomial change parameters. \emph{Psychological Methods, 16}(4),
#'   391--405. \doi{10.1037/a0023352}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapters 11 and 15.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_pcm}}, \code{\link{ss_aipe_mixed_effects_sensitivity}}
#'
#' @examples
#' # Every replication simulates two full groups of subjects, fits a
#' # slope for each subject, and forms a confidence interval on the
#' # difference in mean slopes, so the sweep is not run at example time.
#' # The G below is far smaller than a reported sensitivity study would
#' # use; the default of 1000 is the realistic setting. The call is:
#' # set.seed(113)
#' # ss_aipe_pcm_sensitivity(
#' #   true_variance_trend       = 0.003,
#' #   true_error_variance       = 0.0262,
#' #   estimated_variance_trend  = 0.003,
#' #   estimated_error_variance  = 0.0262,
#' #   duration  = 4, frequency = 1,
#' #   width     = 0.05,
#' #   G = 20, print_iter = FALSE
#' # )
#'
#' @keywords design multivariate
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export
ss_aipe_pcm_sensitivity <- function(true_variance_trend = NULL,
                                    true_error_variance = NULL,
                                    estimated_variance_trend = NULL,
                                    estimated_error_variance = NULL,
                                    duration, frequency, width,
                                    n_per_group = NULL,
                                    conf_level = 0.95,
                                    assurance = NULL,
                                    G = 1000, print_iter = FALSE,
                                    save = FALSE,
                                    filename = "ss_aipe_pcm_sensitivity_result.csv") {
  est_supplied <- !is.null(estimated_variance_trend) && !is.null(estimated_error_variance)
  if (!est_supplied && is.null(n_per_group))
    stop("Supply (estimated_variance_trend, estimated_error_variance) or 'n_per_group'.", call. = FALSE)
  if (est_supplied && !is.null(n_per_group))
    stop("Supply estimated-variance args or 'n_per_group', not both.", call. = FALSE)
  if (is.null(true_variance_trend) || true_variance_trend <= 0)
    stop("'true_variance_trend' must be a positive number.", call. = FALSE)
  if (is.null(true_error_variance) || true_error_variance <= 0)
    stop("'true_error_variance' must be a positive number.", call. = FALSE)

  if (est_supplied) {
    plan <- ss_aipe_pcm(variance_trend = estimated_variance_trend,
                        error_variance = estimated_error_variance,
                        duration = duration, frequency = frequency,
                        width = width, conf_level = conf_level,
                        trend = "linear", assurance = assurance)
    n <- plan$value[plan$term == "necessary_n_per_group"]
  } else {
    n <- n_per_group
  }
  n <- as.integer(n)
  if (n < 4) stop("Resolved sample size is < 4.", call. = FALSE)

  # Time grid for each subject: M = f*D + 1 equally spaced points on [0, D].
  M <- as.integer(frequency * duration + 1L)
  if (M < 3) stop("Resolved number of timepoints (frequency*duration + 1) is < 3.", call. = FALSE)
  t_grid <- seq(0, duration, length.out = M)

  crit <- stats::qt(1 - (1 - conf_level) / 2, df = 2 * n - 2)

  beta_hat <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_w     <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  sd_pi <- sqrt(true_variance_trend)
  sd_e  <- sqrt(true_error_variance)

  # One group's subject-level OLS slopes. Each subject's true slope is
  # drawn from N(0, true_variance_trend); its observed slope additionally
  # carries the OLS estimation error implied by true_error_variance on the
  # actual time grid, so the observed subject slope has variance
  # (true_variance_trend + V), with V the per-unit-time slope error
  # variance that ss_aipe_pcm() targets.
  group_slopes <- function() {
    pi_i <- stats::rnorm(n, mean = 0, sd = sd_pi)
    vapply(seq_len(n), function(i) {
      y <- pi_i[i] * t_grid + stats::rnorm(M, mean = 0, sd = sd_e)
      unname(stats::coef(stats::lm(y ~ t_grid))[2L])
    }, numeric(1))
  }

  # ss_aipe_pcm() plans the width of a CI on the BETWEEN-group difference
  # in mean slopes (the change parameter beta_m1), so the simulator must
  # form that two-group difference: two independent groups of n subjects,
  # a pooled-variance standard error, and 2n - 2 degrees of freedom. The
  # population difference is 0 here (both groups share the null mean
  # slope), so coverage is assessed against 0.
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    slopes_c <- group_slopes()
    slopes_t <- group_slopes()
    b <- mean(slopes_t) - mean(slopes_c)
    s2_pooled <- ((n - 1) * stats::var(slopes_c) +
                  (n - 1) * stats::var(slopes_t)) / (2 * n - 2)
    se <- sqrt(s2_pooled * (1 / n + 1 / n))
    lo <- b - crit * se
    hi <- b + crit * se
    beta_hat[g] <- b
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_w[g]     <- hi - lo
    tI_lower[g] <- 0 < lo
    tI_upper[g] <- 0 > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(slope_diff = beta_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_w,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_slope_diff", "median_slope_diff", "sd_slope_diff",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_per_group", "n_timepoints",
              "true_variance_trend", "true_error_variance",
              "estimated_variance_trend", "estimated_error_variance",
              "duration", "frequency", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(beta_hat, na.rm = TRUE),
              stats::median(beta_hat, na.rm = TRUE),
              stats::sd(beta_hat, na.rm = TRUE),
              mean(ci_w, na.rm = TRUE),
              stats::median(ci_w, na.rm = TRUE),
              stats::sd(ci_w, na.rm = TRUE),
              mean(ci_w <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, M,
              true_variance_trend, true_error_variance,
              if (is.null(estimated_variance_trend)) NA_real_ else estimated_variance_trend,
              if (is.null(estimated_error_variance)) NA_real_ else estimated_error_variance,
              duration, frequency, width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
