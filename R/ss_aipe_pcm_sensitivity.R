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
#' true linear slope drawn from a normal distribution with mean 0 and
#' variance \code{true_var_trend}, and within-subject
#' observations have residual variance \code{true_var_error}.
#' Subject-level OLS slopes are computed in each group, the
#' between-group difference in mean slopes (the change parameter
#' \eqn{\beta_{m1}} that \code{\link{ss_aipe_pcm}} plans for) is
#' estimated, and a two-group \emph{t}-confidence interval on that
#' difference (pooled standard error, \eqn{2n - 2} degrees of freedom)
#' is recorded. The function only handles \code{trend = "linear"} in
#' the simulator; for quadratic / cubic trends, the planner's
#' closed-form solution is still available via \code{\link{ss_aipe_pcm}}.
#'
#' @param true_var_trend Population between-subject variance of
#'   the polynomial change coefficient (the data generating
#'   \eqn{\sigma^2_{\upsilon_m}} of Kelley & Rausch, 2011).
#' @param true_var_error Population within-subject error variance
#'   (\eqn{\sigma^2_\epsilon}).
#' @param supposed_var_trend Planning value of
#'   \code{true_var_trend} passed to \code{\link{ss_aipe_pcm}}; supply
#'   this and \code{supposed_var_error}, or supply
#'   \code{specified_n_per_group}.
#' @param supposed_var_error Planning value of
#'   \code{true_var_error} passed to \code{\link{ss_aipe_pcm}}.
#' @param duration Study duration (in time units).
#' @param frequency Number of measurements per unit time. The number of
#'   measurement occasions is \eqn{f \times D + 1} (the row
#'   \code{n_occasions}).
#' @param width Desired full width of the CI on the between-group
#'   difference in change parameters (\eqn{\beta_{m1}}).
#' @param specified_n_per_group Number of subjects per group to evaluate
#'   (incompatible with the supposed-variance arguments).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to the
#'   planner.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param filename Optional path to a CSV file; when supplied, the
#'   per-replication results (the slope difference, its interval limits and
#'   width, and the two tail misses) are written there, appended when the
#'   file already exists, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}; the default \code{NULL} writes
#'   nothing.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
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
#' # Every replication simulates two full groups of subjects, fits a slope
#' # for each subject, and forms a confidence interval on the difference in
#' # mean slopes. G = 10 keeps the example quick; a reported sensitivity
#' # study deserves the default G = 1000. With the planning values equal to
#' # the population values, the realized mean width should sit close to
#' # the target.
#' set.seed(113)
#' ss_aipe_pcm_sensitivity(
#'   true_var_trend = 0.003, true_var_error = 0.0262,
#'   supposed_var_trend = 0.003, supposed_var_error = 0.0262,
#'   duration = 4, frequency = 1, width = 0.05,
#'   G = 10, print_iter = FALSE
#' )
#'
#' @keywords design multivariate
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_pcm_sensitivity <- function(true_var_trend = NULL,
                                    true_var_error = NULL,
                                    supposed_var_trend = NULL,
                                    supposed_var_error = NULL,
                                    duration, frequency, width,
                                    specified_n_per_group = NULL,
                                    conf_level = 0.95,
                                    assurance = NULL,
                                    G = 1000, print_iter = FALSE,
                                    filename = NULL,
                                    seed = NULL) {
  est_supplied <- !is.null(supposed_var_trend) && !is.null(supposed_var_error)
  if (!est_supplied && is.null(specified_n_per_group))
    .dmar_refuse_argument(
      "Supply (supposed_var_trend, supposed_var_error) or 'specified_n_per_group'.",
      arg = "specified_n_per_group")
  if (est_supplied && !is.null(specified_n_per_group))
    .dmar_refuse_argument(
      "Supply the supposed variances ('supposed_var_trend', 'supposed_var_error') or 'specified_n_per_group', not both.",
      arg = "specified_n_per_group")
  if (is.null(true_var_trend) || true_var_trend <= 0)
    .dmar_refuse_argument(
      "'true_var_trend' must be a positive number.",
      arg = "true_var_trend")
  if (is.null(true_var_error) || true_var_error <= 0)
    .dmar_refuse_argument(
      "'true_var_error' must be a positive number.",
      arg = "true_var_error")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  if (est_supplied) {
    plan <- ss_aipe_pcm(true_var_trend = supposed_var_trend,
                        true_var_error = supposed_var_error,
                        duration = duration, frequency = frequency,
                        width = width, conf_level = conf_level,
                        trend = "linear", assurance = assurance)
    n <- plan$value[plan$term == "necessary_n_per_group"]
  } else {
    n <- specified_n_per_group
  }
  n <- as.integer(n)
  if (n < 4) .dmar_refuse_argument("Resolved sample size is < 4.", arg = NULL)

  # Time grid for each subject: M = f*D + 1 equally spaced points on [0, D].
  M <- as.integer(frequency * duration + 1L)
  if (M < 3) .dmar_refuse_argument(
    "Resolved number of timepoints (frequency*duration + 1) is < 3.",
    arg = NULL)
  t_grid <- seq(0, duration, length.out = M)

  crit <- stats::qt(1 - (1 - conf_level) / 2, df = 2 * n - 2)

  beta_hat <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_w     <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  sd_pi <- sqrt(true_var_trend)
  sd_e  <- sqrt(true_var_error)

  # One group's subject-level OLS slopes. Each subject's true slope is
  # drawn from N(0, true_var_trend); its observed slope additionally
  # carries the OLS estimation error implied by true_var_error on the
  # actual time grid, so the observed subject slope has variance
  # (true_var_trend + V), with V the per-unit-time slope error
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

  if (!is.null(filename)) {
    per_rep <- data.frame(slope_diff = beta_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_w,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_slope_diff", "median_slope_diff", "sd_slope_diff",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_per_group", "n_occasions",
              "true_var_trend", "true_var_error",
              "supposed_var_trend", "supposed_var_error",
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
              true_var_trend, true_var_error,
              if (is.null(supposed_var_trend)) NA_real_ else supposed_var_trend,
              if (is.null(supposed_var_error)) NA_real_ else supposed_var_error,
              duration, frequency, width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(ci_w), conf_level)
}
