# Monte Carlo sensitivity of the crossed units by raters planner.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for the Mean of a Crossed Design of Units and Raters
#'
#' @description
#' Shows by Monte Carlo simulation what a plan from
#' \code{\link{ss_aipe_crossed}} delivers when the true variance
#' components differ from the planning values. The design is either
#' planned from the supposed components, as a researcher would plan it,
#' or fixed by the user; complete arrays of that size are then drawn from
#' the population with the true components, the interval of
#' \code{\link{ci_crossed}} is formed from each, and the realized widths
#' and the coverage of the true mean are summarized. A planning value of
#' the rater component that is too small is the costly mistake in this
#' design, since that component sets the floor below which no number of
#' units can bring the interval.
#'
#' @param true_var_unit,true_var_rater,true_var_residual The data
#'   generating variance components of the unit effects, the rater
#'   effects, and the residuals (see \code{\link{ss_aipe_crossed}}). The
#'   first two must be nonnegative and the third positive.
#' @param supposed_var_unit,supposed_var_rater,supposed_var_residual
#'   The supposed variance components, the planning values the design is
#'   planned from (not estimates computed from data), all three together.
#'   Supply these or \code{specified_n_units} and
#'   \code{specified_n_raters}, but not both.
#' @param specified_n_units,specified_n_raters A design to evaluate
#'   (whole numbers of at least 2, both together) in place of a planned
#'   one; the planner is then bypassed.
#' @param width Desired full width of the confidence interval for the
#'   mean, against which the realized widths are judged.
#' @param conf_level Confidence level of the interval (default
#'   \code{0.95}).
#' @param assurance Optional assurance, strictly between 0 and 1, passed
#'   to \code{\link{ss_aipe_crossed}} when the design is planned; echoed
#'   in the output when supplied, with a specified design as well.
#' @param n_raters,unit_cost,rater_cost,judgment_cost,quantile,max_units,max_raters
#'   Passed to \code{\link{ss_aipe_crossed}} when the design is planned,
#'   and unused with a specified design.
#' @param true_grand_mean The data generating mean (default \code{0}).
#'   The widths and the coverage do not depend on it; it sets the scale
#'   of the reported estimates.
#' @param G Number of Monte Carlo replications (default \code{1000}).
#' @param print_iter Logical. Print the replication index as the
#'   simulation runs; default \code{FALSE}.
#' @param filename Optional path of a CSV file to receive the
#'   per-replication results (the grand mean, the two limits, the width,
#'   the degrees of freedom, and the two tail misses), appended when the
#'   file already holds a table and created otherwise; the default
#'   \code{NULL} writes nothing, and a throwaway run that wants the file
#'   should point it at \code{tempfile(fileext = ".csv")}.
#' @param seed Optional integer seed. The default \code{NULL} draws from,
#'   and advances, the current random number stream; an integer makes the
#'   call reproducible, and the caller's random number state is restored
#'   when the call returns.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{mean_grand_mean}, \code{median_grand_mean}, and
#'   \code{sd_grand_mean} (summaries of the \eqn{G} estimates);
#'   \code{mean_ci_width}, \code{median_ci_width}, and
#'   \code{sd_ci_width} (of the realized widths); \code{pct_ci_less_width}
#'   (the proportion of intervals no wider than \code{width});
#'   \code{pct_ci_miss_low} (the proportion lying entirely above the true
#'   mean), \code{pct_ci_miss_high} (entirely below it), and
#'   \code{total_type_I_error} (their sum); \code{mean_df} (the mean
#'   Satterthwaite degrees of freedom); and the echoes \code{n_units} and
#'   \code{n_raters} (the design evaluated), the three \code{true_var_*}
#'   and three \code{supposed_var_*} components (the latter \code{NA}
#'   for a specified design), \code{true_grand_mean}, \code{width},
#'   \code{conf_level}, and \code{assurance} when one was supplied. The
#'   proportions are on the 0 to 1 scale.
#'
#' @details
#' Each replication draws \eqn{n_u} unit effects, \eqn{n_r} rater
#' effects, and \eqn{n_u n_r} residuals from normal distributions with
#' the true variances, forms the complete array, and computes the
#' interval exactly as \code{\link{ci_crossed}} does: the ANOVA estimates
#' of the components, negative ones truncated at zero, plugged into the
#' variance of the grand mean, with Satterthwaite (1946) degrees of
#' freedom. The summary separates two sources of error: the planning
#' values (compare a well specified run with a misspecified one) and the
#' approximations of the plan itself (a well specified run whose
#' \code{pct_ci_less_width} is well below one half, or below the assurance,
#' or whose coverage is far from \code{conf_level}; a plan without an
#' assurance meets its width somewhat more than half the time, as
#' \code{\link{ss_aipe_crossed}} explains). A plan with
#' \code{quantile = "normal"} shows the second kind, since its target
#' ignores the \emph{t} quantile the interval uses; with a rater
#' component at or near zero it can fall far short.
#'
#' @references
#' Satterthwaite, F. E. (1946). An approximate distribution of
#'   estimates of variance components. \emph{Biometrics Bulletin, 2}(6),
#'   110--114. \doi{10.2307/3002019}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_crossed}}, \code{\link{ci_crossed}}
#'
#' @examples
#' # G = 200 keeps the example quick; a reported study deserves the
#' # default G = 1000 or more. The first plan on the help page of
#' # ss_aipe_crossed(), evaluated when its planning values are right:
#' ss_aipe_crossed_sensitivity(
#'   true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
#'   supposed_var_unit = 0.50, supposed_var_rater = 0.10,
#'   supposed_var_residual = 0.40,
#'   width = 0.40, G = 200, seed = 113)
#'
#' # The same plan when the raters in fact differ twice as much as
#' # planned: nearly every interval is now wider than the target, since
#' # the design has too few raters for the larger rater component.
#' ss_aipe_crossed_sensitivity(
#'   true_var_unit = 0.50, true_var_rater = 0.20, true_var_residual = 0.40,
#'   supposed_var_unit = 0.50, supposed_var_rater = 0.10,
#'   supposed_var_residual = 0.40,
#'   width = 0.40, G = 200, seed = 113)
#'
#' @keywords design htest
#'
#' @family AIPE sample size planning
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_crossed_sensitivity <- function(true_var_unit, true_var_rater,
                                        true_var_residual,
                                        supposed_var_unit = NULL,
                                        supposed_var_rater = NULL,
                                        supposed_var_residual = NULL,
                                        specified_n_units = NULL,
                                        specified_n_raters = NULL,
                                        width, conf_level = 0.95,
                                        assurance = NULL, n_raters = NULL,
                                        unit_cost = 0, rater_cost = 0,
                                        judgment_cost = 1,
                                        quantile = c("t", "normal"),
                                        max_units = 100000,
                                        max_raters = 500,
                                        true_grand_mean = 0, G = 1000,
                                        print_iter = FALSE, filename = NULL,
                                        seed = NULL) {
  quantile <- match.arg(quantile)
  estimated <- list(supposed_var_unit, supposed_var_rater,
                    supposed_var_residual)
  specified <- list(specified_n_units, specified_n_raters)
  n_est  <- sum(!vapply(estimated, is.null, logical(1)))
  n_spec <- sum(!vapply(specified, is.null, logical(1)))
  if (n_est %in% c(1L, 2L))
    .dmar_refuse_argument(c(
      "Supply all three of 'supposed_var_unit', 'supposed_var_rater', ",
      "and 'supposed_var_residual', or none of them."),
      arg = c("supposed_var_unit", "supposed_var_rater", "supposed_var_residual"))
  if (n_spec == 1L)
    .dmar_refuse_argument(c(
      "Supply both 'specified_n_units' and 'specified_n_raters', or ",
      "neither."), arg = c("specified_n_units", "specified_n_raters"))
  if (n_est == 0L && n_spec == 0L)
    .dmar_refuse_argument(c(
      "You must specify either the 'supposed_var_*' planning values or ",
      "the design 'specified_n_units' and 'specified_n_raters'."),
      arg = c("specified_n_units", "specified_n_raters"))
  if (n_est == 3L && n_spec == 2L)
    .dmar_refuse_argument(c(
      "You must specify the 'supposed_var_*' planning values or the ",
      "design 'specified_n_units' and 'specified_n_raters', but not both."),
      arg = c("specified_n_units", "specified_n_raters"))
  .crossed_check_number(true_var_unit, "true_var_unit", lower = 0)
  .crossed_check_number(true_var_rater, "true_var_rater", lower = 0)
  .crossed_check_number(true_var_residual, "true_var_residual", lower = 0,
                        open = TRUE)
  .crossed_check_number(width, "width", lower = 0, open = TRUE)
  .crossed_check_number(conf_level, "conf_level", lower = 0, upper = 1,
                        open = TRUE)
  .crossed_check_number(true_grand_mean, "true_grand_mean")
  # Checked here as well as in the planner, so that an assurance echoed
  # beside a specified design is a probability too, and the value column
  # stays numeric.
  if (!is.null(assurance))
    .crossed_check_number(assurance, "assurance", lower = 0, upper = 1,
                          open = TRUE)
  .check_filename(filename)
  # Two replications at least: the summaries include standard deviations.
  G <- .dmar_check_count(G, "G", min = 2L)

  if (n_est == 3L) {
    plan <- ss_aipe_crossed(true_var_unit = supposed_var_unit,
                            true_var_rater = supposed_var_rater,
                            true_var_residual = supposed_var_residual,
                            width = width, conf_level = conf_level,
                            assurance = assurance, n_raters = n_raters,
                            unit_cost = unit_cost, rater_cost = rater_cost,
                            judgment_cost = judgment_cost,
                            quantile = quantile, max_units = max_units,
                            max_raters = max_raters)
    v   <- stats::setNames(plan$value, plan$term)
    n_u <- v[["necessary_n_units"]]
    n_r <- if (is.null(n_raters)) v[["necessary_n_raters"]] else
      v[["specified_n_raters"]]
  } else {
    if (!is.null(n_raters))
      .dmar_refuse_argument(c(
        "'n_raters' fixes the raters of a planned design; with a ",
        "specified design give the raters as 'specified_n_raters'."),
        arg = c("n_raters", "specified_n_raters"))
    .crossed_check_count(specified_n_units, "specified_n_units")
    .crossed_check_count(specified_n_raters, "specified_n_raters")
    n_u <- specified_n_units
    n_r <- specified_n_raters
  }

  .dmar_local_seed(seed)

  # Each replication draws the unit effects, then the rater effects, then
  # the residuals, and keeps the array's mean squares; the G intervals are
  # then formed in one vectorized pass through the same numerics
  # ci_crossed() uses.
  sd_u <- sqrt(true_var_unit)
  sd_r <- sqrt(true_var_rater)
  sd_e <- sqrt(true_var_residual)
  ms <- matrix(NA_real_, G, 4L,
               dimnames = list(NULL, c("grand_mean", "ms_unit", "ms_rater",
                                       "ms_residual")))
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    x <- true_grand_mean +
      outer(stats::rnorm(n_u, 0, sd_u), stats::rnorm(n_r, 0, sd_r), "+") +
      matrix(stats::rnorm(n_u * n_r, 0, sd_e), n_u, n_r)
    ms[g, ] <- .crossed_ms(x)[1:4]
  }
  ci <- .crossed_full_ci(ms[, "grand_mean"], ms[, "ms_unit"],
                         ms[, "ms_rater"], ms[, "ms_residual"],
                         n_u, n_r, conf_level)
  ci_width  <- ci$ci_upper - ci$ci_lower
  miss_low  <- true_grand_mean < ci$ci_lower
  miss_high <- true_grand_mean > ci$ci_upper

  if (!is.null(filename)) {
    .write_sensitivity_csv(data.frame(grand_mean = ci$grand_mean,
                                      ci_lower = ci$ci_lower,
                                      ci_upper = ci$ci_upper,
                                      ci_width = ci_width, df = ci$df,
                                      type_I_lower = miss_low,
                                      type_I_upper = miss_high),
                           filename)
  }

  est <- function(x) if (is.null(x)) NA_real_ else x
  out <- data.frame(
    term  = c("mean_grand_mean", "median_grand_mean", "sd_grand_mean",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high",
              "total_type_I_error", "mean_df",
              "n_units", "n_raters",
              "true_var_unit", "true_var_rater", "true_var_residual",
              "supposed_var_unit", "supposed_var_rater",
              "supposed_var_residual", "true_grand_mean",
              "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(ci$grand_mean), stats::median(ci$grand_mean),
              stats::sd(ci$grand_mean),
              mean(ci_width), stats::median(ci_width), stats::sd(ci_width),
              mean(ci_width <= width), mean(miss_low), mean(miss_high),
              mean(miss_low | miss_high), mean(ci$df),
              n_u, n_r,
              true_var_unit, true_var_rater, true_var_residual,
              est(supposed_var_unit), est(supposed_var_rater),
              est(supposed_var_residual), true_grand_mean,
              width, conf_level, assurance),
    stringsAsFactors = FALSE
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(ci_width), conf_level)
}
