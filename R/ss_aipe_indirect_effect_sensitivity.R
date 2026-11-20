# Sensitivity analysis for AIPE on a mediation indirect effect (a*b).
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Indirect Effect
#'
#' @description
#' Quantifies how much misspecification of the population mediation path
#' coefficients \eqn{a} and \eqn{b} distorts an AIPE-based sample size
#' plan for the indirect effect \eqn{ab}. On each replication the
#' function simulates a three-variable mediation system
#' \eqn{X \to M \to Y} of size \emph{n} with population standardized
#' path coefficients \code{true_std_reg_coef_a} and
#' \code{true_std_reg_coef_b}, fits the two
#' regressions of \emph{M} on \emph{X} and \emph{Y} on \emph{M} and
#' \emph{X}, computes the sample indirect effect \eqn{\hat a \hat b},
#' and forms the interval the plan targeted: the symmetric Wald interval
#' from the delta method standard error
#' (\code{method = "closed_form"}) or the Monte Carlo interval
#' (\code{method = "monte_carlo"}), matching
#' \code{\link{ss_aipe_indirect_effect}}.
#'
#' @param true_std_reg_coef_a Population standardized regression
#'   coefficient of path \emph{a} (from \emph{X} to \emph{M}); the data
#'   generating value.
#' @param true_std_reg_coef_b Population standardized regression
#'   coefficient of path \emph{b} (from \emph{M} to \emph{Y} after
#'   controlling for \emph{X}); the data generating value.
#' @param supposed_std_reg_coef_a,supposed_std_reg_coef_b The standardized
#'   path coefficients supposed when the study is planned (passed to
#'   \code{\link{ss_aipe_indirect_effect}} as its
#'   \code{true_std_reg_coef_a} and \code{true_std_reg_coef_b}). They are
#'   planning assumptions, not estimates from data, and may differ from the
#'   data generating values. Supply both or neither (if neither, supply
#'   \code{specified_N}).
#' @param width Desired full width of the CI on \eqn{ab}.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{supposed_std_reg_coef_a} and \code{supposed_std_reg_coef_b}).
#' @param method One of \code{"closed_form"} (default) or
#'   \code{"monte_carlo"}; the interval computed on each replication,
#'   also forwarded to the planner when the sample size is planned from
#'   \code{supposed_std_reg_coef_a} and \code{supposed_std_reg_coef_b}. A
#'   planning call with \code{method = "monte_carlo"} runs the planner's
#'   a priori Monte
#'   Carlo search at its default \code{G}, so it takes a few seconds.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param B Number of Monte Carlo draws used for the indirect-effect CI
#'   when \code{method = "monte_carlo"} (default 5000). The interval is
#'   read from the quantiles of the draws, so \code{B} follows the rule
#'   of \code{\link{dmar_bootstrap}}: at least 400 at the 95 percent
#'   level, unless \code{options(dmar.bootstrap_exploratory = TRUE)} is
#'   set.
#' @param G Number of outer simulation replications (default 1000).
#' @param print_iter Logical.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample indirect effect \eqn{\hat a \hat b} in
#'   the column \code{std_indirect_effect}, the two confidence limits, the
#'   interval width, and two indicators of whether the interval missed
#'   \code{true_std_reg_coef_a * true_std_reg_coef_b} below or above):
#'   nothing is written when \code{filename} is \code{NULL} (the default),
#'   a new file with a header row is created otherwise, an existing file
#'   at that path is appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with rows for the mean, median, and SD of
#'   the sample indirect effect \eqn{\hat a \hat b}
#'   (\code{mean_std_indirect_effect}, \code{median_std_indirect_effect},
#'   \code{sd_std_indirect_effect}) and of the CI width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of the population value \code{true_std_indirect_effect}
#'   (\code{true_std_reg_coef_a * true_std_reg_coef_b}), and the input
#'   echoes \code{total_N}, \code{true_std_reg_coef_a},
#'   \code{true_std_reg_coef_b}, \code{supposed_std_reg_coef_a},
#'   \code{supposed_std_reg_coef_b} (\code{NA} when \code{specified_N} was
#'   given), \code{width}, and \code{conf_level}. Under
#'   \code{method = "monte_carlo"} the result carries the
#'   attributes \code{B}, \code{B_used}, \code{n_dropped},
#'   \code{B_required}, and \code{bootstrap_status} of
#'   \code{\link{dmar_bootstrap}}; they record the request each replication
#'   ran (\code{B_used} equals \code{B}), and \code{bootstrap_status} is
#'   \code{"exploratory"} when \code{B} is below \code{B_required}.
#'
#' @references
#' Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
#'   models: Quantitative strategies for communicating indirect effects.
#'   \emph{Psychological Methods, 16}(2), 93--115. \doi{10.1037/a0022658}
#'
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
#'   analysis: Introducing the model-based constrained optimization
#'   procedure. \emph{Psychological Methods, 25}(4), 496--515.
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
#'   true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.3,
#'   supposed_std_reg_coef_a = 0.4, supposed_std_reg_coef_b = 0.3,
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
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_indirect_effect_sensitivity <- function(true_std_reg_coef_a = NULL,
                                                true_std_reg_coef_b = NULL,
                                                supposed_std_reg_coef_a = NULL,
                                                supposed_std_reg_coef_b = NULL,
                                                width,
                                                specified_N = NULL,
                                                method = c("closed_form",
                                                           "monte_carlo"),
                                                conf_level = 0.95,
                                                B = 5000L,
                                                G = 1000, print_iter = FALSE,
                                                filename = NULL,
                                                seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  method <- match.arg(method)
  est_supplied <- !is.null(supposed_std_reg_coef_a) &&
    !is.null(supposed_std_reg_coef_b)
  if (!est_supplied && is.null(specified_N))
    .dmar_refuse_argument(c(
      "Supply either both 'supposed_std_reg_coef_a' and ",
      "'supposed_std_reg_coef_b', or 'specified_N'."),
      arg = c("supposed_std_reg_coef_a", "supposed_std_reg_coef_b",
              "specified_N"))
  if (est_supplied && !is.null(specified_N))
    .dmar_refuse_argument(
      "Supply supposed paths or 'specified_N', but not both.",
      arg = "specified_N")
  if (is.null(true_std_reg_coef_a) || is.null(true_std_reg_coef_b) ||
      !is.numeric(true_std_reg_coef_a) || !is.numeric(true_std_reg_coef_b))
    .dmar_refuse_argument(
      "'true_std_reg_coef_a' and 'true_std_reg_coef_b' must both be numeric.",
      arg = c("true_std_reg_coef_a", "true_std_reg_coef_b"))
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  # The Monte Carlo interval of each replication is read from the quantiles
  # of B draws, so B follows the one rule for replications
  # (R/dmar_bootstrap.R).
  if (method == "monte_carlo") {
    if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
      .dmar_refuse_argument(
        "'conf_level' must be in (0, 1).", arg = "conf_level")
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }
  .dmar_local_seed(seed)

  if (est_supplied) {
    # The mark of an exploratory run is given once, by this function at
    # its end, and not also by the planner it calls.
    plan <- withCallingHandlers(
      ss_aipe_indirect_effect(true_std_reg_coef_a = supposed_std_reg_coef_a,
                              true_std_reg_coef_b = supposed_std_reg_coef_b,
                              width = width, method = method,
                              conf_level = conf_level, B = B),
      dmar_bootstrap_exploratory = function(w) invokeRestart("muffleWarning"))
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < 10) .dmar_refuse_argument(
    "Resolved sample size is < 10; too small for a stable mediation fit.",
    arg = NULL)

  true_std_indirect_effect <- true_std_reg_coef_a * true_std_reg_coef_b
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
    eM <- stats::rnorm(n, sd = sqrt(max(1 - true_std_reg_coef_a^2,
                                        .Machine$double.eps)))
    M  <- true_std_reg_coef_a * X + eM
    # Y given X and M; in the partial setup Y = b*M + c'*X + e_Y; let c' = 0.
    eY <- stats::rnorm(n, sd = sqrt(max(1 - true_std_reg_coef_b^2,
                                        .Machine$double.eps)))
    Y  <- true_std_reg_coef_b * M + eY

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
    tI_lower[g] <- true_std_indirect_effect < lo
    tI_upper[g] <- true_std_indirect_effect > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(std_indirect_effect = ab_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_std_indirect_effect", "median_std_indirect_effect",
              "sd_std_indirect_effect",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N",
              "true_std_reg_coef_a", "true_std_reg_coef_b",
              "true_std_indirect_effect",
              "supposed_std_reg_coef_a", "supposed_std_reg_coef_b",
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
              true_std_reg_coef_a, true_std_reg_coef_b, true_std_indirect_effect,
              if (is.null(supposed_std_reg_coef_a)) NA_real_ else
                supposed_std_reg_coef_a,
              if (is.null(supposed_std_reg_coef_b)) NA_real_ else
                supposed_std_reg_coef_b,
              width, conf_level),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  out <- .ss_sens_record(out, G, !is.na(ci_width), conf_level)
  # The request each replication ran is recorded, and a run below the
  # required number (options(dmar.bootstrap_exploratory = TRUE)) is marked
  # once for the call.
  if (method == "monte_carlo") {
    out <- .dmar_boot_mark_planner(
      out, B, tail = (1 - conf_level) / 2, conf_level = conf_level,
      fun = "ss_aipe_indirect_effect_sensitivity",
      quantity = "the indirect effect", replications = "Monte Carlo draws")
  }
  out
}
