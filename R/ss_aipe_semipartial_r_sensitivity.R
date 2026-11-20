# Sensitivity analysis for AIPE on a semipartial correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Semipartial Correlation
#'
#' @description
#' Quantifies how much misspecification of the population semipartial
#' correlation, or of the full-model squared multiple correlation,
#' distorts an AIPE-based sample size plan. The function constructs a
#' \eqn{(p + 1) \times (p + 1)} population covariance matrix for \emph{Y},
#' \emph{X}, and \eqn{Z_1, \ldots, Z_{p-1}} whose semipartial correlation
#' between \emph{Y} and \emph{X} (partialing \eqn{Z_1, \ldots, Z_{p-1}} out of
#' \emph{X} only, not out of \emph{Y}) equals \code{true_r_sp} and whose
#' full-model \eqn{R^2} equals \code{true_R2_full}. On each replication it
#' draws an \emph{N}-row sample, computes the sample semipartial
#' correlation and the sample full-model \eqn{R^2}, and forms the Wald
#' interval that \code{\link{ss_aipe_semipartial_r}} plans for: the
#' estimate plus and minus the normal critical value times the square root
#' of the first-order variance (\code{\link{var_semipartial_r}} with
#' \code{under_null = FALSE}) at the two sample values.
#'
#' @param true_r_sp Population semipartial correlation \eqn{\rho_{sp}}; must lie in
#'   \eqn{(-1, 1)}.
#' @param supposed_r_sp Planning value passed to
#'   \code{\link{ss_aipe_semipartial_r}}; supply this or
#'   \code{specified_N} but not both.
#' @param p Number of predictors in the regression of \emph{Y} on \emph{X}
#'   and the variables partialed out, so \code{p - 1} variables are
#'   partialed out and \code{p} is at least 2. A categorical variable with
#'   \emph{a} levels counts as \emph{a} - 1 predictors. It is the \code{p}
#'   of \code{\link{ss_power_R2}}, \code{\link{ci_R2}}, and
#'   \code{\link{ci_reg_coef}}: the correlation concerns one predictor of a
#'   regression with \code{p} predictors. In the notation of a model
#'   comparison, \eqn{q} is the number of predictors removed from the full
#'   model to give the restricted model (\eqn{q = p_{full} - p_{restricted}},
#'   the predictors being tested), and the number of variables partialed out
#'   is \eqn{p - q}; this correlation concerns one predictor, so \eqn{q = 1}
#'   and \eqn{p - 1} variables are partialed out. The simulated data have
#'   \eqn{p + 1} columns.
#' @param width Desired full width of the CI on the semipartial
#'   correlation.
#' @param specified_N Sample size to evaluate; a whole number of at least
#'   \eqn{p + 2}.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to
#'   \code{\link{ss_aipe_semipartial_r}}; used only with
#'   \code{supposed_r_sp}.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print the iteration index of each
#'   replication.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample semipartial correlation, the two confidence
#'   limits, the interval width, and two indicators of whether the
#'   interval missed \code{true_r_sp} below or above): nothing is
#'   written when \code{filename} is \code{NULL} (the default), a new
#'   file with a header row is created otherwise, an existing file at
#'   that path is appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param true_R2_full Population squared multiple correlation of the full
#'   model (\emph{Y} on \emph{X} and \eqn{Z_1, \ldots, Z_{p-1}}); at least
#'   \code{true_r_sp} squared and below 1. When \code{NULL} (the default)
#'   it is \code{true_r_sp} squared, a population in which the variables
#'   partialed out explain none of the outcome.
#' @param supposed_R2_full Optional planning value of the full-model
#'   \eqn{R^2}, passed to \code{\link{ss_aipe_semipartial_r}} as
#'   \code{true_R2_full}; when \code{NULL} that function plans for the largest
#'   variance over every possible value. It cannot be given with
#'   \code{specified_N}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the run reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns. With
#'   \code{assurance}, the seed also makes the simulated plan of
#'   \code{\link{ss_aipe_semipartial_r}} reproducible.
#'
#' @details
#' Under multivariate normality the sampling distribution of the sample
#' semipartial correlation depends on the population only through the
#' semipartial correlation, the full-model \eqn{R^2}, \eqn{p}, and the
#' sample size, because the statistic is unchanged by any nonsingular
#' linear transformation of \eqn{Z_1, \ldots, Z_{p-1}}, by adding any linear
#' function of them to \emph{X}, and by positive rescaling of \emph{X} and
#' \emph{Y}. The population simulated here (uncorrelated standardized
#' \emph{X} and \eqn{Z_1, \ldots, Z_{p-1}}, with \emph{Y} related to \emph{X}
#' and to \eqn{Z_1}) therefore stands for every normal population with the
#' same two values.
#'
#' The confidence limits are not truncated at \eqn{-1} and \eqn{1}, so
#' that the widths are those the planner plans; with a large
#' \eqn{|\rho_{sp}|} and a small sample a limit can pass 1 in absolute
#' value.
#'
#' The function simulates \eqn{p - 1} variables partialed out, the count
#' \code{\link{ss_aipe_semipartial_r}} plans for.
#'
#' With \code{assurance}, the sample size comes from the simulated
#' assurance step of \code{\link{ss_aipe_semipartial_r}} (at its default
#' number of draws), so two unseeded calls can simulate sample sizes that
#' differ by a case or two.
#'
#' A value of \code{true_R2_full} or \code{supposed_R2_full} that falls
#' below the squared semipartial correlation by no more than rounding to
#' two decimals allows is taken to mean that the variables partialed out
#' explain none of the outcome, with a warning of class
#' \code{dmar_semipartial_r_rounded}; a smaller value is refused.
#'
#' @return A \code{data.frame} with rows for the realized
#'   semipartial correlation, the interval width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of \code{true_r_sp}, and the input echoes, including
#'   \code{true_R2_full} (the value simulated), \code{supposed_R2_full}
#'   (\code{NA} when none was supplied), and \code{assurance} (present only
#'   when an assurance was supplied).
#'
#' @references
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_semipartial_r}}, \code{\link{ss_aipe_partial_r_sensitivity}}
#'
#' @examples
#' # G = 50 keeps the example quick; a reported check deserves the
#' # default G = 1000 or more.
#' ss_aipe_semipartial_r_sensitivity(
#'   true_r_sp = 0.30, supposed_r_sp = 0.30, p = 4, width = 0.20,
#'   G = 50, print_iter = FALSE,
#'   true_R2_full = 0.50, supposed_R2_full = 0.50, seed = 113
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
ss_aipe_semipartial_r_sensitivity <- function(true_r_sp = NULL,
                                              supposed_r_sp = NULL,
                                              p, width,
                                              specified_N = NULL,
                                              conf_level = 0.95,
                                              assurance = NULL,
                                              G = 1000, print_iter = FALSE,
                                              filename = NULL,
                                              true_R2_full = NULL,
                                              supposed_R2_full = NULL,
                                              seed = NULL) {
  .dmar_check_numbers("ss_aipe_semipartial_r_sensitivity", "supposed_r_sp")
  if (missing(p))
    .dmar_refuse_argument(
      "'p' (the number of predictors) is required.", arg = "p")
  p <- .check_partial_p(p)
  # The arithmetic below is written in J, the number of variables partialed
  # out: N - J - 1 = N - p, the smallest N is J + 3 = p + 2, and the
  # population has J + 2 = p + 1 variables.
  J <- p - 1
  if (is.null(supposed_r_sp) && is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify either 'supposed_r_sp' or 'specified_N'.",
      arg = c("supposed_r_sp", "specified_N"))
  if (!is.null(supposed_r_sp) && !is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify 'supposed_r_sp' or 'specified_N', but not both.",
      arg = c("supposed_r_sp", "specified_N"))
  if (is.null(true_r_sp) || !is.numeric(true_r_sp) ||
      length(true_r_sp) != 1L || is.na(true_r_sp) || abs(true_r_sp) >= 1)
    .dmar_refuse_argument(
      "'true_r_sp' must be a single value in (-1, 1).", arg = "true_r_sp")
  G <- .dmar_check_count(G, "G")
  if (!is.null(supposed_R2_full) && !is.null(specified_N))
    .dmar_refuse_argument(c(
      "'supposed_R2_full' is a planning value and cannot be given with ",
      "'specified_N'."), arg = c("supposed_R2_full", "specified_N"))
  if (is.null(true_R2_full)) {
    true_R2_full <- true_r_sp^2
  } else {
    .check_semipartial_R2_full(true_r_sp, true_R2_full,
                               arg = "true_R2_full", r_arg = "true_r_sp")
    if (true_R2_full >= 1)
      .dmar_refuse_argument(
        "'true_R2_full' must be below 1.", arg = "true_R2_full")
  }
  .check_filename(filename)

  # The seed is set before the plan, whose assurance step draws.
  .dmar_local_seed(seed)

  if (!is.null(supposed_r_sp)) {
    # The planner checks its true_R2_full itself; its messages name
    # `true_R2_full`, so the pair is checked here first under the names of
    # this function.
    # A pair that passes with the rounding warning would warn again inside
    # the planner, so the planner is given the value it will use.
    planning_R2_full <- supposed_R2_full
    if (!is.null(supposed_R2_full)) {
      .check_semipartial_R2_full(supposed_r_sp, supposed_R2_full,
                                 arg = "supposed_R2_full",
                                 r_arg = "supposed_r_sp")
      planning_R2_full <- max(supposed_R2_full, supposed_r_sp^2)
    }
    plan <- ss_aipe_semipartial_r(true_r_sp = supposed_r_sp, p = p,
                                  width = width, conf_level = conf_level,
                                  assurance = assurance,
                                  true_R2_full = planning_R2_full)
    n <- as.integer(plan$value[plan$term == "necessary_N"])
    if (n < J + 3)
      .dmar_refuse_argument(c(
        "The planned sample size (", n, ") is below p + 2, the smallest ",
        "this function simulates."), arg = NULL)
  } else {
    n <- .check_whole_n(specified_N, "specified_N", J + 3L)
  }

  # Population of (Y, X, Z_1, ..., Z_J), all standardized, with X and the
  # Z's mutually uncorrelated: Cor(Y, X) = true_r_sp is then the
  # semipartial correlation, and Cor(Y, Z_1) = sqrt(q) gives Y a squared
  # multiple correlation of q with the Z's, so the full-model R2 is
  # true_r_sp^2 + q. Under normality the distribution of the sample
  # semipartial correlation depends on the population through
  # (true_r_sp, true_R2_full, J, n) only, so this form loses no generality.
  q <- max(true_R2_full - true_r_sp^2, 0)
  Sigma <- diag(J + 2)
  Sigma[1, 2] <- Sigma[2, 1] <- true_r_sp
  if (q > 0) Sigma[1, 3] <- Sigma[3, 1] <- sqrt(q)

  if (!requireNamespace("MASS", quietly = TRUE))
    .dmar_refuse(
      "Package 'MASS' is required for the multivariate normal sampling step.",
      "package_required", package = "MASS")

  sp_hat   <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_width <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    dat <- MASS::mvrnorm(n = n, mu = rep(0, J + 2), Sigma = Sigma)
    Y <- dat[, 1]; X <- dat[, 2]
    Z <- dat[, 3:(J + 2), drop = FALSE]
    # Semipartial: corr(Y, X residualized on the Z's). The sample
    # full-model R2 is its square plus the R2 of Y on the Z's alone.
    e_x  <- stats::residuals(stats::lm(X ~ Z))
    r_sp <- stats::cor(Y, e_x)
    R2_full_hat <- r_sp^2 + summary(stats::lm(Y ~ Z))$r.squared
    # The Wald interval the planner plans for, with the first-order
    # variance at the two sample values; the limits are not truncated.
    se  <- sqrt(.semipartial_r_var_numerator(r_sp, R2_full_hat) / (n - J - 1))
    lo  <- r_sp - crit * se
    hi  <- r_sp + crit * se
    sp_hat[g]   <- r_sp
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_r_sp < lo
    tI_upper[g] <- true_r_sp > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(r_sp = sp_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_r_sp", "median_r_sp", "sd_r_sp",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "p",
              "true_r_sp", "supposed_r_sp", "width", "conf_level",
              "true_R2_full", "supposed_R2_full",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(sp_hat, na.rm = TRUE),
              stats::median(sp_hat, na.rm = TRUE),
              stats::sd(sp_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, p, true_r_sp,
              if (is.null(supposed_r_sp)) NA_real_ else supposed_r_sp,
              width, conf_level,
              true_R2_full,
              if (is.null(supposed_R2_full)) NA_real_ else supposed_R2_full,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  .ss_sens_record(out, G, !is.na(ci_width), conf_level)
}
