# Sensitivity analysis for AIPE on a partial correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Partial Correlation
#'
#' @description
#' Simulates a sample size plan for a partial correlation to show what it
#' delivers when the population value differs from the planning value, or
#' at a sample size chosen in advance. It takes the population partial
#' correlation, either a planning value or a sample size, the number of
#' predictors (\emph{X} and the variables partialed out), and a target
#' width, and returns summaries of the estimate, the interval width, the proportion of intervals no wider
#' than the target, and coverage over \code{G} simulated samples.
#'
#' @param true_cor Population partial correlation \eqn{\rho} between \emph{Y} and
#'   \emph{X} controlling for \eqn{Z_1, \ldots, Z_{p-1}}; must lie in
#'   \eqn{(-1, 1)}.
#' @param supposed_cor Planning value of the partial correlation
#'   passed to \code{\link{ss_aipe_partial_r}}; supply this or
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
#' @param width Desired full width of the CI on the partial
#'   correlation.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{supposed_cor}); a whole number of at least \eqn{p + 2}, or
#'   \eqn{p + 3} with \code{fisher_z = TRUE}.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to
#'   \code{\link{ss_aipe_partial_r}}; used only with
#'   \code{supposed_cor}.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print iteration index per replication.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample partial correlation, the two confidence
#'   limits, the interval width, and two indicators of whether the
#'   interval missed \code{true_cor} below or above): nothing is
#'   written when \code{filename} is \code{NULL} (the default), a new
#'   file with a header row is created otherwise, an existing file at
#'   that path is appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#'
#' @param fisher_z Logical. The interval that is planned and simulated.
#'   \code{FALSE} (the default) is the raw-scale interval,
#'   \eqn{r \pm z_{1 - \alpha/2} (1 - r^2) / \sqrt{N - p}};
#'   \code{TRUE} is the interval formed on Fisher's \emph{Z} scale with
#'   standard error \eqn{1/\sqrt{N - p - 2}} and back-transformed. The
#'   value is passed to \code{\link{ss_aipe_partial_r}}, so the interval
#'   simulated is always the one planned.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the run reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' \strong{Estimand and interval.} The population partial correlation
#' between \emph{Y} and \emph{X} given \eqn{Z_1, \ldots, Z_{p-1}}, estimated
#' from \emph{N} independent multivariate normal observations by the
#' correlation of the two sets of least squares residuals. The interval is
#' two-sided at \code{conf_level} and is the one \code{fisher_z} names.
#' The raw-scale limits are not truncated at \eqn{-1} and \eqn{1}, so that
#' the widths are those the planner plans; with a large \eqn{|\rho|} and a
#' small sample a limit can pass 1 in absolute value.
#'
#' \strong{What is simulated.} The function constructs a
#' \eqn{(p + 1) \times (p + 1)} population covariance matrix for
#' \emph{Y}, \emph{X}, and \eqn{Z_1, \ldots, Z_{p-1}}, all standardized and
#' uncorrelated but for a correlation of \code{true_cor} between \emph{Y}
#' and \emph{X}, so that the partial correlation is \code{true_cor}. Under
#' normality the distribution of the sample partial correlation depends on
#' the population only through that value, \eqn{p}, and \emph{N}, so this
#' population stands for every normal population with the same partial
#' correlation. With \code{supposed_cor}, \emph{N} is the plan of
#' \code{\link{ss_aipe_partial_r}} for that value (with \code{assurance}
#' and \code{fisher_z} passed on); with \code{specified_N} it is that
#' number. Each of \code{G} replications draws \emph{N} rows, computes the
#' sample partial correlation, and forms the interval. The \eqn{p - 1}
#' variables partialed out are the count \code{\link{ss_aipe_partial_r}}
#' plans for and \code{\link{var_partial_r}} uses. A sample with nothing
#' partialed out is the case of \code{\link{ss_aipe_r_sensitivity}}.
#'
#' \strong{Reading the result.} \code{pct_ci_less_width} is the proportion of
#' the \code{G} intervals no wider than \code{width}: about one half for
#' a plan without \code{assurance} whose planning value is correct, and at
#' least \code{assurance}, within Monte Carlo error, for a plan with it.
#' The Monte Carlo standard error of a proportion \eqn{A} is
#' \eqn{\sqrt{A (1 - A) / G}}. \code{pct_ci_miss_low} counts intervals
#' that lay entirely above \code{true_cor} and \code{pct_ci_miss_high}
#' those entirely below; their sum is \code{total_type_I_error}. Every
#' replication has an interval, so all rows are over \code{G}
#' replications.
#'
#' @return A \code{data.frame} with the columns \code{term} and
#'   \code{value}: \code{mean_partial_r}, \code{median_partial_r}, and
#'   \code{sd_partial_r} (the sample partial correlations);
#'   \code{mean_ci_width}, \code{median_ci_width}, and \code{sd_ci_width}
#'   (the full widths of the intervals); \code{pct_ci_less_width} (the
#'   proportion of intervals no wider than \code{width});
#'   \code{pct_ci_miss_low}, \code{pct_ci_miss_high}, and
#'   \code{total_type_I_error} (non-coverage of \code{true_cor}, by side
#'   and in total); and the echoes \code{total_N} (the sample size
#'   simulated), \code{p}, \code{true_cor}, \code{supposed_cor}
#'   (\code{NA} when \code{specified_N} was given), \code{width},
#'   \code{conf_level}, and \code{assurance} (present only when an
#'   assurance was supplied). The proportions are on the 0 to 1 scale. The
#'   attribute \code{interval} is \code{"raw"} or \code{"fisher_z"}, and
#'   the printed table names the interval beneath it.
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
#' @seealso \code{\link{ss_aipe_partial_r}}, \code{\link{ss_aipe_semipartial_r_sensitivity}}
#'
#' @examples
#' # A plan for a full width of .40 at a planning value of .40 that is
#' # correct. G = 50 keeps the example quick; a reported check deserves
#' # the default G = 1000 or more.
#' ss_aipe_partial_r_sensitivity(
#'   true_cor = 0.40, supposed_cor = 0.40, p = 4, width = 0.40,
#'   G = 50, seed = 113
#' )
#'
#' # The same plan when the population value is smaller than planned: the
#' # intervals are wider and fewer meet the target.
#' ss_aipe_partial_r_sensitivity(
#'   true_cor = 0.20, supposed_cor = 0.40, p = 4, width = 0.40,
#'   G = 50, seed = 113
#' )
#'
#' # Fisher's Z interval, planned and simulated.
#' ss_aipe_partial_r_sensitivity(
#'   true_cor = 0.40, supposed_cor = 0.40, p = 4, width = 0.40,
#'   G = 50, fisher_z = TRUE, seed = 113
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
ss_aipe_partial_r_sensitivity <- function(true_cor = NULL,
                                          supposed_cor = NULL,
                                          p, width,
                                          specified_N = NULL,
                                          conf_level = 0.95,
                                          assurance = NULL,
                                          G = 1000, print_iter = FALSE,
                                          filename = NULL,
                                          fisher_z = FALSE,
                                          seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  if (missing(p))
    .dmar_refuse_argument(
      "'p' (the number of predictors) is required.", arg = "p")
  p <- .check_partial_p(p)
  # The arithmetic below is written in J, the number of variables partialed
  # out: N - J - 1 = N - p, N - J - 3 = N - p - 2, and the population has
  # J + 2 = p + 1 variables.
  J <- p - 1
  if (is.null(supposed_cor) && is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify either 'supposed_cor' or 'specified_N'.",
      arg = c("supposed_cor", "specified_N"))
  if (!is.null(supposed_cor) && !is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify 'supposed_cor' or 'specified_N', but not both.",
      arg = c("supposed_cor", "specified_N"))
  if (is.null(true_cor) || !is.numeric(true_cor) || length(true_cor) != 1L ||
      is.na(true_cor) || abs(true_cor) >= 1)
    .dmar_refuse_argument(
      "'true_cor' must be a single value in (-1, 1).", arg = "true_cor")
  G <- .dmar_check_count(G, "G")
  .check_filename(filename)
  if (!is.logical(fisher_z) || length(fisher_z) != 1L || is.na(fisher_z))
    .dmar_refuse_argument("'fisher_z' must be TRUE or FALSE.", arg = "fisher_z")

  # The raw-scale interval needs n - J - 1 > 0 and the residuals a degree
  # of freedom to spare, so n is at least J + 3; Fisher's Z interval needs
  # n - J - 3 > 0, so n is at least J + 4.
  n_min <- J + if (fisher_z) 4L else 3L
  if (!is.null(supposed_cor)) {
    plan <- ss_aipe_partial_r(true_cor = supposed_cor, p = p, width = width,
                              conf_level = conf_level, fisher_z = fisher_z,
                              assurance = assurance)
    n <- as.integer(plan$value[plan$term == "necessary_N"])
    if (n < n_min)
      .dmar_refuse_argument(c(
        "The planned sample size (", n, ") is below p + ", n_min - p,
        ", the smallest at which the ",
        if (fisher_z) "Fisher's Z" else "raw-scale",
        " interval is simulated."), arg = NULL)
  } else {
    n <- .check_whole_n(specified_N, "specified_N", n_min)
  }

  # Population of (Y, X, Z_1, ..., Z_J), all standardized and uncorrelated
  # but for Cor(Y, X) = true_cor. With the Z's unrelated to both, the
  # partial correlation of Y and X given Z_1, ..., Z_J is true_cor.
  Sigma <- diag(J + 2)
  Sigma[1, 2] <- Sigma[2, 1] <- true_cor

  if (!requireNamespace("MASS", quietly = TRUE))
    .dmar_refuse(
      "Package 'MASS' is required for the multivariate normal sampling step.",
      "package_required", package = "MASS")

  .dmar_local_seed(seed)

  pr_hat   <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_width <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    dat <- MASS::mvrnorm(n = n, mu = rep(0, J + 2), Sigma = Sigma)
    Y  <- dat[, 1]; X1 <- dat[, 2]
    Xrest <- dat[, 3:(J + 2), drop = FALSE]
    e_y  <- stats::residuals(stats::lm(Y  ~ Xrest))
    e_x  <- stats::residuals(stats::lm(X1 ~ Xrest))
    r_pr <- stats::cor(e_y, e_x)
    if (fisher_z) {
      # Fisher's Z interval for the partial correlation with J variables
      # partialed out: the variance of Z is 1 / (n - J - 3).
      z   <- atanh(r_pr)
      se  <- 1 / sqrt(n - J - 3)
      lo  <- tanh(z - crit * se)
      hi  <- tanh(z + crit * se)
    } else {
      # The raw-scale interval from the first-order variance
      # (1 - r^2)^2 / (n - J - 1); the limits are not truncated.
      half <- crit * (1 - r_pr^2) / sqrt(n - J - 1)
      lo  <- r_pr - half
      hi  <- r_pr + half
    }
    pr_hat[g]   <- r_pr
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_cor < lo
    tI_upper[g] <- true_cor > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(partial_r = pr_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_partial_r", "median_partial_r", "sd_partial_r",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "p",
              "true_cor", "supposed_cor", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(pr_hat, na.rm = TRUE),
              stats::median(pr_hat, na.rm = TRUE),
              stats::sd(pr_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, p, true_cor,
              if (is.null(supposed_cor)) NA_real_ else supposed_cor,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      subclass = "dmar_ss_aipe_partial_r_sensitivity")
  attr(out, "interval") <- if (fisher_z) "fisher_z" else "raw"
  .ss_sens_record(out, G, !is.na(ci_width), conf_level)
}

# The line beneath a table of ss_aipe_partial_r_sensitivity() that names
# the interval the replications formed (the `interval` attribute), added to
# the lines every table carries (.dmar_tbl_footer() in
# R/dmar_output_helpers.R).
#' @exportS3Method .dmar_tbl_footer dmar_ss_aipe_partial_r_sensitivity
#' @noRd
.dmar_tbl_footer.dmar_ss_aipe_partial_r_sensitivity <- function(x, ...) {
  interval <- attr(x, "interval", exact = TRUE)
  line <- if (identical(interval, "fisher_z")) {
    "Interval: Fisher's Z scale, back-transformed (standard error 1 / sqrt(N - p - 2))"
  } else if (identical(interval, "raw")) {
    "Interval: raw scale, r plus or minus z (1 - r^2) / sqrt(N - p)"
  }
  .dmar_footer_add(NextMethod(), line, "interval")
}
