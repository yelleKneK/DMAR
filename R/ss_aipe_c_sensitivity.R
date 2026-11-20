# Sensitivity analysis for AIPE on an unstandardized contrast of means.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Unstandardized Contrast
#'
#' @description
#' Quantifies how much misspecification of the population error variance
#' distorts an AIPE-based sample size plan for an unstandardized contrast
#' of means. Because the half-width of the confidence interval on
#' \eqn{\psi = \sum_j c_j \mu_j} depends on the error variance, the
#' contrast weights, and the per-group sample size, but not on the
#' value of \eqn{\psi} itself, this sensitivity analysis varies the
#' planning value of the error variance. On each replication the
#' function simulates \eqn{n} observations per group from a normal
#' population with variance \code{true_var_error}, builds the
#' confidence interval via \code{\link{ci_c}}, and summarizes the
#' realized widths and coverage of \code{true_contrast}.
#'
#' @param true_var_error Population error variance (the
#'   data generating value). Must be positive.
#' @param supposed_var_error Error variance used to plan the study
#'   (the value passed to \code{\link{ss_aipe_c}}). Supply this or
#'   \code{specified_n_per_group} but not both.
#' @param c_weights Contrast weight vector. Must sum to zero.
#' @param width Desired full width of the confidence interval on the
#'   unstandardized contrast.
#' @param true_contrast Population value of the contrast; the simulator places
#'   group means such that \eqn{\sum_j c_j \mu_j = }\code{true_contrast}. The
#'   width of the interval does not depend on this value but the realized
#'   coverage of \code{true_contrast} does. Default \code{0}.
#' @param specified_n_per_group Per-group sample size to evaluate (incompatible
#'   with \code{supposed_var_error}); when used, the planner is
#'   bypassed.
#' @param conf_level Confidence level (default \code{0.95}). A level above .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param assurance Optional probability that the realized interval is
#'   no wider than \code{width}; passed to \code{\link{ss_aipe_c}} when
#'   resolving the planned sample size.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print the iteration index after each
#'   replication (helpful for long runs); default \code{FALSE}.
#' @param filename Optional path of a CSV file to receive the
#'   per-replication results (the contrast estimate, the confidence limits,
#'   the interval width, and the two tail misses), appended when the file
#'   already exists and created otherwise; the default \code{NULL} writes
#'   nothing, and a throwaway run that wants the file should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized contrast (\code{mean_contrast}, \code{median_contrast},
#'   \code{sd_contrast}) and of the interval width, the proportion of
#'   intervals at or below \code{width}, the tail-specific and overall
#'   empirical non-coverage of \code{true_contrast}, and the input echoes
#'   (per-group sample size, total sample size, the true and supposed
#'   error variances \code{true_var_error} and \code{supposed_var_error},
#'   \code{true_contrast}, width, confidence level, and, when one was
#'   supplied, assurance).
#'
#' @references
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_c}}, \code{\link{ci_c}}, \code{\link{ss_aipe_sc_sensitivity}}
#'
#' @examples
#' # G = 50 keeps the example quick; a reported analysis deserves the
#' # default G of 1000. First the well-specified case: the planner used an
#' # error variance of 4, and the truth is 4.
#' set.seed(113)
#' ss_aipe_c_sensitivity(
#'   true_var_error     = 4,
#'   supposed_var_error = 4,
#'   c_weights = c(-1, 0, 1),
#'   width = 1, G = 50, print_iter = FALSE
#' )
#'
#' # Misspecified: planner used 4, truth is 9. Realized widths inflate.
#' set.seed(113)
#' ss_aipe_c_sensitivity(
#'   true_var_error     = 9,
#'   supposed_var_error = 4,
#'   c_weights = c(-1, 0, 1),
#'   width = 1, G = 50, print_iter = FALSE
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
ss_aipe_c_sensitivity <- function(true_var_error = NULL,
                                  supposed_var_error = NULL,
                                  c_weights, width,
                                  true_contrast = 0,
                                  specified_n_per_group = NULL,
                                  conf_level = 0.95,
                                  assurance = NULL,
                                  G = 1000, print_iter = FALSE,
                                  filename = NULL,
                                  seed = NULL) {
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_numbers("ss_aipe_c_sensitivity", "true_contrast")
  if (is.null(supposed_var_error) && is.null(specified_n_per_group))
    .dmar_refuse_argument(
      "You must specify either 'supposed_var_error' or 'specified_n_per_group'.",
      arg = c("supposed_var_error", "specified_n_per_group"))
  if (!is.null(supposed_var_error) && !is.null(specified_n_per_group))
    .dmar_refuse_argument(
      "You must specify 'supposed_var_error' or 'specified_n_per_group', but not both.",
      arg = c("supposed_var_error", "specified_n_per_group"))
  if (is.null(true_var_error) || !is.numeric(true_var_error) || true_var_error <= 0)
    .dmar_refuse_argument(
      "'true_var_error' must be a single positive number.",
      arg = "true_var_error")
  if (!is.numeric(c_weights) || abs(sum(c_weights)) > 1e-8)
    .dmar_refuse_argument(
      "'c_weights' must be numeric and sum to zero.", arg = "c_weights")
  # All-zero weights define no contrast; they gave a summary of zeros.
  .dmar_check_contrast_nonzero(c_weights)
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single positive number.", arg = "width")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  J <- length(c_weights)

  if (!is.null(supposed_var_error)) {
    plan <- ss_aipe_c(true_var_error = supposed_var_error,
                      c_weights = c_weights, width = width,
                      conf_level = conf_level, assurance = assurance)
    n <- plan$value[plan$term == "necessary_n_per_group"]
  } else {
    n <- specified_n_per_group
  }
  n <- as.integer(n)
  if (n < 2) .dmar_refuse_argument(
    "Resolved sample size per group is < 2.", arg = NULL)

  # Place group means so that sum(c_weights * mu) = true_contrast, with mu_j = 0
  # for any group whose weight is zero. Distribute the contrast across the
  # nonzero-weight groups via the minimum-norm solution mu = c * true_contrast / sum(c^2).
  ssq <- sum(c_weights^2)
  mu  <- if (ssq == 0) rep(0, J) else c_weights * true_contrast / ssq

  psi_hat   <- numeric(G)
  ci_lo     <- numeric(G)
  ci_hi     <- numeric(G)
  ci_width  <- numeric(G)
  tI_lower  <- logical(G)
  tI_upper  <- logical(G)

  sd_pop <- sqrt(true_var_error)
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    y <- unlist(lapply(seq_len(J), function(j)
      stats::rnorm(n, mean = mu[j], sd = sd_pop)))
    grp <- factor(rep(seq_len(J), each = n))
    means <- tapply(y, grp, mean)
    s_anova <- sqrt(stats::var(stats::residuals(stats::lm(y ~ grp))) *
                    (n * J - 1) / (n * J - J))

    ci <- ci_c(means = as.numeric(means), s_anova = s_anova,
               c_weights = c_weights, n_per_group = rep(n, J), N = n * J,
               conf_level = conf_level)
    lo <- ci$value[ci$term == "ci_lower"]
    hi <- ci$value[ci$term == "ci_upper"]
    psi_hat[g]  <- ci$value[ci$term == "contrast"]
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_contrast < lo
    tI_upper[g] <- true_contrast > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(contrast = psi_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_contrast", "median_contrast", "sd_contrast",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_per_group", "total_N",
              "true_var_error", "supposed_var_error",
              "true_contrast", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(psi_hat, na.rm = TRUE), stats::median(psi_hat, na.rm = TRUE),
              stats::sd(psi_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE), stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE), mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, n * J,
              true_var_error,
              if (is.null(supposed_var_error)) NA_real_ else supposed_var_error,
              true_contrast, width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(ci_width), conf_level)
}

# The one place the sensitivity family validates its `filename` argument:
# NULL means nothing is written; anything else must be a single non-empty
# string. Called before any replication runs so a bad path fails at once.
.check_filename <- function(filename) {
  if (is.null(filename)) return(invisible(NULL))
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) ||
      !nzchar(filename)) {
    .dmar_refuse_argument(
      "'filename' must be NULL or a single character string naming the file to write.",
      arg = "filename")
  }
  invisible(filename)
}

# Shared writer for the ss_aipe_*_sensitivity() family when a filename is
# supplied: the per-replication rows are appended without a header when the
# file already holds a table and written with a header otherwise. Base R
# only, so no Suggests dependency is needed at write time.
.write_sensitivity_csv <- function(per_rep, filename) {
  suppressWarnings(file_exist <- try(utils::read.csv(filename), silent = TRUE))
  if (!is.null(dim(file_exist))) {
    utils::write.table(per_rep, filename, sep = ",", row.names = FALSE,
                       col.names = FALSE, append = TRUE)
  } else {
    utils::write.table(per_rep, filename, sep = ",", row.names = FALSE,
                       append = FALSE)
  }
}
