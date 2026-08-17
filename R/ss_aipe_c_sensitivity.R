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
#' population with variance \code{true_error_variance}, builds the
#' confidence interval via \code{\link{ci_c}}, and summarizes the
#' realized widths and coverage of \code{true_psi}.
#'
#' @param true_error_variance Population error variance (the
#'   data generating value). Must be positive.
#' @param estimated_error_variance Error variance used to plan the study
#'   (the value passed to \code{\link{ss_aipe_c}}). Supply this or
#'   \code{n_per_group} but not both.
#' @param c_weights Contrast weight vector. Must sum to zero.
#' @param width Desired full width of the confidence interval on the
#'   unstandardized contrast.
#' @param true_psi Population value of the contrast; the simulator places
#'   group means such that \eqn{\sum_j c_j \mu_j = }\code{true_psi}. The
#'   width of the interval does not depend on this value but the realized
#'   coverage of \code{true_psi} does. Default \code{0}.
#' @param n_per_group Per-group sample size to evaluate (incompatible
#'   with \code{estimated_error_variance}); when used, the planner is
#'   bypassed.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional probability that the realized interval is
#'   no wider than \code{width}; passed to \code{\link{ss_aipe_c}} when
#'   resolving the planned sample size.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print the iteration index after each
#'   replication (helpful for long runs); default \code{FALSE}.
#' @param save Logical. If \code{TRUE} the per-replication results are
#'   appended to \code{filename}; default \code{FALSE}.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized estimator and interval width, the proportion of
#'   intervals at or below \code{width}, the tail-specific and overall
#'   empirical non-coverage of \code{true_psi}, and the input echoes
#'   (per-group sample size, total sample size, true and estimated
#'   error variances, width, confidence level, and, when one was
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
#' # Monte Carlo sweep; G is small here so the example runs quickly.
#' # Well-specified: planner used error_variance = 4, truth is 4.
#' set.seed(113)
#' ss_aipe_c_sensitivity(
#'   true_error_variance      = 4,
#'   estimated_error_variance = 4,
#'   c_weights = c(-1, 0, 1),
#'   width = 1, G = 50, print_iter = FALSE
#' )
#'
#' # Misspecified: planner used 4, truth is 9. Realized widths inflate.
#' set.seed(113)
#' ss_aipe_c_sensitivity(
#'   true_error_variance      = 9,
#'   estimated_error_variance = 4,
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
#' @export
ss_aipe_c_sensitivity <- function(true_error_variance = NULL,
                                  estimated_error_variance = NULL,
                                  c_weights, width,
                                  true_psi = 0,
                                  n_per_group = NULL,
                                  conf_level = 0.95,
                                  assurance = NULL,
                                  G = 1000, print_iter = FALSE,
                                  save = FALSE,
                                  filename = "ss_aipe_c_sensitivity_result.csv") {
  if (is.null(estimated_error_variance) && is.null(n_per_group))
    stop("You must specify either 'estimated_error_variance' or 'n_per_group'.", call. = FALSE)
  if (!is.null(estimated_error_variance) && !is.null(n_per_group))
    stop("You must specify 'estimated_error_variance' or 'n_per_group', but not both.", call. = FALSE)
  if (is.null(true_error_variance) || !is.numeric(true_error_variance) || true_error_variance <= 0)
    stop("'true_error_variance' must be a single positive number.", call. = FALSE)
  if (!is.numeric(c_weights) || abs(sum(c_weights)) > 1e-8)
    stop("'c_weights' must be numeric and sum to zero.", call. = FALSE)
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.", call. = FALSE)

  J <- length(c_weights)

  if (!is.null(estimated_error_variance)) {
    plan <- ss_aipe_c(error_variance = estimated_error_variance,
                      c_weights = c_weights, width = width,
                      conf_level = conf_level, assurance = assurance)
    n <- plan$value[plan$term == "necessary_n_per_group"]
  } else {
    n <- n_per_group
  }
  n <- as.integer(n)
  if (n < 2) stop("Resolved sample size per group is < 2.", call. = FALSE)

  # Place group means so that sum(c_weights * mu) = true_psi, with mu_j = 0
  # for any group whose weight is zero. Distribute the contrast across the
  # nonzero-weight groups via the minimum-norm solution mu = c * true_psi / sum(c^2).
  ssq <- sum(c_weights^2)
  mu  <- if (ssq == 0) rep(0, J) else c_weights * true_psi / ssq

  psi_hat   <- numeric(G)
  ci_lo     <- numeric(G)
  ci_hi     <- numeric(G)
  ci_width  <- numeric(G)
  tI_lower  <- logical(G)
  tI_upper  <- logical(G)

  sd_pop <- sqrt(true_error_variance)
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    y <- unlist(lapply(seq_len(J), function(j)
      stats::rnorm(n, mean = mu[j], sd = sd_pop)))
    grp <- factor(rep(seq_len(J), each = n))
    means <- tapply(y, grp, mean)
    s_anova <- sqrt(stats::var(stats::residuals(stats::lm(y ~ grp))) *
                    (n * J - 1) / (n * J - J))

    ci <- ci_c(means = as.numeric(means), s_anova = s_anova,
               c_weights = c_weights, n = rep(n, J), N = n * J,
               conf_level = conf_level)
    lo <- ci$value[ci$term == "lower_limit"]
    hi <- ci$value[ci$term == "upper_limit"]
    psi_hat[g]  <- ci$value[ci$term == "contrast"]
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_psi < lo
    tI_upper[g] <- true_psi > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(psi_hat = psi_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_psi", "median_psi", "sd_psi",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_per_group", "total_N",
              "true_error_variance", "estimated_error_variance",
              "true_psi", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(psi_hat, na.rm = TRUE), stats::median(psi_hat, na.rm = TRUE),
              stats::sd(psi_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE), stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE), mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, n * J,
              true_error_variance,
              if (is.null(estimated_error_variance)) NA_real_ else estimated_error_variance,
              true_psi, width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}

# Shared writer for *_sensitivity CSV saves. Uses base R so no
# Suggests dependency is needed at write time.
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
