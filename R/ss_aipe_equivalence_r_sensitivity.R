# Sensitivity analysis for AIPE on the equivalence-test (TOST) correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test Correlation
#'
#' @description
#' Quantifies how much misspecification of the population correlation
#' distorts an AIPE-based sample size plan for the two-one-sided-tests
#' (TOST) confidence interval on the Pearson correlation. On each
#' replication the function simulates \emph{N} bivariate normal pairs
#' with population correlation \code{true_r}, computes the sample
#' correlation and its Fisher's \eqn{Z} confidence interval via
#' \code{\link{ci_r}}, and summarizes the realized widths and the
#' proportion of replications in which the computed interval falls
#' entirely inside (\code{equivalent}) the specified equivalence
#' bounds.
#'
#' @param true_r Population correlation (the data generating value).
#'   Defaults to \code{0} (no association, the exact equivalence
#'   case).
#' @param estimated_r Planning value of the population correlation
#'   passed to \code{\link{ss_aipe_equivalence_r}}; supply this or
#'   \code{specified_N} but not both.
#' @param width Desired full width of the two-sided CI on the
#'   correlation.
#' @param rho_lower,rho_upper Equivalence bounds on the correlation,
#'   as positive magnitudes with the same meaning as in
#'   \code{\link{equivalence_r}}: the region is
#'   \eqn{(-\rho_L, +\rho_U)}. \code{rho_upper} is required;
#'   \code{rho_lower} defaults to \code{rho_upper} (a symmetric
#'   region). The simulator records whether the realized CI falls
#'   entirely inside the region.
#' @param specified_N Sample size to evaluate.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized correlation and CI width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of \code{true_r}, the proportion of intervals
#'   classified as \code{equivalent} (CI fully inside the bounds),
#'   and the input echoes, including \code{assurance} (present only
#'   when an assurance was supplied).
#'
#' @references
#' Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for
#'   comparing correlation and regression coefficients. \emph{British
#'   Journal of Mathematical and Statistical Psychology, 68}(2),
#'   292--309. \doi{10.1111/bmsp.12045}
#'
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
#'   association: An equivalence testing approach. \emph{British
#'   Journal of Mathematical and Statistical Psychology, 63}(3),
#'   527--537. \doi{10.1348/000711009X475853}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_equivalence_r}}, \code{\link{equivalence_r}},
#'   \code{\link{ss_aipe_r_sensitivity}},
#'   \code{\link{ss_aipe_equivalence_smd_sensitivity}}
#'
#' @examples
#' # Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
#' set.seed(113)
#' ss_aipe_equivalence_r_sensitivity(
#'   true_r      = 0.0,
#'   estimated_r = 0.0,
#'   width       = 0.30,
#'   rho_upper   = 0.20,
#'   G = 50, print_iter = FALSE
#' )
#'
#' @keywords design htest
#'
#' @family AIPE sample size planning
#'
#' @export
ss_aipe_equivalence_r_sensitivity <- function(true_r = 0,
                                       estimated_r = NULL,
                                       width,
                                       rho_lower = NULL,
                                       rho_upper = NULL,
                                       specified_N = NULL,
                                       conf_level = 0.95,
                                       assurance = NULL,
                                       G = 1000, print_iter = FALSE,
                                       save = FALSE,
                                       filename = "ss_aipe_equivalence_r_sensitivity_result.csv") {
  if (is.null(estimated_r) && is.null(specified_N))
    stop("You must specify either 'estimated_r' or 'specified_N'.",
         call. = FALSE)
  if (!is.null(estimated_r) && !is.null(specified_N))
    stop("You must specify 'estimated_r' or 'specified_N', but not both.",
         call. = FALSE)
  if (!is.numeric(true_r) || length(true_r) != 1L || abs(true_r) >= 1)
    stop("'true_r' must be a single correlation with |r| < 1.",
         call. = FALSE)

  if (is.null(rho_upper))
    stop("'rho_upper' must be specified.", call. = FALSE)
  if (!is.numeric(rho_upper) || rho_upper <= 0 || rho_upper >= 1)
    stop("'rho_upper' must be in (0, 1).", call. = FALSE)
  if (is.null(rho_lower)) rho_lower <- rho_upper
  if (!is.numeric(rho_lower) || rho_lower <= 0 || rho_lower >= 1)
    stop("'rho_lower' must be in (0, 1).", call. = FALSE)

  # ss_aipe_equivalence_r() plans the width of a (1 - 2 * alpha) confidence
  # interval, and the simulation below evaluates ci_r() at conf_level, so
  # the two match only when alpha = (1 - conf_level) / 2 (the same rule
  # ss_aipe_equivalence_smd_sensitivity() follows).
  alpha <- (1 - conf_level) / 2
  if (!is.null(estimated_r)) {
    plan <- ss_aipe_equivalence_r(population_r = estimated_r, width = width,
                           alpha_level = alpha, assurance = assurance)
    N <- plan$value[plan$term == "necessary_N"]
  } else {
    N <- specified_N
  }
  N <- as.integer(N)
  if (N < 4) stop("Resolved sample size is < 4; the Fisher's Z interval ",
                  "is based on the variance 1/(N - 3).", call. = FALSE)

  r_hat      <- numeric(G)
  ci_lo      <- numeric(G)
  ci_hi      <- numeric(G)
  ci_width   <- numeric(G)
  equivalent <- logical(G)
  tI_lower   <- logical(G)
  tI_upper   <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    x <- stats::rnorm(N)
    y <- true_r * x + sqrt(1 - true_r^2) * stats::rnorm(N)
    r <- stats::cor(x, y)

    ci <- ci_r(r = r, n = N, conf_level = conf_level)
    lo <- ci$value[ci$term == "lower_limit"]
    hi <- ci$value[ci$term == "upper_limit"]
    r_hat[g]      <- r
    ci_lo[g]      <- lo
    ci_hi[g]      <- hi
    ci_width[g]   <- hi - lo
    equivalent[g] <- (lo > -rho_lower) && (hi < rho_upper)
    tI_lower[g]   <- true_r < lo
    tI_upper[g]   <- true_r > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(r = r_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          equivalent = equivalent,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_r", "median_r", "sd_r",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w", "pct_equivalent",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N",
              "true_r", "estimated_r",
              "width", "conf_level",
              "rho_lower", "rho_upper",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(r_hat, na.rm = TRUE),
              stats::median(r_hat, na.rm = TRUE),
              stats::sd(r_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(equivalent, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              N,
              true_r,
              if (is.null(estimated_r)) NA_real_ else estimated_r,
              width, conf_level,
              -rho_lower, rho_upper,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
