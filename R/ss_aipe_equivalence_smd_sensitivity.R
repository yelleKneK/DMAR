# Sensitivity analysis for AIPE on the equivalence-test (TOST) SMD.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test SMD
#'
#' @description
#' Quantifies how much misspecification of the population standardized
#' mean difference distorts an AIPE-based sample size plan for the
#' two-one-sided-tests (TOST) confidence interval on the SMD. On each
#' replication the function simulates two normal groups of size \emph{n}
#' per group with population standardized mean difference
#' \code{true_smd}, computes the SMD and its noncentral \emph{t}
#' confidence interval via \code{\link{ci_smd}}, and summarizes the
#' realized widths and the proportion of replications in which the
#' computed interval falls entirely inside (\code{equivalent}) the
#' specified equivalence bounds.
#'
#' @param true_smd Population standardized mean difference (the
#'   data generating value). Defaults to \code{0} (perfect
#'   equivalence).
#' @param estimated_smd Planning value of the population SMD passed to
#'   \code{\link{ss_aipe_equivalence_smd}}; supply this or \code{n_per_group}
#'   but not both.
#' @param width Desired full width of the two-sided CI on the SMD.
#' @param delta_lower,delta_upper Equivalence bounds on the SMD, as
#'   positive magnitudes with the same meaning as in
#'   \code{\link{equivalence_smd}}: the region is
#'   \eqn{(-\code{delta_lower}, \code{delta_upper})}.
#'   \code{delta_upper} is required; \code{delta_lower} defaults to
#'   \code{delta_upper} (a symmetric region). The simulator records
#'   whether the realized CI falls entirely inside the region.
#' @param n_per_group Per-group sample size to evaluate.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized SMD and CI width, the proportion of intervals at or
#'   below \code{width}, tail-specific and overall non-coverage of
#'   \code{true_smd}, the proportion of intervals classified as
#'   \code{equivalent} (CI fully inside the bounds), and the input
#'   echoes, including \code{assurance} (present only when an
#' assurance was supplied).
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11},
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_equivalence_smd}}, \code{\link{equivalence_smd}}, \code{\link{ss_aipe_smd_sensitivity}}
#'
#' @examples
#' # Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
#' set.seed(113)
#' ss_aipe_equivalence_smd_sensitivity(
#'   true_smd      = 0.0,
#'   estimated_smd = 0.0,
#'   width         = 0.30,
#'   delta_upper   = 0.20,
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
ss_aipe_equivalence_smd_sensitivity <- function(true_smd = 0,
                                         estimated_smd = NULL,
                                         width,
                                         delta_lower = NULL,
                                         delta_upper = NULL,
                                         n_per_group = NULL,
                                         conf_level = 0.95,
                                         assurance = NULL,
                                         G = 1000, print_iter = FALSE,
                                         save = FALSE,
                                         filename = "ss_aipe_equivalence_smd_sensitivity_result.csv") {
  if (is.null(estimated_smd) && is.null(n_per_group))
    stop("You must specify either 'estimated_smd' or 'n_per_group'.", call. = FALSE)
  if (!is.null(estimated_smd) && !is.null(n_per_group))
    stop("You must specify 'estimated_smd' or 'n_per_group', but not both.", call. = FALSE)
  if (!is.numeric(true_smd) || length(true_smd) != 1L)
    stop("'true_smd' must be a single numeric value.", call. = FALSE)

  if (is.null(delta_upper))
    stop("'delta_upper' must be specified (the upper equivalence bound).",
         call. = FALSE)
  if (!is.numeric(delta_upper) || delta_upper <= 0)
    stop("'delta_upper' must be a positive number.", call. = FALSE)
  if (is.null(delta_lower)) delta_lower <- delta_upper
  if (!is.numeric(delta_lower) || delta_lower <= 0)
    stop("'delta_lower' must be a positive number.", call. = FALSE)

  # ss_aipe_equivalence_smd() plans the width of a (1 - 2 * alpha) confidence interval,
  # and the simulation below evaluates ci_smd() at conf_level, so the two match
  # only when alpha = (1 - conf_level) / 2. Using 1 - conf_level planned a 90%
  # width while the 95% interval was scored, underplanning the sample size.
  alpha <- (1 - conf_level) / 2
  if (!is.null(estimated_smd)) {
    plan <- ss_aipe_equivalence_smd(population_smd = estimated_smd, width = width,
                             alpha_level = alpha, assurance = assurance,
                             balanced = TRUE)
    n_per <- plan$value[plan$term == "necessary_n_per_group"]
  } else {
    n_per <- n_per_group
  }
  n_per <- as.integer(n_per)
  if (n_per < 2) stop("Resolved per-group sample size is < 2.", call. = FALSE)

  smd_hat    <- numeric(G)
  ci_lo      <- numeric(G)
  ci_hi      <- numeric(G)
  ci_width   <- numeric(G)
  equivalent <- logical(G)
  tI_lower   <- logical(G)
  tI_upper   <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    x1 <- stats::rnorm(n_per, mean = true_smd, sd = 1)
    x2 <- stats::rnorm(n_per, mean = 0,        sd = 1)
    m_diff <- mean(x1) - mean(x2)
    s_pool <- sqrt(((n_per - 1) * stats::var(x1) +
                    (n_per - 1) * stats::var(x2)) / (2 * n_per - 2))
    d <- m_diff / s_pool

    ci <- ci_smd(smd = d, n_1 = n_per, n_2 = n_per, conf_level = conf_level)
    lo <- ci$value[ci$term == "lower_limit"]
    hi <- ci$value[ci$term == "upper_limit"]
    smd_hat[g]    <- d
    ci_lo[g]      <- lo
    ci_hi[g]      <- hi
    ci_width[g]   <- hi - lo
    equivalent[g] <- (lo > -delta_lower) && (hi < delta_upper)
    tI_lower[g]   <- true_smd < lo
    tI_upper[g]   <- true_smd > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(smd = smd_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          equivalent = equivalent,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_smd", "median_smd", "sd_smd",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w", "pct_equivalent",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_per_group", "total_N",
              "true_smd", "estimated_smd",
              "width", "conf_level",
              "delta_lower", "delta_upper",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(smd_hat, na.rm = TRUE),
              stats::median(smd_hat, na.rm = TRUE),
              stats::sd(smd_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(equivalent, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n_per, 2 * n_per,
              true_smd,
              if (is.null(estimated_smd)) NA_real_ else estimated_smd,
              width, conf_level,
              -delta_lower, delta_upper,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
