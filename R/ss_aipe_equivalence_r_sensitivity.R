# Sensitivity analysis for AIPE on the equivalence-test (TOST) correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test Correlation
#'
#' @description
#' Quantifies how much misspecification of the population correlation
#' distorts an AIPE-based sample size plan for the two-one-sided-tests
#' (TOST) confidence interval on the Pearson correlation. On each
#' replication the function simulates \emph{N} bivariate normal pairs
#' with population correlation \code{true_cor}, computes the sample
#' correlation and its Fisher's \eqn{Z} confidence interval via
#' \code{\link{ci_r}}, and summarizes the realized widths and the
#' proportion of replications in which the computed interval falls
#' entirely inside (\code{equivalent}) the specified equivalence
#' bounds.
#'
#' @param true_cor Population correlation (the data generating value).
#'   Defaults to \code{0} (no association, the exact equivalence
#'   case).
#' @param supposed_cor Planning value of the population correlation
#'   passed to \code{\link{ss_aipe_equivalence_r}}; supply this or
#'   \code{specified_N} but not both.
#' @param width Desired full width of the two-sided CI on the
#'   correlation.
#' @param cor_lower,cor_upper Equivalence bounds on the correlation,
#'   as positive magnitudes with the same meaning as in
#'   \code{\link{equivalence_r}}: the region is
#'   \eqn{(-\rho_L, +\rho_U)}, where \eqn{\rho_L} is \code{cor_lower} and
#'   \eqn{\rho_U} is \code{cor_upper}. \code{cor_upper} is required;
#'   \code{cor_lower} defaults to \code{cor_upper} (a symmetric
#'   region). The simulator records whether the realized CI falls
#'   entirely inside the region.
#' @param specified_N Sample size to evaluate.
#' @param conf_level Confidence level (default \code{0.95}). A level above .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param assurance Optional assurance probability.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample correlation, the two confidence limits, the
#'   interval width, whether the interval fell inside the equivalence
#'   region, and two indicators of whether the interval missed
#'   \code{true_cor} below or above): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a
#'   header row is created otherwise, an existing file at that path is
#'   appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized correlation and CI width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of \code{true_cor}, the proportion of intervals
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
#'   true_cor     = 0.0,
#'   supposed_cor = 0.0,
#'   width        = 0.30,
#'   cor_upper    = 0.20,
#'   G = 50, print_iter = FALSE
#' )
#'
#' @keywords design htest
#'
#' @family AIPE sample size planning
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_equivalence_r_sensitivity <- function(true_cor = 0,
                                       supposed_cor = NULL,
                                       width,
                                       cor_lower = NULL,
                                       cor_upper = NULL,
                                       specified_N = NULL,
                                       conf_level = 0.95,
                                       assurance = NULL,
                                       G = 1000, print_iter = FALSE,
                                       filename = NULL,
                                       seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance(assurance)
  if (is.null(supposed_cor) && is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify either 'supposed_cor' or 'specified_N'.",
      arg = c("supposed_cor", "specified_N"))
  if (!is.null(supposed_cor) && !is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify 'supposed_cor' or 'specified_N', but not both.",
      arg = c("supposed_cor", "specified_N"))
  if (!is.numeric(true_cor) || length(true_cor) != 1L || abs(true_cor) >= 1)
    .dmar_refuse_argument(
      "'true_cor' must be a single correlation with |r| < 1.", arg = "true_cor")

  if (is.null(cor_upper))
    .dmar_refuse_argument("'cor_upper' must be specified.", arg = "cor_upper")
  if (!is.numeric(cor_upper) || cor_upper <= 0 || cor_upper >= 1)
    .dmar_refuse_argument("'cor_upper' must be in (0, 1).", arg = "cor_upper")
  if (is.null(cor_lower)) cor_lower <- cor_upper
  if (!is.numeric(cor_lower) || cor_lower <= 0 || cor_lower >= 1)
    .dmar_refuse_argument("'cor_lower' must be in (0, 1).", arg = "cor_lower")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  # ss_aipe_equivalence_r() plans the width of a (1 - 2 * alpha) confidence
  # interval, and the simulation below evaluates ci_r() at conf_level, so
  # the two match only when alpha = (1 - conf_level) / 2 (the same rule
  # ss_aipe_equivalence_smd_sensitivity() follows).
  alpha <- (1 - conf_level) / 2
  if (!is.null(supposed_cor)) {
    plan <- ss_aipe_equivalence_r(true_cor = supposed_cor, width = width,
                           alpha_level = alpha, assurance = assurance)
    N <- plan$value[plan$term == "necessary_N"]
  } else {
    N <- specified_N
  }
  N <- as.integer(N)
  if (N < 4) .dmar_refuse_argument(c(
    "Resolved sample size is < 4; the Fisher's Z interval ",
    "is based on the variance 1/(N - 3)."), arg = NULL)

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
    y <- true_cor * x + sqrt(1 - true_cor^2) * stats::rnorm(N)
    r <- stats::cor(x, y)

    ci <- ci_r(r = r, N = N, conf_level = conf_level)
    lo <- ci$value[ci$term == "ci_lower"]
    hi <- ci$value[ci$term == "ci_upper"]
    r_hat[g]      <- r
    ci_lo[g]      <- lo
    ci_hi[g]      <- hi
    ci_width[g]   <- hi - lo
    equivalent[g] <- (lo > -cor_lower) && (hi < cor_upper)
    tI_lower[g]   <- true_cor < lo
    tI_upper[g]   <- true_cor > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(r = r_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          equivalent = equivalent,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_r", "median_r", "sd_r",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width", "pct_equivalent",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N",
              "true_cor", "supposed_cor",
              "width", "conf_level",
              "cor_lower", "cor_upper",
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
              true_cor,
              if (is.null(supposed_cor)) NA_real_ else supposed_cor,
              width, conf_level,
              -cor_lower, cor_upper,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(ci_width), conf_level)
}
