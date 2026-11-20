# Sensitivity analysis for AIPE on a Pearson correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Pearson Correlation
#'
#' @description
#' Quantifies how much misspecification of the population Pearson
#' correlation distorts an AIPE-based sample size plan. On each
#' replication the function draws an \emph{n}-row sample from a
#' bivariate normal distribution with correlation \code{true_cor} and
#' computes the sample correlation and its Fisher's \eqn{Z} CI, the
#' interval \code{\link{ss_aipe_r}} plans for and
#' \code{\link{correlations_test}} reports. Because the back-transformed
#' width is largest at \eqn{\rho = 0} and shrinks as \eqn{|\rho|} grows,
#' a planning value whose magnitude overstates the population
#' correlation yields realized intervals wider than planned, and the
#' summary rows report by how much.
#'
#' @param true_cor Population Pearson correlation \eqn{\rho}; must lie in
#'   \eqn{(-1, 1)}.
#' @param supposed_cor Planning value of the correlation passed to
#'   \code{\link{ss_aipe_r}}; supply this or \code{specified_N} but not
#'   both.
#' @param width Desired full width of the CI on the correlation.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{supposed_cor}).
#' @param conf_level Confidence level (default \code{0.95}). A level above .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param assurance Optional assurance probability passed to
#'   \code{\link{ss_aipe_r}}.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print iteration index per replication.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample correlation, the two confidence
#'   limits, the interval width, and two indicators of whether the
#'   interval missed \code{true_cor} below or above): nothing is
#'   written when \code{filename} is \code{NULL} (the default), a new
#'   file with a header row is created otherwise, an existing file at
#'   that path is appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with rows for the realized correlation,
#'   the interval width, the proportion of intervals at or below
#'   \code{width}, tail-specific and overall non-coverage of
#'   \code{true_cor}, and the input echoes, including \code{assurance} (present only when an
#'   assurance was supplied).
#'
#' @references
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
#'   estimating Pearson, Kendall and Spearman correlations.
#'   \emph{Psychometrika, 65}(1), 23--28. \doi{10.1007/BF02294183}
#'
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
#' @seealso \code{\link{ss_aipe_r}}, \code{\link{ss_aipe_partial_r_sensitivity}}
#'
#' @examples
#' # Reduced replications and a wide target interval keep this fast.
#' set.seed(113)
#' ss_aipe_r_sensitivity(
#'   true_cor = 0.30, supposed_cor = 0.30, width = 0.40,
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
ss_aipe_r_sensitivity <- function(true_cor = NULL,
                                  supposed_cor = NULL,
                                  width,
                                  specified_N = NULL,
                                  conf_level = 0.95,
                                  assurance = NULL,
                                  G = 1000, print_iter = FALSE,
                                  filename = NULL,
                                  seed = NULL) {
  .dmar_check_conf_level(conf_level)
  # The floor on a tail area (R/dmar_conditions.R), which ci_r() and
  # ss_aipe_r() keep; checked here as well, since the simulated interval
  # below is formed without them when 'specified_N' is given.
  .dmar_check_tail_floor((1 - conf_level) / 2, (1 - conf_level) / 2,
                         conf_level = conf_level)
  .dmar_check_assurance_single(assurance)
  if (is.null(supposed_cor) && is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify either 'supposed_cor' or 'specified_N'.",
      arg = c("supposed_cor", "specified_N"))
  if (!is.null(supposed_cor) && !is.null(specified_N))
    .dmar_refuse_argument(
      "You must specify 'supposed_cor' or 'specified_N', but not both.",
      arg = c("supposed_cor", "specified_N"))
  if (is.null(true_cor) || !is.numeric(true_cor) || abs(true_cor) >= 1)
    .dmar_refuse_argument(
      "'true_cor' must be a single value in (-1, 1).", arg = "true_cor")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  if (!is.null(supposed_cor)) {
    plan <- ss_aipe_r(true_cor = supposed_cor, width = width,
                      conf_level = conf_level, assurance = assurance)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < 4L)
    .dmar_refuse_argument(
      "Resolved sample size must be at least 4 for the Fisher's Z interval.",
      arg = NULL)

  # Population covariance matrix: standard bivariate normal with
  # correlation true_cor between Y and X.
  Sigma <- matrix(c(1, true_cor, true_cor, 1), nrow = 2, ncol = 2)

  if (!requireNamespace("MASS", quietly = TRUE))
    .dmar_refuse(
      "Package 'MASS' is required for the multivariate normal sampling step.",
      "package_required", package = "MASS")

  r_hat    <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_width <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  se   <- 1 / sqrt(n - 3)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    dat <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
    r <- stats::cor(dat[, 1], dat[, 2])
    # Fisher's Z CI for the Pearson correlation: SE(z) = 1 / sqrt(n - 3)
    # (Bonett & Wright, 2000, Equation 2), back-transformed through tanh.
    z  <- atanh(r)
    lo <- tanh(z - crit * se)
    hi <- tanh(z + crit * se)
    r_hat[g]    <- r
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_cor < lo
    tI_upper[g] <- true_cor > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(r = r_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_r", "median_r", "sd_r",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N",
              "true_cor", "supposed_cor", "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(r_hat, na.rm = TRUE),
              stats::median(r_hat, na.rm = TRUE),
              stats::sd(r_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n, true_cor,
              if (is.null(supposed_cor)) NA_real_ else supposed_cor,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(ci_width), conf_level)
}
