# Sensitivity analysis for AIPE on a Pearson correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Pearson Correlation
#'
#' @description
#' Quantifies how much misspecification of the population Pearson
#' correlation distorts an AIPE-based sample size plan. On each
#' replication the function draws an \emph{n}-row sample from a
#' bivariate normal distribution with correlation \code{true_rho} and
#' computes the sample correlation and its Fisher's \eqn{Z} CI, the
#' interval \code{\link{ss_aipe_r}} plans for and
#' \code{\link{correlations_test}} reports. Because the back-transformed
#' width is largest at \eqn{\rho = 0} and shrinks as \eqn{|\rho|} grows,
#' a planning value whose magnitude overstates the population
#' correlation yields realized intervals wider than planned, and the
#' summary rows report by how much.
#'
#' @param true_rho Population Pearson correlation; must lie in
#'   \eqn{(-1, 1)}.
#' @param estimated_rho Planning value of the correlation passed to
#'   \code{\link{ss_aipe_r}}; supply this or \code{specified_N} but not
#'   both.
#' @param width Desired full width of the CI on the correlation.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{estimated_rho}).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to
#'   \code{\link{ss_aipe_r}}.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print iteration index per replication.
#' @param save Logical. If \code{TRUE} write per-replication results to
#'   \code{filename}.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for the realized correlation,
#'   the interval width, the proportion of intervals at or below
#'   \code{width}, tail-specific and overall non-coverage of
#'   \code{true_rho}, and the input echoes, including \code{assurance} (present only when an
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
#'   true_rho = 0.30, estimated_rho = 0.30, width = 0.40,
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
ss_aipe_r_sensitivity <- function(true_rho = NULL,
                                  estimated_rho = NULL,
                                  width,
                                  specified_N = NULL,
                                  conf_level = 0.95,
                                  assurance = NULL,
                                  G = 1000, print_iter = FALSE,
                                  save = FALSE,
                                  filename = "ss_aipe_r_sensitivity_result.csv") {
  if (is.null(estimated_rho) && is.null(specified_N))
    stop("You must specify either 'estimated_rho' or 'specified_N'.", call. = FALSE)
  if (!is.null(estimated_rho) && !is.null(specified_N))
    stop("You must specify 'estimated_rho' or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_rho) || !is.numeric(true_rho) || abs(true_rho) >= 1)
    stop("'true_rho' must be a single value in (-1, 1).", call. = FALSE)

  if (!is.null(estimated_rho)) {
    plan <- ss_aipe_r(rho = estimated_rho, width = width,
                      conf_level = conf_level, assurance = assurance)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < 4L)
    stop("Resolved sample size must be at least 4 for the Fisher's Z interval.", call. = FALSE)

  # Population covariance matrix: standard bivariate normal with
  # correlation true_rho between Y and X.
  Sigma <- matrix(c(1, true_rho, true_rho, 1), nrow = 2, ncol = 2)

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for the multivariate normal sampling step.")

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
    tI_lower[g] <- true_rho < lo
    tI_upper[g] <- true_rho > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(r = r_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_r", "median_r", "sd_r",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N",
              "true_rho", "estimated_rho", "width", "conf_level",
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
              n, true_rho,
              if (is.null(estimated_rho)) NA_real_ else estimated_rho,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
