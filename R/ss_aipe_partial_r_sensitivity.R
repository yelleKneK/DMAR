# Sensitivity analysis for AIPE on a partial correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Partial Correlation
#'
#' @description
#' Quantifies how much misspecification of the population partial
#' correlation distorts an AIPE-based sample size plan. The function
#' constructs an \eqn{(J + 1) \times (J + 1)} population covariance
#' matrix whose implied partial correlation between \eqn{Y} and
#' \eqn{X_1} (controlling for \eqn{X_2, \ldots, X_J}) equals
#' \code{true_rho}, then on each replication draws an \emph{n}-row
#' sample from the corresponding multivariate normal distribution and
#' computes the sample partial correlation and its Fisher's \eqn{Z} CI.
#'
#' @param true_rho Population partial correlation between \emph{Y} and
#'   \eqn{X_1} controlling for \eqn{X_2, \ldots, X_J}; must lie in
#'   \eqn{(-1, 1)}.
#' @param estimated_rho Planning value of the partial correlation
#'   passed to \code{\link{ss_aipe_partial_r}}; supply this or
#'   \code{specified_N} but not both.
#' @param J Total number of predictors (so the partial correlation is
#'   between \emph{Y} and one of the \emph{J} predictors, partialing
#'   out the other \eqn{J - 1}). Must be at least 1.
#' @param width Desired full width of the CI on the partial
#'   correlation.
#' @param specified_N Sample size to evaluate (incompatible with
#'   \code{estimated_rho}).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to
#'   \code{\link{ss_aipe_partial_r}}.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print iteration index per replication.
#' @param save Logical. If \code{TRUE} write per-replication results to
#'   \code{filename}.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for the realized partial
#'   correlation, the interval width, the proportion of intervals at
#'   or below \code{width}, tail-specific and overall non-coverage of
#'   \code{true_rho}, and the input echoes, including \code{assurance} (present only when an
#'   assurance was supplied).
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
#' # Reduced replications and a wide target interval keep this fast.
#' set.seed(113)
#' ss_aipe_partial_r_sensitivity(
#'   true_rho = 0.40, estimated_rho = 0.40, J = 3, width = 0.40,
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
ss_aipe_partial_r_sensitivity <- function(true_rho = NULL,
                                          estimated_rho = NULL,
                                          J, width,
                                          specified_N = NULL,
                                          conf_level = 0.95,
                                          assurance = NULL,
                                          G = 1000, print_iter = FALSE,
                                          save = FALSE,
                                          filename = "ss_aipe_partial_r_sensitivity_result.csv") {
  if (is.null(estimated_rho) && is.null(specified_N))
    stop("You must specify either 'estimated_rho' or 'specified_N'.", call. = FALSE)
  if (!is.null(estimated_rho) && !is.null(specified_N))
    stop("You must specify 'estimated_rho' or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_rho) || !is.numeric(true_rho) || abs(true_rho) >= 1)
    stop("'true_rho' must be a single value in (-1, 1).", call. = FALSE)
  if (!is.numeric(J) || length(J) != 1L || J < 1)
    stop("'J' must be a positive integer.", call. = FALSE)

  if (!is.null(estimated_rho)) {
    plan <- ss_aipe_partial_r(rho = estimated_rho, J = J, width = width,
                              conf_level = conf_level, assurance = assurance)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < J + 3)
    stop("Resolved sample size is too small for J predictors.", call. = FALSE)

  # Population covariance matrix: X is J-dim with independent components,
  # Y has partial corr true_rho with X_1 and zero partial corr with X_2..X_J,
  # achieved by setting Y = true_rho * X_1 + sqrt(1 - true_rho^2) * eps.
  Sigma <- matrix(0, nrow = J + 1, ncol = J + 1)
  diag(Sigma) <- 1
  Sigma[1, 2] <- Sigma[2, 1] <- true_rho

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for the multivariate normal sampling step.")

  pr_hat   <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_width <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    dat <- MASS::mvrnorm(n = n, mu = rep(0, J + 1), Sigma = Sigma)
    Y  <- dat[, 1]; X1 <- dat[, 2]
    if (J == 1L) {
      # No covariates to partial out; the partial correlation reduces to
      # the simple correlation cor(Y, X1).
      r_pr <- stats::cor(Y, X1)
    } else {
      Xrest <- dat[, 3:(J + 1), drop = FALSE]
      e_y  <- stats::residuals(stats::lm(Y  ~ Xrest))
      e_x  <- stats::residuals(stats::lm(X1 ~ Xrest))
      r_pr <- stats::cor(e_y, e_x)
    }
    # Fisher's Z CI for the partial correlation with df = n - J.
    z   <- atanh(r_pr)
    se  <- 1 / sqrt(n - J - 3)
    crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    lo  <- tanh(z - crit * se)
    hi  <- tanh(z + crit * se)
    pr_hat[g]   <- r_pr
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_rho < lo
    tI_upper[g] <- true_rho > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(partial_r = pr_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_partial_r", "median_partial_r", "sd_partial_r",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "J",
              "true_rho", "estimated_rho", "width", "conf_level",
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
              n, J, true_rho,
              if (is.null(estimated_rho)) NA_real_ else estimated_rho,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
