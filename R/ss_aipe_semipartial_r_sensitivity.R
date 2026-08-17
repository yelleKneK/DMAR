# Sensitivity analysis for AIPE on a semipartial correlation.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Semipartial Correlation
#'
#' @description
#' Quantifies how much misspecification of the population semipartial
#' correlation distorts an AIPE-based sample size plan. The function
#' constructs a population covariance matrix whose implied semipartial
#' correlation between \emph{Y} and \eqn{X_1} (partialing \eqn{X_2,
#' \ldots, X_J} out of \eqn{X_1} only, not out of \emph{Y}) equals
#' \code{true_r_sp}, then on each replication draws an \emph{n}-row
#' sample, computes the sample semipartial correlation, and forms a
#' Fisher's \eqn{Z}-style CI scaled by the standardized regression
#' coefficient.
#'
#' @param true_r_sp Population semipartial correlation; must lie in
#'   \eqn{(-1, 1)}.
#' @param estimated_r_sp Planning value passed to
#'   \code{\link{ss_aipe_semipartial_r}}; supply this or
#'   \code{specified_N} but not both.
#' @param J Total number of predictors. Must be at least 1.
#' @param width Desired full width of the CI on the semipartial
#'   correlation.
#' @param specified_N Sample size to evaluate.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for the realized
#'   semipartial correlation, the interval width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of \code{true_r_sp}, and the input echoes, including \code{assurance} (present only when an
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
#' @seealso \code{\link{ss_aipe_semipartial_r}}, \code{\link{ss_aipe_partial_r_sensitivity}}
#'
#' @examples
#' set.seed(113)
#' ss_aipe_semipartial_r_sensitivity(
#'   true_r_sp = 0.30, estimated_r_sp = 0.30, J = 3, width = 0.20,
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
ss_aipe_semipartial_r_sensitivity <- function(true_r_sp = NULL,
                                              estimated_r_sp = NULL,
                                              J, width,
                                              specified_N = NULL,
                                              conf_level = 0.95,
                                              assurance = NULL,
                                              G = 1000, print_iter = FALSE,
                                              save = FALSE,
                                              filename = "ss_aipe_semipartial_r_sensitivity_result.csv") {
  if (is.null(estimated_r_sp) && is.null(specified_N))
    stop("You must specify either 'estimated_r_sp' or 'specified_N'.", call. = FALSE)
  if (!is.null(estimated_r_sp) && !is.null(specified_N))
    stop("You must specify 'estimated_r_sp' or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_r_sp) || !is.numeric(true_r_sp) || abs(true_r_sp) >= 1)
    stop("'true_r_sp' must be a single value in (-1, 1).", call. = FALSE)
  if (!is.numeric(J) || J < 1)
    stop("'J' must be a positive integer.", call. = FALSE)

  if (!is.null(estimated_r_sp)) {
    plan <- ss_aipe_semipartial_r(r_sp = estimated_r_sp, J = J, width = width,
                                  conf_level = conf_level, assurance = assurance)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    n <- specified_N
  }
  n <- as.integer(n)
  if (n < J + 3)
    stop("Resolved sample size is too small for J predictors.", call. = FALSE)

  # Population covariance: same construction as for partial correlation.
  # The semipartial differs from the partial only in scaling, so the
  # data generating mechanism is identical; we report the semipartial
  # statistic in the simulator.
  Sigma <- matrix(0, nrow = J + 1, ncol = J + 1)
  diag(Sigma) <- 1
  Sigma[1, 2] <- Sigma[2, 1] <- true_r_sp

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for the multivariate normal sampling step.")

  sp_hat   <- numeric(G)
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
      r_sp <- stats::cor(Y, X1)
    } else {
      Xrest <- dat[, 3:(J + 1), drop = FALSE]
      e_x  <- stats::residuals(stats::lm(X1 ~ Xrest))
      # Semipartial: corr(Y, residualized X_1).
      r_sp <- stats::cor(Y, e_x)
    }
    # Approximate CI via Fisher's Z transform applied on the semipartial
    # scale, the partial-correlation convention carried over by analogy.
    # Use the standard one-sample SE = 1/sqrt(n - J - 2)
    # as a first-order approximation; this matches the conservative form
    # used inside the planner.
    z   <- atanh(r_sp)
    se  <- 1 / sqrt(n - J - 2)
    crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    lo  <- tanh(z - crit * se)
    hi  <- tanh(z + crit * se)
    sp_hat[g]   <- r_sp
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_width[g] <- hi - lo
    tI_lower[g] <- true_r_sp < lo
    tI_upper[g] <- true_r_sp > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(r_sp = sp_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_r_sp", "median_r_sp", "sd_r_sp",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "J",
              "true_r_sp", "estimated_r_sp", "width", "conf_level",
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
              n, J, true_r_sp,
              if (is.null(estimated_r_sp)) NA_real_ else estimated_r_sp,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
