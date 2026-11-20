#' Model Implied Covariance Matrix of a One-Factor Confirmatory Factor Model
#'
#' Computes the covariance matrix of the indicators that a single-factor
#' confirmatory factor model implies, from the factor loadings and the
#' residual variances of the indicators:
#' \deqn{\Sigma = \lambda \lambda^\top + \Theta, \qquad
#'   \Theta = \mathrm{diag}(\theta_1, \ldots, \theta_p).}{Sigma = lambda
#'   lambda' + Theta, Theta = diag(theta_1, ..., theta_p).}
#' The loadings and residual variances may be population values, for a
#' covariance matrix to plan a study or to simulate from, or estimates from
#' a fitted model, for the covariance matrix the fitted model implies. The
#' calculation is the same in both cases, and the result is the model
#' implied covariance matrix: the population covariance matrix of the
#' model when the inputs are population values, and the fitted model's
#' covariance matrix when they are estimates.
#'
#' @param factor_loading A numeric vector of factor loadings, one per
#'   indicator (at least two), on the metric of a factor with variance 1.
#'   A single-row or single-column matrix is accepted and read as a vector.
#'   Population loadings or estimates from a fitted model.
#' @param residual_var A numeric vector of residual variances (unique or
#'   error variances), one per indicator, or a single value used for every
#'   indicator (the equal residual variance model). Population values or
#'   estimates from a fitted model.
#' @param tol_det The value below which the determinant of the implied
#'   covariance matrix triggers a warning that the matrix may not be
#'   positive definite (default \code{1e-05}).
#'
#' @return A list with the single element \code{model_implied_cov}: the
#'   model implied covariance matrix of the indicators, a symmetric
#'   \eqn{p \times p} matrix.
#'
#' @details
#' Under the single-factor model, each indicator score is
#' \eqn{x_i = \lambda_i \xi + \delta_i}, where \eqn{\xi} is the common
#' factor and \eqn{\delta_i} is the residual of indicator \emph{i}, with
#' variance \eqn{\theta_i}, uncorrelated with the factor and with the other
#' residuals. The general form of the model implied covariance matrix,
#' \eqn{\Lambda \Phi \Lambda^\top + \Theta}, reduces to
#' \eqn{\lambda \lambda^\top + \Theta} here: \eqn{\Lambda} is the
#' \eqn{p \times 1} column of loadings \eqn{\lambda =
#' (\lambda_1, \ldots, \lambda_p)^\top} (\code{factor_loading}), the factor
#' variance \eqn{\Phi} is fixed at 1, and \eqn{\Theta} is the diagonal
#' matrix of the residual variances (\code{residual_var}). Because
#' \eqn{\Theta} is diagonal, the residuals are uncorrelated by
#' construction; a model with correlated residuals needs a more general
#' formulation than this function provides.
#'
#' Estimates from a fit that fixed one loading at 1 and estimated the
#' factor variance \eqn{\phi} are put on this metric by multiplying each
#' loading by \eqn{\sqrt{\phi}}; the residual variances are unchanged, and
#' so is the implied matrix. \code{\link{cfa_1}} fixes the factor variance
#' at 1, so its loadings and residual variances can be given as they are.
#'
#' Estimated residual variances at or below 0 (a Heywood case) are used as
#' given; the function does not check their sign. The determinant check
#' warns only when the determinant of the implied matrix is below
#' \code{tol_det}, which a Heywood case need not produce: with loadings .9,
#' .8, and 1.1 and residual variances .2, .3, and -.01, the implied matrix
#' is positive definite and no warning is given. A determinant above
#' \code{tol_det} does not show that the matrix is positive definite
#' either, because an even number of negative eigenvalues gives a positive
#' determinant. Check the residual variances themselves for a Heywood case,
#' and the eigenvalues of the result for positive definiteness.
#'
#' The \code{cfa} spelling matches the \code{\link{cfa_1}} naming.
#'
#' @examples
#' # Five indicators with equal loadings and equal residual variances
#' covmat_from_cfa(factor_loading = rep(0.7, 5), residual_var = rep(0.51, 5))
#'
#' # Unequal loadings
#' covmat_from_cfa(factor_loading = c(0.5, 0.6, 0.7, 0.8),
#'                 residual_var   = c(0.75, 0.64, 0.51, 0.36))
#'
#' # From the estimates of a fitted one-factor model, whose factor variance
#' # cfa_1() fixes at 1: the covariance matrix the fitted model implies for
#' # four verbal tests of holzinger_swineford. The tests have no missing
#' # scores, so the model is fit to their covariance matrix, which gives
#' # the estimates the raw scores give.
#' items <- c("t5_general_information", "t6_paragraph_comprehension",
#'            "t7_sentence", "t9_word_meaning")
#' fit <- cfa_1(S = cov(holzinger_swineford[, items]), N = 301)
#' covmat_from_cfa(
#'   factor_loading = fit$estimate[grepl("=~", fit$syntax)],
#'   residual_var = fit$estimate[fit$syntax %in% paste(items, "~~", items)])
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cfa_1}}, \code{\link{ss_aipe_reliability}}
#'
#' @keywords multivariate
#'
#' @export

covmat_from_cfa <- function(factor_loading, residual_var, tol_det = 1e-05) {
  absent <- c("factor_loading", "residual_var")[
    c(missing(factor_loading), missing(residual_var))]
  if (length(absent) > 0L) {
    .dmar_refuse_argument(c(
      paste0("'", absent, "'", collapse = " and "),
      if (length(absent) == 1L) " is" else " are",
      " missing. covmat_from_cfa() needs the factor loadings as ",
      "'factor_loading' and the residual variances as 'residual_var'."),
      arg = absent)
  }
  if (!is.numeric(tol_det) || length(tol_det) != 1L || is.na(tol_det)) {
    .dmar_refuse_argument(c(
      "'tol_det' must be a single number; ", .dmar_given(tol_det), "."),
      arg = "tol_det", value = tol_det)
  }

  .dmar_check_numbers("covmat_from_cfa", c("factor_loading", "residual_var"))
  # A one-row or one-column matrix of loadings is read as a vector.
  factor_loading <- as.numeric(factor_loading)
  p <- length(factor_loading)
  if (p < 2L) {
    .dmar_refuse_argument(
      "'factor_loading' must contain at least two factor loadings (one per indicator).",
      arg = "factor_loading")
  }
  if (!all(is.finite(factor_loading))) {
    .dmar_refuse_argument(
      "Every value of 'factor_loading' must be a finite number.",
      arg = "factor_loading", value = factor_loading)
  }

  residual_var <- as.numeric(residual_var)
  if (length(residual_var) == 1L) {
    residual_var <- rep(residual_var, p)
  }
  if (length(residual_var) != p) {
    .dmar_refuse_argument(c(
      "'residual_var' must be the same length as 'factor_loading' ",
      "(or length 1 for the equal residual variance model)."),
      arg = c("residual_var", "factor_loading"))
  }
  if (!all(is.finite(residual_var))) {
    .dmar_refuse_argument(
      "Every value of 'residual_var' must be a finite number.",
      arg = "residual_var", value = residual_var)
  }

  # Sigma = lambda lambda^T + diag(theta) for the single-factor model with
  # factor variance 1.
  Sigma <- tcrossprod(factor_loading) + diag(residual_var, nrow = p)

  if (det(Sigma) < tol_det) {
    .dmar_warn(c(
      "The determinant of the implied covariance matrix is very small ",
      "(< tol_det); the matrix may not be positive-definite."),
      "result_qualified")
  }

  list(model_implied_cov = Sigma)
}
