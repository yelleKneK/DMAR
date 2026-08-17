#' Generate a Population Covariance Matrix From a One-Factor Confirmatory Factor Model
#'
#' Given a vector of factor loadings (\eqn{\lambda}) and the corresponding
#' vector of unique (error) variances (\eqn{\psi^2}), this function computes
#' the implied population covariance matrix under a single-factor
#' confirmatory factor model:
#' \deqn{\Sigma = \Lambda \Lambda^\top + \Psi.}
#' This builds the model implied covariance for a confirmatory factor model
#' with uncorrelated errors. The unique-variances matrix \eqn{\Psi} is strictly
#' diagonal, so no residual covariances (correlated uniquenesses) are allowed; a
#' model with correlated errors is outside the scope of this function.
#' Because the model is a single common factor, the loadings matrix
#' \eqn{\Lambda} reduces to a column vector \eqn{\lambda} of length
#' \eqn{p} (one per indicator), and the unique-variances matrix
#' \eqn{\Psi} is the diagonal \eqn{\mathrm{diag}(\psi^2)} of a length-\eqn{p}
#' vector. The formals are named in the lowercase vector form
#' (\code{lambda} and \code{psi_squared}) to reflect this; the
#' matrix-form symbols (\eqn{\Lambda}, \eqn{\Psi}) remain in the
#' mathematical exposition above.
#'
#' @param lambda A numeric vector of factor loadings (one per indicator).
#'   Can also be supplied as a single-row or single-column matrix and will
#'   be coerced to a vector.
#' @param psi_squared A numeric vector of unique (error) variances, one
#'   per indicator. Recycled to match the length of \code{lambda} when a
#'   single value is supplied (equal-error-variance model).
#' @param \dots Optional advanced controls. Currently the only recognized
#'   passthrough is \code{tol_det}: the tolerance below which the
#'   determinant of the implied covariance matrix triggers a positive-
#'   definite warning (default \code{1e-05}). The argument is hidden in
#'   \code{\dots} because the default is appropriate for almost every
#'   application; users who pass an unrecognized name through
#'   \code{\dots} (for example, a misspelling) are notified with a
#'   warning so the typo is not silently ignored.
#'
#' @return A list with the single element \code{population_cov}: the
#'   implied population covariance matrix of the manifest indicators
#'   (\eqn{p \times p}, symmetric).
#'
#' @details
#' Under the single-factor common-factor model each indicator score is
#' \eqn{x_i = \lambda_i \xi + \delta_i}, where \eqn{\xi} is the
#' (standardized) latent factor and \eqn{\delta_i} is the indicator-specific
#' residual with variance \eqn{\psi_i^2}. The population covariance among
#' the manifest indicators is therefore
#' \deqn{\Sigma = \Lambda \Lambda^\top + \Psi,}
#' where \eqn{\Lambda} is the \eqn{p \times 1} column of loadings
#' \eqn{(\lambda_1, \ldots, \lambda_p)^\top} and
#' \eqn{\Psi = \mathrm{diag}(\psi_1^2, \ldots, \psi_p^2)}. In code we work
#' with the vectors \code{lambda} and \code{psi_squared} directly. Because
#' \eqn{\Psi} is built with \code{diag()} from a length-\eqn{p} vector, the
#' errors are uncorrelated by construction: there is no way to specify a
#' residual covariance between two indicators. A confirmatory factor model
#' with correlated errors (correlated uniquenesses) requires a more general
#' formulation than this function provides.
#'
#' The \code{cfa} spelling matches the \code{\link{cfa_1}} naming.
#'
#' @examples
#' # Five indicators with equal loadings and equal error variances
#' covmat_from_cfa(lambda = rep(0.7, 5), psi_squared = rep(0.51, 5))
#'
#' # Unequal loadings
#' covmat_from_cfa(lambda      = c(0.5, 0.6, 0.7, 0.8),
#'                 psi_squared = c(0.75, 0.64, 0.51, 0.36))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cfa_1}}, \code{\link{ss_aipe_reliability}}
#'
#' @keywords multivariate
#'
#' @export

covmat_from_cfa <- function(lambda, psi_squared, ...) {
  # Pull advanced controls out of `...`. Recognized passthroughs:
  #   tol_det -- determinant threshold for the positive-definite check.
  # Anything else in `...` is almost certainly a typo; warn so it isn't
  # silently dropped.
  dots <- list(...)
  recognized <- "tol_det"
  unknown    <- setdiff(names(dots), recognized)
  if (length(unknown) > 0L) {
    warning(
      "Unrecognized argument(s) passed via `...`: ",
      paste(shQuote(unknown), collapse = ", "),
      ". Recognized passthrough names are: ",
      paste(shQuote(recognized), collapse = ", "),
      ". The unrecognized arguments will be ignored.",
      call. = FALSE
    )
  }
  tol_det <- if ("tol_det" %in% names(dots)) dots$tol_det else 1e-05

  # Coerce lambda to a numeric vector (handles matrix input).
  lambda <- as.numeric(lambda)
  p      <- length(lambda)
  if (p < 2L) {
    stop("'lambda' must contain at least two factor loadings (one per indicator).",
         call. = FALSE)
  }

  psi_squared <- as.numeric(psi_squared)
  if (length(psi_squared) == 1L) {
    psi_squared <- rep(psi_squared, p)
  }
  if (length(psi_squared) != p) {
    stop("'psi_squared' must be the same length as 'lambda' ",
         "(or length 1 for the equal-error-variance model).",
         call. = FALSE)
  }

  # Sigma = lambda lambda^T + diag(psi_squared) for the single-factor model.
  Sigma <- tcrossprod(lambda) + diag(psi_squared, nrow = p)

  if (det(Sigma) < tol_det) {
    warning(
      "The determinant of the implied covariance matrix is very small ",
      "(< tol_det); the matrix may not be positive-definite.",
      call. = FALSE
    )
  }

  list(population_cov = Sigma)
}

