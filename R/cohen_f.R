#' Cohen's \emph{f} Effect Size
#'
#' Computes Cohen's \emph{f} = \eqn{\sigma_m / \sigma}, the population standard
#' deviation of means relative to the within-group standard deviation, by any
#' of three equivalent specifications:
#' \enumerate{
#'   \item raw population means and within-group variance,
#'   \item the population proportion of variance accounted for, \eqn{\eta^2},
#'   \item \eqn{\sigma_m} and \eqn{\sigma} directly.
#' }
#' \emph{Cohen's f is a population quantity}; supplied with population
#' parameters it returns the population value, supplied with sample estimates
#' it returns the corresponding sample value.
#'
#' @param mu Numeric vector of population means (one per group). Use together
#'   with \code{sigma_squared}.
#' @param sigma_squared The within-group variance. Use together with \code{mu}.
#' @param n Optional. Per-group sample sizes (a single number for equal group
#'   sizes, or a vector of length \code{length(mu)} for unequal). When
#'   \code{NULL}, equal weighting across groups is used (i.e., the population
#'   variance of \code{mu} is computed with the \eqn{1/k} divisor).
#' @param eta_squared The population proportion of variance accounted for.
#'   Use this argument alone.
#' @param sigma_m The population standard deviation of the means
#'   (\eqn{\sigma_m}). Use together with \code{sigma}.
#' @param sigma The within-group standard deviation (\eqn{\sigma}). Use
#'   together with \code{sigma_m}.
#'
#' @details
#' All three calling modes return the same value when applied to compatible
#' inputs (Cohen 1988, eq. 8.2.1):
#' \itemize{
#'   \item Raw form: \eqn{f = \sqrt{\sum n_j (\mu_j - \bar\mu)^2 / N \cdot 1/\sigma^2}}.
#'   \item From \eqn{\eta^2}: \eqn{f = \sqrt{\eta^2 / (1 - \eta^2)}}.
#'   \item From the variance ratio: \eqn{f = \sigma_m / \sigma}.
#' }
#'
#' @return A 1-row \code{data.frame} with columns \code{term} and \code{value};
#'   \code{term} is \code{"cohen_f"} and \code{value} is the computed value.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' @seealso \code{\link{ci_srsnr}}, \code{\link{ci_snr}},
#'   \code{\link{ss_power_R2}}
#'
#' @examples
#' # (1) From raw means and within-group variance:
#' cohen_f(mu = c(94, 91, 92, 83), sigma_squared = 67.375)
#'
#' # Equal n weights are the default; equivalent with explicit equal n:
#' cohen_f(mu = c(94, 91, 92, 83), sigma_squared = 67.375, n = 6)
#'
#' # Unequal n:
#' cohen_f(mu = c(94, 91, 92, 83), sigma_squared = 67.375, n = c(4, 6, 5, 5))
#'
#' # (2) From eta_squared:
#' cohen_f(eta_squared = 0.10)
#'
#' # (3) From sigma_m and sigma directly:
#' cohen_f(sigma_m = 4, sigma = 8)
#'
#' @export
cohen_f <- function(mu = NULL, sigma_squared = NULL, n = NULL,
                    eta_squared = NULL, sigma_m = NULL, sigma = NULL) {
  raw_mode   <- !is.null(mu)          && !is.null(sigma_squared)
  eta_mode   <- !is.null(eta_squared)
  ratio_mode <- !is.null(sigma_m)     && !is.null(sigma)

  if (sum(c(raw_mode, eta_mode, ratio_mode)) != 1L) {
    stop("Specify exactly one of: (mu + sigma_squared), eta_squared, or (sigma_m + sigma).",
         call. = FALSE)
  }

  if (raw_mode) {
    if (sigma_squared <= 0) stop("'sigma_squared' must be positive.", call. = FALSE)
    if (length(mu) < 2L)    stop("'mu' must contain at least two group means.", call. = FALSE)
    if (is.null(n)) {
      # Equal weighting across groups: population variance of mu, 1/k divisor.
      sigma_m_sq <- mean((mu - mean(mu))^2)
    } else {
      if (length(n) == 1L) n <- rep(n, length(mu))
      if (length(n) != length(mu)) {
        stop("'n' must be either a single number or a vector of the same length as 'mu'.",
             call. = FALSE)
      }
      if (any(n <= 0)) stop("All 'n' values must be positive.", call. = FALSE)
      N <- sum(n)
      mu_bar <- sum(n * mu) / N
      sigma_m_sq <- sum(n * (mu - mu_bar)^2) / N
    }
    f <- sqrt(sigma_m_sq / sigma_squared)
  } else if (eta_mode) {
    if (eta_squared < 0 || eta_squared >= 1) {
      stop("'eta_squared' must be in [0, 1).", call. = FALSE)
    }
    f <- sqrt(eta_squared / (1 - eta_squared))
  } else {
    if (sigma <= 0) stop("'sigma' must be positive.", call. = FALSE)
    if (sigma_m < 0) stop("'sigma_m' must be non-negative.", call. = FALSE)
    f <- sigma_m / sigma
  }

  .as_dmar_tbl(data.frame(term = "cohen_f", value = f))
}
