#' Moments of the Noncentral Chi Square Distribution
#'
#' Returns the mean, variance, standard deviation, skewness, and excess
#' kurtosis of a noncentral chi square distribution with \code{df} degrees of
#' freedom and noncentrality parameter \code{ncp}. The noncentral chi square is
#' the distribution of a sum of squared independent normals with nonzero means
#' (\eqn{\sum (Z_i + \mu_i)^2}, with \eqn{\lambda = \sum \mu_i^2}); it is the
#' building block of the noncentral \emph{F} (whose numerator is a noncentral
#' chi square) and the reference distribution for likelihood ratio and Wald
#' statistics under the alternative. Unlike the noncentral \emph{t} and
#' \emph{F}, every moment exists, so none of the returned values is ever
#' \code{NA}.
#'
#' @param df Degrees of freedom, a single positive number (need not be a whole
#'   number).
#' @param ncp Noncentrality parameter \eqn{\lambda}, a single non-negative
#'   number. Defaults to 0, the central chi square.
#'
#' @details
#' The cumulants of the noncentral chi square are
#' \eqn{\kappa_n = 2^{n-1}(n-1)!\,(\nu + n\lambda)} for \eqn{n \ge 1}, from
#' which the moments follow directly: the mean is \eqn{\kappa_1 = \nu + \lambda},
#' the variance is \eqn{\kappa_2 = 2(\nu + 2\lambda)}, the skewness is
#' \eqn{\kappa_3 / \kappa_2^{3/2} = \sqrt{8}\,(\nu + 3\lambda)/(\nu + 2\lambda)^{3/2}},
#' and the excess kurtosis is
#' \eqn{\kappa_4 / \kappa_2^{2} = 12(\nu + 4\lambda)/(\nu + 2\lambda)^{2}}. At
#' \eqn{\lambda = 0} these reduce to the central chi square values: mean
#' \eqn{\nu}, variance \eqn{2\nu}, skewness \eqn{\sqrt{8/\nu}}, and excess
#' kurtosis \eqn{12/\nu}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) in
#'   \code{term} / \code{value} layout with the \code{mean}, \code{variance},
#'   \code{sd}, \code{skewness}, and \code{excess_kurtosis}, followed by the
#'   \code{df} and \code{ncp} that produced them.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). \emph{Continuous
#'   univariate distributions} (Vol. 2, 2nd ed., Chapter 29). Wiley.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{moments_ncf}} (whose numerator is a noncentral chi
#'   square) and \code{\link{moments_nct}} for the other noncentral moments;
#'   \code{\link{conf_limits_nc_chisq}} for the noncentral chi square
#'   confidence limits; \code{\link[stats]{dchisq}} for the density.
#'
#' @family noncentral distribution moments
#'
#' @keywords distribution
#'
#' @examples
#' # A noncentral chi square with 5 df and noncentrality 3.
#' moments_nc_chisq(df = 5, ncp = 3)
#'
#' # ncp = 0 is the central chi square: mean df, variance 2 * df.
#' moments_nc_chisq(df = 5)
#'
#' # Every moment exists for any positive df, so nothing is ever NA.
#' anyNA(moments_nc_chisq(df = 1, ncp = 10)$value)
#'
#' @export
moments_nc_chisq <- function(df, ncp = 0) {
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("'df' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(ncp) || length(ncp) != 1L || is.na(ncp) || ncp < 0) {
    stop("'ncp' must be a single non-negative number.", call. = FALSE)
  }

  # Cumulants kappa_n = 2^(n-1) (n-1)! (df + n*ncp). Every moment of the
  # noncentral chi square exists, so (unlike the t and F) none is ever NA.
  k2 <- 2 * (df + 2 * ncp)
  k3 <- 8 * (df + 3 * ncp)
  k4 <- 48 * (df + 4 * ncp)

  mean_x   <- df + ncp
  var_x    <- k2
  sd_x     <- sqrt(k2)
  skew_x   <- k3 / k2^(3 / 2)
  exkurt_x <- k4 / k2^2

  .as_dmar_tbl(data.frame(
    term  = c("mean", "variance", "sd", "skewness", "excess_kurtosis",
              "df", "ncp"),
    value = c(mean_x, var_x, sd_x, skew_x, exkurt_x, df, ncp),
    stringsAsFactors = FALSE
  ))
}
