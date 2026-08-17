#' Convert a Standard Normal \emph{z} Value to the Corresponding Value on a Normal Distribution
#'
#' This function maps a value on the standard normal distribution (the \emph{z}-distribution, with mean 0 and variance 1) to the equivalent point on a normal distribution with arbitrary mean and standard deviation, \eqn{N(mean, sd^2)}.
#'
#' @param z A value on the standard normal distribution (with mean 0 and variance 1).
#' @param mean The mean of the target normal distribution.
#' @param sd The standard deviation of the target normal distribution.
#'
#' @return A 1-row \code{data.frame} with columns \code{term} and
#'   \code{value}. The \code{term} is \code{"value_from_z"} and
#'   \code{value} is the point on \eqn{N(mean, sd^2)} that lies at the same
#'   percentile as \code{z} does on the standard normal distribution.
#'
#' @details The conversion is \code{value = mean + z * sd}, which places the
#'   returned value at the same percentile of \eqn{N(mean, sd^2)} that
#'   \code{z} occupies on the standard normal distribution. Equivalently,
#'   \code{value = qnorm(pnorm(z), mean, sd)}. With the defaults
#'   (\code{mean = 0}, \code{sd = 1}) the value is returned unchanged, since the
#'   target distribution is then the standard normal distribution itself.
#'
#' @examples
#' # A z value of 1.96 on the standard normal distribution maps to the
#' # corresponding point on a normal distribution with mean 100 and sd 15.
#' convert_z_normal(z = 1.96, mean = 100, sd = 15)
#'
#' # With the default standard normal target, the value is returned unchanged.
#' convert_z_normal(z = 1.96)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{cv_z}}
#'
#' @family parameterization conversions
#'
#' @export
#' @import stats

convert_z_normal <- function(z, mean = 0, sd = 1) {
  if (!is.numeric(z) || length(z) != 1L || is.na(z)) {
    stop("'z' must be a single number. For a vector of z-scores, ",
         "compute mean + z * sd directly.", call. = FALSE)
  }
  if (!is.numeric(mean) || length(mean) != 1L || is.na(mean)) {
    stop("'mean' must be a single number.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || is.na(sd)) {
    stop("'sd' must be a single number.", call. = FALSE)
  }
  if (sd < 0) {
    stop("'sd' must be nonnegative.", call. = FALSE)
  }
  term <- "value_from_z"
  value <- mean + z * sd
  .as_dmar_tbl(data.frame(term, value))
}
