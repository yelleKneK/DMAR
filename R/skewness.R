#' Bias-Corrected Sample Skewness
#'
#' Computes the sample skewness of a numeric vector using the bias-corrected
#' (SAS/SPSS Type 2) formula. Skewness measures asymmetry of the
#' distribution: zero is symmetric, positive values indicate a right-tail
#' heavier than the left, negative values the reverse.
#'
#' @param x A numeric vector.
#' @param na_rm Logical. If \code{TRUE} (the default), missing values are
#'   removed before computation. If \code{FALSE}, the result is \code{NA}
#'   when \code{x} contains any \code{NA}.
#'
#' @return A single numeric value: the bias-corrected sample skewness, or
#'   \code{NA_real_} when fewer than three non-missing observations are
#'   available or when the sample standard deviation is zero.
#'
#' @details
#' The reported value is
#' \deqn{\hat\gamma_1^{(2)} = \frac{n}{(n-1)(n-2)}\sum_{i=1}^{n}\left(\frac{x_i - \bar{x}}{s}\right)^3,}
#' where \eqn{s} is the (divisor-\eqn{n-1}) sample standard deviation. This
#' is sometimes called the \dQuote{Type 2} or SAS/SPSS-default form; it is
#' approximately unbiased under normality.
#'
#' \strong{Three formulas.} Three measures of sample skewness are in use.
#' With \eqn{m_r = \sum_{i=1}^{n} (x_i - \bar{x})^r / n}, they are
#' \eqn{g_1 = m_3 / m_2^{3/2}} (Type 1),
#' \eqn{G_1 = g_1 \sqrt{n(n-1)} / (n-2)} (Type 2, the value above), and
#' \eqn{b_1 = g_1 ((n-1)/n)^{3/2}} (Type 3). The three differ by factors
#' that depend only on \eqn{n} and approach 1 as \eqn{n} grows.
#' \code{skewness()} reports \eqn{G_1}, which is also the skewness that
#' \code{\link{descriptives}} reports. \code{\link{kurtosis}} describes the
#' three corresponding measures of kurtosis.
#'
#' @references
#' Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample
#'   skewness and kurtosis. \emph{The Statistician, 47}(1), 183--189.
#'   \doi{10.1111/1467-9884.00122}
#'
#' @examples
#' # Symmetric data: skewness near zero.
#' set.seed(113)
#' skewness(rnorm(1000))
#'
#' # Right-skewed data: positive value.
#' skewness(rexp(1000, rate = 1))
#'
#' # The classic 1:5 example: exactly symmetric (returns 0).
#' skewness(1:5)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{kurtosis}}, \code{\link{descriptives}}
#'
#' @family descriptive statistics
#'
#' @keywords univar
#'
#' @export
#' @import stats
skewness <- function(x, na_rm = TRUE) {
  if (!is.numeric(x)) .dmar_refuse_argument(
    "'x' must be a numeric vector.", arg = "x")
  if (isTRUE(na_rm)) {
    x <- x[!is.na(x)]
  } else if (anyNA(x)) {
    return(NA_real_)
  }
  n <- length(x)
  if (n < 3L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  z <- (x - mean(x)) / s
  (n / ((n - 1) * (n - 2))) * sum(z^3)
}
