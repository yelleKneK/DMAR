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
#' \strong{Why isn't this in base R?} R Core has historically deferred
#' higher-order moment statistics to contributed packages, in part because
#' three popular formulas exist (biased Type 1, bias-corrected Type 2, and
#' Minitab Type 3) and choosing a default would be opinionated. DMAR
#' adopts Type 2, which is the form most often used in psychometric
#' reporting and the one already used internally by
#' \code{\link{descriptives}}.
#'
#' \strong{Diagnostic interpretation.} As a rough rule of thumb,
#' \eqn{|\mathrm{skewness}| > 2} is sometimes flagged as indicative of
#' departures from normality large enough to threaten normal-theory
#' inference (e.g., maximum likelihood estimation in factor analysis or
#' structural equation modeling).
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
  if (!is.numeric(x)) stop("'x' must be a numeric vector.", call. = FALSE)
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
