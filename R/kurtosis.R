#' Bias-Corrected Sample Excess Kurtosis
#'
#' Computes the sample excess kurtosis of a numeric vector using the
#' bias-corrected (SAS/SPSS Type 2) formula. Excess kurtosis measures
#' tailedness relative to the normal distribution: zero matches a normal,
#' positive values indicate heavier tails (\dQuote{leptokurtic}), negative
#' values indicate lighter tails (\dQuote{platykurtic}).
#'
#' @param x A numeric vector.
#' @param na_rm Logical. If \code{TRUE} (the default), missing values are
#'   removed before computation. If \code{FALSE}, the result is \code{NA}
#'   when \code{x} contains any \code{NA}.
#'
#' @return A single numeric value: the bias-corrected sample excess
#'   kurtosis, or \code{NA_real_} when fewer than four non-missing
#'   observations are available or when the sample standard deviation is
#'   zero.
#'
#' @details
#' The reported value is
#' \deqn{\hat\gamma_2^{(2)} = \frac{n(n+1)}{(n-1)(n-2)(n-3)}\sum_{i=1}^{n}\left(\frac{x_i - \bar{x}}{s}\right)^4 - \frac{3(n-1)^2}{(n-2)(n-3)},}
#' where \eqn{s} is the (divisor-\eqn{n-1}) sample standard deviation.
#' Subtracting the asymptotic correction \eqn{3(n-1)^2/((n-2)(n-3))}
#' centers the statistic at 0 for a normal distribution; that is, \emph{excess}
#' kurtosis is reported (rather than \dQuote{raw} kurtosis, which centers
#' at 3).
#'
#' \strong{Why isn't this in base R?} See the same Details in
#' \code{\link{skewness}}: R Core defers higher-order moment statistics to
#' contributed packages, partly because multiple formulas (biased Type 1,
#' bias-corrected Type 2, Minitab Type 3) coexist. DMAR adopts Type 2,
#' the form most common in psychometric reporting and used internally by
#' \code{\link{descriptives}}.
#'
#' \strong{Diagnostic interpretation.} As a rough rule of thumb,
#' \eqn{|\mathrm{kurtosis}| > 7} is sometimes flagged as indicative of
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
#' # Normal data: excess kurtosis near zero.
#' set.seed(113)
#' kurtosis(rnorm(1000))
#'
#' # Heavy-tailed data: positive excess kurtosis.
#' kurtosis(rt(1000, df = 4))
#'
#' # The classic 1:5 example: bias-corrected excess kurtosis = -1.2.
#' kurtosis(1:5)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{skewness}}, \code{\link{descriptives}}
#'
#' @family descriptive statistics
#'
#' @keywords univar
#'
#' @export
#' @import stats
kurtosis <- function(x, na_rm = TRUE) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector.", call. = FALSE)
  if (isTRUE(na_rm)) {
    x <- x[!is.na(x)]
  } else if (anyNA(x)) {
    return(NA_real_)
  }
  n <- length(x)
  if (n < 4L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  z <- (x - mean(x)) / s
  term1      <- (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))
  correction <- 3 * (n - 1)^2 / ((n - 2) * (n - 3))
  term1 * sum(z^4) - correction
}
