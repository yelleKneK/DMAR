#' Unbiased Estimate of the Population Standard Deviation Under Normality
#'
#' Returns the unbiased estimate of the population standard deviation under
#' normality. Although the sample variance \eqn{s^2} is unbiased for
#' \eqn{\sigma^2} when computed with \eqn{N - 1} in the denominator, its
#' square root \eqn{s} is biased downward for \eqn{\sigma} because the
#' square root function is concave (Jensen's inequality). \code{sd_unbiased}
#' applies the classical Holtzman (1950) correction factor that exactly
#' removes that bias under normality.
#'
#' @param s The usual estimate of the standard deviation (the square root
#'   of the unbiased variance \eqn{s^2}). Must not be negative.
#' @param N The sample size on which \code{s} is based. Must exceed 1; it
#'   need not be a whole number.
#' @param X Optional vector of raw scores from which \code{s} and \code{N}
#'   are inferred (\code{s = sd(X)}, \code{N = length(X)}). Mutually
#'   exclusive with \code{s} and \code{N}.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} value is \code{"sd"} and \code{value} is the unbiased
#' estimate of \eqn{\sigma}.
#'
#' @details
#' The sample variance computed with \eqn{N - 1} in the denominator is
#' unbiased for the population variance \eqn{\sigma^2}, but its square
#' root \eqn{s} is biased downward for \eqn{\sigma}. Under normality, the
#' multiplicative bias is
#' \deqn{E[s] \;=\; \sigma \cdot c_N^{-1},
#'   \qquad
#'   c_N \;=\; \sqrt{(N - 1)/2}\,
#'   \cdot \Gamma((N - 1)/2)\,/\,\Gamma(N/2),}
#' (Holtzman, 1950). Multiplying \eqn{s} by \eqn{c_N} therefore yields an
#' unbiased estimator of \eqn{\sigma}. The correction is non-trivial in
#' small samples: \eqn{c_N} is about 1.064 at \eqn{N = 5}, 1.028 at
#' \eqn{N = 10}, 1.009 at \eqn{N = 30}, and is essentially 1 by
#' \eqn{N = 100} (about 1.003). For most applied work the bias of \eqn{s} is small
#' enough to ignore, but it matters when \eqn{s} feeds into downstream
#' quantities (variance components, standardizers, planning calculations)
#' at small \eqn{N}.
#'
#' Implementation note: \eqn{c_N} is computed as the reciprocal of
#' \eqn{c_4(N) = E[s]/\sigma}, and \eqn{c_4(N)} as the exponential of its
#' logarithm. For \eqn{N \ge 101} that logarithm is its asymptotic
#' expansion, \eqn{-1/[4(N - 1)] + 1/[24(N - 1)^3] - \cdots}; for a
#' smaller \eqn{N} it is the expansion at the first of
#' \eqn{N + 2, N + 4, \ldots} that reaches 101, carried back down by the
#' exact recurrence
#' \eqn{\log c_4(N) = \log c_4(N + 2) - \log(1 + 1/[(N - 1)(N + 1)])/2},
#' whose terms all have one sign. Nothing cancels at any \eqn{N}, so the
#' computation avoids both the overflow of \code{gamma} above
#' \eqn{N \approx 340} and the loss of digits in a difference of two
#' \code{lgamma} values, which grows with \eqn{N}. \eqn{c_N} is accurate
#' to two units in the last place at any \eqn{N}, and it is never below 1.
#'
#' @references
#' Holtzman, W. H. (1950). The unbiased estimate of the population variance and standard deviation.
#' \emph{American Journal of Psychology}, \emph{63}, 615--617.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[stats]{sd}}, \code{\link{cv}}
#'
#' @examples
#' set.seed(113)
#' X <- rnorm(10, 100, 15)
#'
#' # The plug-in estimate sqrt(s^2) is biased downward for sigma.
#' var(X)^.5
#'
#' # Holtzman (1950) bias-corrected estimate, supplying s and N.
#' sd_unbiased(s = var(X)^.5, N = length(X))
#'
#' # Equivalent call from the raw vector.
#' sd_unbiased(X = X)
#'
#' @keywords design htest
#'
#' @export


sd_unbiased <- function(s = NULL, N = NULL, X = NULL) {
  if (!is.null(X)) {
    if (!is.null(s) || !is.null(N)) .dmar_refuse_argument(
      "Since 'X' was specified, do not specify 's' or 'N'.",
      arg = c("X", "s", "N"))
    if (anyNA(X)) .dmar_refuse_argument(
      "Missing data in 'X' is not allowed.", arg = "X")
    if (!is.numeric(X) || length(X) < 2L || any(!is.finite(X)))
      .dmar_refuse_argument(
        "'X' must be a numeric vector of at least 2 finite observations.",
        arg = "X")
    s <- sd(X)
    N <- length(X)
  } else {
    if (is.null(s) || is.null(N))
      .dmar_refuse_argument(
        "Specify 's' and 'N', or the raw scores 'X'.", arg = c("s", "N", "X"))
    # In a vector of more than one study a missing element gives a missing
    # row, as it always has; a single missing value, and an inadmissible
    # value anywhere, refuse the call.
    s_seen <- s[!is.na(s)]
    N_seen <- N[!is.na(N)]
    if (!is.numeric(s) || length(s) < 1L || (length(s) == 1L && anyNA(s)) ||
        any(!is.finite(s_seen)) || any(s_seen < 0))
      .dmar_refuse_argument(
        "'s' must be a finite standard deviation that is not negative.",
        arg = "s")
    if (!is.numeric(N) || length(N) < 1L || (length(N) == 1L && anyNA(N)) ||
        any(N_seen <= 1))
      .dmar_refuse_argument(
        "'N' must be a sample size greater than 1.", arg = "N")
    if (length(s) != length(N) && length(s) != 1L && length(N) != 1L)
      .dmar_refuse_argument(
        "'s' and 'N' must have the same length, or one of them length 1.",
        arg = c("s", "N"))
  }
  # The classical bias-correction factor (Holtzman 1950) is
  #   c_N = sqrt((N - 1) / 2) * gamma((N - 1) / 2) / gamma(N / 2) = 1 / c4(N),
  # the reciprocal of E(s) / sigma. gamma() overflows for N above about 340,
  # and a difference of two lgamma() values loses digits as N grows, so
  # c4(N) comes from .c4() (R/c4_internals.R), which forms log c4(N) from
  # its expansion, carried below N = 101 by an exact recurrence, and keeps
  # c_N to two units in the last place at every N (the difference of
  # lgamma() values returned exactly 1 at N = 1e8).
  .as_dmar_tbl(data.frame(term = "sd", value = s / .c4(N)))
}
