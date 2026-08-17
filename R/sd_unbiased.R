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
#'   of the unbiased variance \eqn{s^2}).
#' @param N The sample size on which \code{s} is based.
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
#' Implementation note: the factor \eqn{c_N} is computed via \code{lgamma}
#' to avoid the overflow of \code{gamma} above \eqn{N \approx 340}, so the
#' function is numerically stable for arbitrarily large \eqn{N}.
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
    if (!is.null(s) || !is.null(N)) stop("Since 'X' was specified, do not specify 's' or 'N'.", call. = FALSE)
    if (anyNA(X)) stop("Missing data in 'X' is not allowed.", call. = FALSE)
    s <- sd(X)
    N <- length(X)
  }
  # The classical bias-correction factor (Holtzman 1950) is
  #   c_N = sqrt((N - 1) / 2) * gamma((N - 1) / 2) / gamma(N / 2),
  # which overflows for N > ~340 because gamma() itself overflows. lgamma()
  # is finite for all positive arguments so c_N is computed via the log
  # form below and is then numerically stable for arbitrarily large N.
  log_c_N <- 0.5 * log((N - 1) / 2) + lgamma((N - 1) / 2) - lgamma(N / 2)
  .as_dmar_tbl(data.frame(term = "sd", value = s * exp(log_c_N)))
}
