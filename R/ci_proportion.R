#' Confidence Interval for a Single Proportion
#'
#' The Wilson (1927) score interval for a binomial proportion, the
#' package's default for proportion inference. Unlike the textbook Wald
#' interval, it stays inside [0, 1] and has positive width when the
#' number of successes is 0 or \emph{N}. For a 95\% interval, its
#' coverage averaged over a uniformly distributed true proportion is .954
#' at \emph{N} = 10 and .953 at \emph{N} = 20, against .769 and .846 for
#' the Wald interval; for some true proportions near 0 or 1 the Wilson
#' coverage falls to about .84. For other intervals for a binomial
#' proportion, see Brown, Cai, and DasGupta (2001). The Wald interval is
#' available for instruction and comparison.
#'
#' @param successes Number of successes, a single non-negative integer.
#' @param N Number of trials, a single positive integer at least
#'   \code{successes}.
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param method \code{"wilson"} (default) or \code{"wald"}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{ci_lower}, \code{proportion}, \code{ci_upper},
#'   \code{successes}, and \code{N}, so the point estimate sits between
#'   its confidence limits.
#'
#' @references
#' Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval estimation
#'   for a binomial proportion. \emph{Statistical Science, 16}(2),
#'   101--133. \doi{10.1214/ss/1009213286}
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#'   statistical inference. \emph{Journal of the American Statistical
#'   Association, 22}(158), 209--212.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{responder_analysis}}, which uses this interval for
#'   each group's responder proportion.
#'
#' @family confidence intervals
#'
#' @keywords htest
#'
#' @examples
#' ci_proportion(successes = 17, N = 50)
#'
#' # The Wilson interval stays inside [0, 1] even at the boundary.
#' ci_proportion(successes = 0, N = 20)
#'
#' @export
#' @importFrom stats qnorm
ci_proportion <- function(successes, N, conf_level = 0.95,
                          method = c("wilson", "wald")) {
  method <- match.arg(method)
  if (!is.numeric(successes) || length(successes) != 1L ||
      !is.finite(successes) || successes < 0 || successes != round(successes)) {
    .dmar_refuse_argument(
      "'successes' must be a single non-negative integer.", arg = "successes")
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 1 ||
      N != round(N) || successes > N) {
    .dmar_refuse_argument(
      "'N' must be a single positive integer at least 'successes'.",
      arg = c("N", "successes"))
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  p <- successes / N
  z <- qnorm(1 - (1 - conf_level) / 2)
  if (method == "wilson") {
    den <- 1 + z^2 / N
    ctr <- (p + z^2 / (2 * N)) / den
    hw  <- z * sqrt(p * (1 - p) / N + z^2 / (4 * N^2)) / den
    lims <- c(ctr - hw, ctr + hw)
    # The Wilson limits lie in [0, 1], and the lower limit is 0 exactly at
    # 0 successes and the upper limit 1 exactly at N, where ctr and hw
    # cancel. In floating point that cancellation left values a few ulps
    # either side (the lower limit at 0 successes ranged from -5.6e-17 to
    # 2.8e-17, and the upper limit at N successes from 1 - 2.2e-16 to
    # 1 + 2.2e-16, over N = 1 to 500 and five confidence levels), so the
    # boundary limits are set exactly and the others kept in [0, 1].
    if (successes == 0) lims[1L] <- 0
    if (successes == N) lims[2L] <- 1
    lims <- pmin(pmax(lims, 0), 1)
  } else {
    hw <- z * sqrt(p * (1 - p) / N)
    lims <- c(p - hw, p + hw)
  }
  .as_dmar_tbl(data.frame(
    term  = c("ci_lower", "proportion", "ci_upper", "successes", "N"),
    value = c(lims[1], p, lims[2], successes, N),
    stringsAsFactors = FALSE
  ), conf_level = conf_level)
}
