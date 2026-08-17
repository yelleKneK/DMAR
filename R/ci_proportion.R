#' Confidence Interval for a Single Proportion
#'
#' The Wilson (1927) score interval for a binomial proportion, the
#' package's default for proportion inference: unlike the textbook Wald
#' interval it cannot escape [0, 1], behaves sensibly at 0 and 1 counts,
#' and holds close to nominal coverage at small \emph{n} (Brown, Cai, &
#' DasGupta, 2001, recommend it for general use). The Wald interval is
#' available for instruction and comparison.
#'
#' @param successes Number of successes, a single non-negative integer.
#' @param n Number of trials, a single positive integer at least
#'   \code{successes}.
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param method \code{"wilson"} (default) or \code{"wald"}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{lower_limit}, \code{proportion}, \code{upper_limit},
#'   \code{successes}, and \code{n}, so the point estimate sits between
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
#' ci_proportion(successes = 17, n = 50)
#'
#' # The Wilson interval stays inside [0, 1] even at the boundary.
#' ci_proportion(successes = 0, n = 20)
#'
#' @export
#' @importFrom stats qnorm
ci_proportion <- function(successes, n, conf_level = 0.95,
                          method = c("wilson", "wald")) {
  method <- match.arg(method)
  if (!is.numeric(successes) || length(successes) != 1L ||
      is.na(successes) || successes < 0 || successes != round(successes)) {
    stop("'successes' must be a single non-negative integer.",
         call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 ||
      n != round(n) || successes > n) {
    stop("'n' must be a single positive integer at least 'successes'.",
         call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  p <- successes / n
  z <- qnorm(1 - (1 - conf_level) / 2)
  if (method == "wilson") {
    den <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / den
    hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
    lims <- c(ctr - hw, ctr + hw)
  } else {
    hw <- z * sqrt(p * (1 - p) / n)
    lims <- c(p - hw, p + hw)
  }
  .as_dmar_tbl(data.frame(
    term  = c("lower_limit", "proportion", "upper_limit", "successes", "n"),
    value = c(lims[1], p, lims[2], successes, n),
    stringsAsFactors = FALSE
  ), conf_level = conf_level)
}
