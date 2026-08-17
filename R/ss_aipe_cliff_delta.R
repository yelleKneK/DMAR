# Sample size for AIPE on Cliff's delta.
#' Sample Size for AIPE on Cliff's \eqn{\delta}
#'
#' Determines the sample size needed for the confidence interval on
#' Cliff's (1993) \eqn{\delta} (and equivalently Vargha-Delaney's
#' \eqn{A = (\delta + 1) / 2}) to have a desired width, using the
#' maximum-variance bound on \eqn{\hat\delta} (Feng & Cliff, 2004,
#' Equation 6, p. 324).
#'
#' @param delta Anticipated population Cliff's \eqn{\delta}; numeric
#'   scalar in \eqn{(-1, 1)}.
#' @param width Desired full width of the CI on \eqn{\delta}.
#' @param which_width \code{"Full"} (default), \code{"Lower"}, or
#'   \code{"Upper"}.
#' @param conf_level Desired confidence level. Default \code{0.95}.
#' @param ratio Ratio \eqn{n_1 / n_2} of the two group sample sizes.
#'   Default \code{1} (balanced).
#' @param assurance Optional. Probability that the realized CI is no
#'   wider than \code{width}.
#'
#' @return A \code{data.frame} with rows for the recommended group
#'   sample sizes \eqn{n_1, n_2}, the expected CI width, and the inputs
#'   echoed back.
#'
#' @details
#' \strong{Maximum-variance bound.} The variance of \eqn{\hat\delta} at a
#' given \eqn{\delta} is largest in the bimodal configuration, where it
#' equals \eqn{(1 - \delta^2)/n_b} with \eqn{n_b} the bimodal group's
#' size; for unequal groups the smaller sample size is used
#' conservatively (Feng & Cliff, 2004, Equation 6 and following text,
#' p. 324):
#' \deqn{\mathrm{Var}(\hat\delta) \;\le\;
#'   \frac{(1 - \delta^2)}{\min(n_1, n_2)}.}
#' Setting the half-width of a Wald-style CI \eqn{z_{1-\alpha/2}
#' \sqrt{\mathrm{Var}(\hat\delta)}} equal to the target half-width and
#' solving gives the recommended per-group sample size. The bound is
#' conservative; the realized CI is generally narrower than the target.
#'
#' \strong{Allocation.} The bound is dominated by \eqn{\min(n_1, n_2)},
#' so balanced allocation (\code{ratio = 1}) is approximately optimal
#' under standard conditions; unbalanced allocations require the larger
#' total \emph{N} to achieve the same precision.
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
#'   ordinal questions. \emph{Psychological Bulletin, 114}(3), 494--509.
#'   \doi{10.1037/0033-2909.114.3.494}
#'
#' Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
#'   improved confidence interval. \emph{Journal of Modern Applied
#'   Statistical Methods, 3}(2), 322--332. \doi{10.22237/jmasm/1099267560}
#'
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of
#'   the CL common language effect size statistics of McGraw and Wong.
#'   \emph{Journal of Educational and Behavioral Statistics, 25}(2),
#'   101--132. \doi{10.3102/10769986025002101}
#'
#' @seealso \code{\link{cliff_delta}}, \code{\link{ss_aipe_partial_r}},
#'   \code{\link{ss_aipe_semipartial_r}}
#'
#' @examples
#' # 1. Plan a balanced design so the 95% CI on delta has full width
#' #        <= 0.20 when anticipating delta = 0.30.
#' ss_aipe_cliff_delta(delta = 0.30, width = 0.20)
#'
#' # 2. Unbalanced: twice as many in group 1.
#' ss_aipe_cliff_delta(delta = 0.30, width = 0.20, ratio = 2)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export

ss_aipe_cliff_delta <- function(delta, width,
                                which_width = c("Full", "Lower", "Upper"),
                                conf_level = 0.95,
                                ratio = 1,
                                assurance = NULL) {
  which_width <- match.arg(which_width)
  if (!is.numeric(delta) || length(delta) != 1L || abs(delta) >= 1)
    stop("'delta' must be a single value in (-1, 1).")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.")
  if (!is.numeric(ratio) || length(ratio) != 1L || ratio <= 0)
    stop("'ratio' (n_1 / n_2) must be positive.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  half_width <- if (which_width == "Full") width / 2 else width
  z_alpha    <- stats::qnorm(1 - (1 - conf_level) / 2)

  # The maximum-variance bound (Feng & Cliff, 2004, Eq. 6):
  # Var(delta) <= (1 - delta^2)/min(n_1, n_2).
  # Setting z_alpha sqrt(Var) <= half_width gives min(n_1, n_2) >=
  #   z_alpha^2 (1 - delta^2) / half_width^2.
  n_min <- ceiling(z_alpha^2 * (1 - delta^2) / half_width^2)
  if (ratio >= 1) {
    n_2 <- n_min
    n_1 <- ceiling(n_min * ratio)
  } else {
    n_1 <- n_min
    n_2 <- ceiling(n_min / ratio)
  }

  if (!is.null(assurance)) {
    if (assurance <= 0.5 || assurance >= 1)
      stop("'assurance' must be in (0.5, 1).")
    df <- max(1, n_1 + n_2 - 2)
    inflate <- stats::qchisq(assurance, df = df) / df
    n_1 <- ceiling(n_1 * inflate)
    n_2 <- ceiling(n_2 * inflate)
  }

  expected_w <- 2 * z_alpha * sqrt((1 - delta^2) / min(n_1, n_2))

  out <- data.frame(
    term  = c("n_1", "n_2", "necessary_N", "expected_width",
              "delta", "ratio", "width_target", "conf_level"),
    value = c(n_1, n_2, n_1 + n_2, expected_w,
              delta, ratio, width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
