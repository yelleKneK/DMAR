#' Moments of the Noncentral \emph{t} Distribution
#'
#' Returns the mean, variance, standard deviation, skewness, and excess
#' kurtosis of a noncentral \emph{t} distribution with \code{df} degrees of
#' freedom and noncentrality parameter \code{ncp}. These are the closed-form
#' moments surveyed by Owen (1968); they are the engine behind the bias and
#' variance of the standardized mean difference (Cohen's \emph{d}), since
#' \emph{d} is a scaled noncentral \emph{t} variate. A central \emph{t}
#' (\code{ncp = 0}) is the special case with mean 0 and the familiar
#' \eqn{\mathit{df}/(\mathit{df}-2)} variance.
#'
#' @param df Degrees of freedom, a single positive number (need not be a whole
#'   number).
#' @param ncp Noncentrality parameter \eqn{\delta}, a single number (may be
#'   negative, which mirrors the distribution about 0). Defaults to 0, the
#'   central \emph{t}.
#'
#' @details
#' Writing the noncentral \emph{t} as \eqn{T = (Z + \delta)/\sqrt{W/\nu}} with
#' \eqn{Z \sim N(0, 1)} and \eqn{W \sim \chi^2_\nu} independent, the raw moments
#' are
#' \deqn{\mathrm{E}[T^k] = \mathrm{E}[(Z + \delta)^k]\,
#'       \Bigl(\tfrac{\nu}{2}\Bigr)^{k/2}\,
#'       \frac{\Gamma\!\bigl((\nu - k)/2\bigr)}{\Gamma(\nu/2)},
#'       \qquad \nu > k,}
#' computed here on the log scale for stability. The mean exists for
#' \eqn{\nu > 1}, the variance for \eqn{\nu > 2}, the skewness for
#' \eqn{\nu > 3}, and the excess kurtosis for \eqn{\nu > 4}; a moment whose
#' degrees of freedom condition is not met is returned as \code{NA}. The mean
#' is \eqn{\delta\sqrt{\nu/2}\,\Gamma((\nu-1)/2)/\Gamma(\nu/2)}, the
#' \eqn{\delta}-scaled reciprocal of the Hedges (1981) bias-correction factor
#' that \code{\link{expected_smd}} and \code{\link{smd}} use; that is why the
#' standardized mean difference is upward biased.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) in
#'   \code{term} / \code{value} layout with the \code{mean}, \code{variance},
#'   \code{sd}, \code{skewness}, and \code{excess_kurtosis} (any of which may be
#'   \code{NA} when the degrees of freedom are too small), followed by the
#'   \code{df} and \code{ncp} that produced them.
#'
#' @references
#' Owen, D. B. (1968). A survey of properties and applications of the
#'   noncentral t-distribution. \emph{Technometrics, 10}(3), 445--478.
#'   \doi{10.1080/00401706.1968.10490590}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#'   size and related estimators. \emph{Journal of Educational Statistics,
#'   6}(2), 107--128.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{moments_ncf}} for the noncentral \emph{F};
#'   \code{\link{expected_smd}} and \code{\link{var_smd}} for the same moments
#'   specialized to Cohen's \emph{d}; \code{\link[stats]{dt}} for the density.
#'
#' @family noncentral distribution moments
#'
#' @keywords distribution
#'
#' @examples
#' # A noncentral t with 20 df and noncentrality 2.5.
#' moments_nct(df = 20, ncp = 2.5)
#'
#' # ncp = 0 is the central t: mean 0, variance df / (df - 2), no skew.
#' moments_nct(df = 10)
#'
#' # The mean is the noncentrality times the Hedges bias factor's reciprocal,
#' # which is why Cohen's d (a scaled noncentral t) is upward biased.
#' m <- moments_nct(df = 18, ncp = 1.2)
#' m$value[m$term == "mean"]
#'
#' @export
moments_nct <- function(df, ncp = 0) {
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("'df' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(ncp) || length(ncp) != 1L || is.na(ncp)) {
    stop("'ncp' must be a single number.", call. = FALSE)
  }

  d <- ncp
  # Raw moments of N(ncp, 1): E[(Z + ncp)^k] for k = 1, 2, 3, 4.
  norm_moment <- c(d, 1 + d^2, d^3 + 3 * d, d^4 + 6 * d^2 + 3)

  # E[T^k] = E[(Z+ncp)^k] * (df/2)^(k/2) * Gamma((df-k)/2)/Gamma(df/2), df > k,
  # evaluated on the log scale. Undefined (NA) when df <= k.
  raw <- rep(NA_real_, 4L)
  for (k in 1:4) {
    if (df > k) {
      raw[k] <- norm_moment[k] *
        exp((k / 2) * log(df / 2) + lgamma((df - k) / 2) - lgamma(df / 2))
    }
  }
  m1 <- raw[1L]; m2 <- raw[2L]; m3 <- raw[3L]; m4 <- raw[4L]

  mean_t   <- m1
  var_t    <- if (df > 2) m2 - m1^2 else NA_real_
  sd_t     <- if (is.na(var_t)) NA_real_ else sqrt(var_t)
  skew_t   <- if (df > 3) (m3 - 3 * m1 * m2 + 2 * m1^3) / var_t^(3 / 2) else NA_real_
  exkurt_t <- if (df > 4) {
    (m4 - 4 * m1 * m3 + 6 * m1^2 * m2 - 3 * m1^4) / var_t^2 - 3
  } else NA_real_

  .as_dmar_tbl(data.frame(
    term  = c("mean", "variance", "sd", "skewness", "excess_kurtosis",
              "df", "ncp"),
    value = c(mean_t, var_t, sd_t, skew_t, exkurt_t, df, ncp),
    stringsAsFactors = FALSE
  ))
}
