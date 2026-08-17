#' Moments of the Noncentral \emph{F} Distribution
#'
#' Returns the mean, variance, standard deviation, skewness, and excess
#' kurtosis of a noncentral \emph{F} distribution with \code{df_1} numerator and
#' \code{df_2} denominator degrees of freedom and noncentrality parameter
#' \code{ncp}. The noncentral \emph{F} is the reference distribution of the
#' \emph{F} statistic when an effect is present, so its moments describe the
#' sampling behavior of \eqn{R^2}, eta squared, and the omnibus \emph{F} test
#' under the alternative. A central \emph{F} (\code{ncp = 0}) is the special
#' case.
#'
#' @param df_1 Numerator degrees of freedom, a single positive number.
#' @param df_2 Denominator degrees of freedom, a single positive number.
#' @param ncp Noncentrality parameter \eqn{\lambda}, a single non-negative
#'   number. Defaults to 0, the central \emph{F}.
#'
#' @details
#' Writing the noncentral \emph{F} as
#' \eqn{F = (X_1/\nu_1)/(X_2/\nu_2)} with \eqn{X_1 \sim \chi^2_{\nu_1}(\lambda)}
#' a noncentral chi square and \eqn{X_2 \sim \chi^2_{\nu_2}} independent, the
#' raw moments are
#' \deqn{\mathrm{E}[F^k] = \Bigl(\tfrac{\nu_2}{\nu_1}\Bigr)^k
#'       \mathrm{E}[X_1^k]\, \prod_{i=1}^{k}\frac{1}{\nu_2 - 2i},
#'       \qquad \nu_2 > 2k,}
#' where the noncentral chi square moments \eqn{\mathrm{E}[X_1^k]} follow from
#' its cumulants \eqn{\kappa_n = 2^{n-1}(n-1)!\,(\nu_1 + n\lambda)}. The mean
#' exists for \eqn{\nu_2 > 2}, the variance for \eqn{\nu_2 > 4}, the skewness
#' for \eqn{\nu_2 > 6}, and the excess kurtosis for \eqn{\nu_2 > 8}; a moment
#' whose denominator degrees of freedom condition is not met is returned as
#' \code{NA}. The mean reduces to the familiar
#' \eqn{\nu_2(\nu_1 + \lambda)/[\nu_1(\nu_2 - 2)]}, and the variance to
#' \eqn{2(\nu_2/\nu_1)^2[(\nu_1 + \lambda)^2 + (\nu_1 + 2\lambda)(\nu_2 - 2)] /
#' [(\nu_2 - 2)^2(\nu_2 - 4)]}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) in
#'   \code{term} / \code{value} layout with the \code{mean}, \code{variance},
#'   \code{sd}, \code{skewness}, and \code{excess_kurtosis} (any of which may be
#'   \code{NA} when \code{df_2} is too small), followed by the \code{df_1},
#'   \code{df_2}, and \code{ncp} that produced them.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). \emph{Continuous
#'   univariate distributions} (Vol. 2, 2nd ed., Chapter 30). Wiley.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{moments_nct}} for the noncentral \emph{t};
#'   \code{\link{conf_limits_ncf}} for the noncentral \emph{F} confidence
#'   limits used in effect size intervals; \code{\link[stats]{df}} for the
#'   density.
#'
#' @family noncentral distribution moments
#'
#' @keywords distribution
#'
#' @examples
#' # A noncentral F with 3 and 40 df and noncentrality 8.
#' moments_ncf(df_1 = 3, df_2 = 40, ncp = 8)
#'
#' # ncp = 0 is the central F: mean df_2 / (df_2 - 2).
#' moments_ncf(df_1 = 3, df_2 = 40)
#'
#' # The variance is undefined for four or fewer denominator df.
#' moments_ncf(df_1 = 2, df_2 = 4, ncp = 5)
#'
#' @export
moments_ncf <- function(df_1, df_2, ncp = 0) {
  if (!is.numeric(df_1) || length(df_1) != 1L || is.na(df_1) || df_1 <= 0) {
    stop("'df_1' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(df_2) || length(df_2) != 1L || is.na(df_2) || df_2 <= 0) {
    stop("'df_2' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(ncp) || length(ncp) != 1L || is.na(ncp) || ncp < 0) {
    stop("'ncp' must be a single non-negative number.", call. = FALSE)
  }

  lam <- ncp
  # Cumulants of the noncentral chi square numerator:
  # kappa_n = 2^(n-1) (n-1)! (df_1 + n*lam).
  k1 <- df_1 + lam
  k2 <- 2 * (df_1 + 2 * lam)
  k3 <- 8 * (df_1 + 3 * lam)
  k4 <- 48 * (df_1 + 4 * lam)
  # Raw moments of the noncentral chi square, E[X1^j] for j = 1..4.
  ex1 <- c(k1,
           k2 + k1^2,
           k3 + 3 * k2 * k1 + k1^3,
           k4 + 4 * k3 * k1 + 3 * k2^2 + 6 * k2 * k1^2 + k1^4)

  # E[F^k] = (df_2/df_1)^k * E[X1^k] * prod_{i=1}^k 1/(df_2 - 2i), df_2 > 2k.
  raw <- rep(NA_real_, 4L)
  for (k in 1:4) {
    if (df_2 > 2 * k) {
      inv_chi <- prod(1 / (df_2 - 2 * seq_len(k)))
      raw[k]  <- (df_2 / df_1)^k * ex1[k] * inv_chi
    }
  }
  m1 <- raw[1L]; m2 <- raw[2L]; m3 <- raw[3L]; m4 <- raw[4L]

  mean_f   <- m1
  var_f    <- if (df_2 > 4) m2 - m1^2 else NA_real_
  sd_f     <- if (is.na(var_f)) NA_real_ else sqrt(var_f)
  skew_f   <- if (df_2 > 6) (m3 - 3 * m1 * m2 + 2 * m1^3) / var_f^(3 / 2) else NA_real_
  exkurt_f <- if (df_2 > 8) {
    (m4 - 4 * m1 * m3 + 6 * m1^2 * m2 - 3 * m1^4) / var_f^2 - 3
  } else NA_real_

  .as_dmar_tbl(data.frame(
    term  = c("mean", "variance", "sd", "skewness", "excess_kurtosis",
              "df_1", "df_2", "ncp"),
    value = c(mean_f, var_f, sd_f, skew_f, exkurt_f, df_1, df_2, ncp),
    stringsAsFactors = FALSE
  ))
}
