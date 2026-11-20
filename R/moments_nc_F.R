#' Moments of the Noncentral \emph{F} Distribution
#'
#' Returns the mean, variance, standard deviation, skewness, and excess
#' kurtosis of a noncentral \emph{F} distribution with \code{df_1} numerator and
#' \code{df_2} denominator degrees of freedom and noncentrality parameter
#' \code{ncp}, in one table. A moment that does not exist for the given
#' \code{df_2} is returned as \code{NA}.
#'
#' The noncentral \emph{F} is the reference distribution of the
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
#' \strong{The distribution.} The noncentral \emph{F} is
#' \eqn{F = (X_1/\nu_1)/(X_2/\nu_2)} with \eqn{X_1 \sim \chi^2_{\nu_1}(\lambda)}
#' a noncentral chi square and \eqn{X_2 \sim \chi^2_{\nu_2}} independent of
#' it. Its raw moments are
#' \deqn{\mathrm{E}[F^k] = \Bigl(\tfrac{\nu_2}{\nu_1}\Bigr)^k
#'       \mathrm{E}[X_1^k]\, \prod_{i=1}^{k}\frac{1}{\nu_2 - 2i},
#'       \qquad \nu_2 > 2k,}
#' where the noncentral chi square moments \eqn{\mathrm{E}[X_1^k]} follow from
#' its cumulants \eqn{\kappa_j = 2^{j-1}(j-1)!\,(\nu_1 + j\lambda)}
#' (Johnson, Kotz, & Balakrishnan, 1995).
#'
#' \strong{What is returned.} The mean is
#' \eqn{\nu_2(\nu_1 + \lambda)/[\nu_1(\nu_2 - 2)]} and the variance is
#' \eqn{2(\nu_2/\nu_1)^2[(\nu_1 + \lambda)^2 + (\nu_1 + 2\lambda)(\nu_2 - 2)] /
#' [(\nu_2 - 2)^2(\nu_2 - 4)]}. The skewness is the third central moment
#' over the variance to the power 3/2, and the excess kurtosis is the fourth
#' central moment over the squared variance, minus 3 (so that a normal
#' distribution has 0). The mean exists for \eqn{\nu_2 > 2}, the variance
#' for \eqn{\nu_2 > 4}, the skewness for \eqn{\nu_2 > 6}, and the excess
#' kurtosis for \eqn{\nu_2 > 8}; a moment whose condition is not met is
#' returned as \code{NA}, and so is the standard deviation when the variance
#' is.
#'
#' \strong{How the moments are computed.} A central moment formed as a
#' difference of raw moments loses digits when the degrees of freedom are
#' large, because the raw moments agree in their leading digits: for
#' \eqn{\nu_1 = \nu_2 = 10^6} and \eqn{\lambda = 0} that form gives an excess
#' kurtosis of \eqn{-5.29 \times 10^{-5}}, where the value is
#' \eqn{6.60006 \times 10^{-5}}. The function therefore works with the
#' two independent factors of \eqn{F}. Write \eqn{F = (\nu_2/\nu_1) X Y} with
#' \eqn{X = X_1} and \eqn{Y = 1/X_2}, and let \eqn{U = X/\mathrm{E}[X] - 1} and
#' \eqn{V = Y/\mathrm{E}[Y] - 1}, so that
#' \eqn{F/\mathrm{E}[F] - 1 = U + V + UV} with \eqn{U} and \eqn{V} independent
#' and of mean 0. The central moments of \eqn{F/\mathrm{E}[F]} are then sums
#' of products of the central moments of \eqn{U} and \eqn{V}, in which every
#' term is nonnegative and nothing cancels. With \eqn{n = \nu_2},
#' \deqn{u_2 = \kappa_2/\kappa_1^2, \quad u_3 = \kappa_3/\kappa_1^3, \quad
#'       c_u = \kappa_4/\kappa_1^4, \quad u_4 = c_u + 3u_2^2,}
#' \deqn{v_2 = \frac{2}{n - 4}, \quad v_3 = \frac{16}{(n - 4)(n - 6)}, \quad
#'       c_v = \frac{48(5n - 22)}{(n - 4)^2(n - 6)(n - 8)}, \quad
#'       v_4 = c_v + 3v_2^2,}
#' the second and third central moments and the fourth cumulant of \eqn{U}
#' and of \eqn{V} (those of \eqn{V} from
#' \eqn{\mathrm{E}[Y^k] = \prod_{i=1}^{k} 1/(n - 2i)}). Then
#' \deqn{m_2 = u_2 + v_2 + u_2 v_2, \qquad
#'       \mathrm{Var}(F) = \mathrm{E}[F]^2\, m_2,}
#' \deqn{m_3 = u_3 + v_3 + u_3 v_3 + 3u_3 v_2 + 3u_2 v_3 + 6u_2 v_2, \qquad
#'       \mathrm{skewness} = m_3 / m_2^{3/2},}
#' and the excess kurtosis is \eqn{e_4 / m_2^2}, where
#' \deqn{e_4 = c_u + c_v + c_u c_v + 3c_u v_2^2 + 3u_2^2 c_v + 6u_2^2 v_2^2
#'       + 4u_4 v_3 + 4u_3 v_4 + 6c_u v_2 + 12u_2^2 v_2 + 6u_2 c_v
#'       + 12u_2 v_2^2 + 12u_3 v_2 + 12u_2 v_3 + 12u_3 v_3}
#' is the fourth central moment of \eqn{U + V + UV} minus \eqn{3m_2^2}, with
#' the terms that cancel removed algebraically. These expressions are the
#' same quantities as the differences of raw moments above, derived for
#' DMAR; they are not taken from the reference. Checked against exact
#' rational arithmetic in 14 cases with degrees of freedom up to
#' \eqn{10^{12}} and noncentrality up to \eqn{10^6}, the four moments are
#' right to a relative error of \eqn{5 \times 10^{-16}}. At ordinary
#' degrees of freedom (200 or fewer) the values agree with those of the
#' difference of raw moments to a relative difference of \eqn{10^{-11}}.
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
#' @seealso \code{\link{moments_nc_t}} for the noncentral \emph{t};
#'   \code{\link{ci_nc_F}} for the noncentral \emph{F} confidence
#'   limits used in effect size intervals; \code{\link[stats]{df}} for the
#'   density.
#'
#' @family noncentral distribution moments
#'
#' @keywords distribution
#'
#' @examples
#' # A noncentral F with 3 and 40 df and noncentrality 8.
#' moments_nc_F(df_1 = 3, df_2 = 40, ncp = 8)
#'
#' # ncp = 0 is the central F: mean df_2 / (df_2 - 2).
#' moments_nc_F(df_1 = 3, df_2 = 40)
#'
#' # The variance is undefined for four or fewer denominator df.
#' moments_nc_F(df_1 = 2, df_2 = 4, ncp = 5)
#'
#' @export
moments_nc_F <- function(df_1, df_2, ncp = 0) {
  if (!is.numeric(df_1) || length(df_1) != 1L || is.na(df_1) || df_1 <= 0) {
    .dmar_refuse_argument(
      "'df_1' must be a single positive number.", arg = "df_1")
  }
  if (!is.numeric(df_2) || length(df_2) != 1L || is.na(df_2) || df_2 <= 0) {
    .dmar_refuse_argument(
      "'df_2' must be a single positive number.", arg = "df_2")
  }
  if (!is.numeric(ncp) || length(ncp) != 1L || is.na(ncp) || ncp < 0) {
    .dmar_refuse_argument(
      "'ncp' must be a single non-negative number.", arg = "ncp")
  }

  lam <- ncp
  n <- df_2
  # Cumulants of the noncentral chi square numerator X:
  # kappa_j = 2^(j-1) (j-1)! (df_1 + j * lam).
  k1 <- df_1 + lam
  k2 <- 2 * (df_1 + 2 * lam)
  k3 <- 8 * (df_1 + 3 * lam)
  k4 <- 48 * (df_1 + 4 * lam)
  # U = X / E(X) - 1: second and third central moments, fourth cumulant
  # (c4u), fourth central moment.
  u2  <- k2 / k1^2
  u3  <- k3 / k1^3
  c4u <- k4 / k1^4
  u4  <- c4u + 3 * u2^2
  # V = Y / E(Y) - 1 with Y the reciprocal of a chi square on n = df_2
  # degrees of freedom, from E(Y^k) = prod_{i = 1..k} 1 / (n - 2 i); each
  # exists for n > 2 k.
  v2  <- if (n > 4) 2 / (n - 4) else NA_real_
  v3  <- if (n > 6) 16 / ((n - 4) * (n - 6)) else NA_real_
  c4v <- if (n > 8) {
    48 * (5 * n - 22) / ((n - 4)^2 * (n - 6) * (n - 8))
  } else NA_real_
  v4  <- c4v + 3 * v2^2

  # F / E(F) - 1 = U + V + U V with U and V independent and of mean 0, so
  # the central moments of F / E(F) are sums of products of those of U and
  # V. Every term below is nonnegative: nothing cancels, where the
  # differences of raw moments used through 1.0.x lost their leading digits
  # at a large df (the excess kurtosis at df_1 = df_2 = 1e6 came out
  # negative).
  mean_f <- if (n > 2) (n / df_1) * k1 / (n - 2) else NA_real_
  m2     <- u2 + v2 + u2 * v2
  var_f  <- mean_f^2 * m2
  sd_f   <- sqrt(var_f)
  m3     <- u3 + v3 + u3 * v3 + 3 * u3 * v2 + 3 * u2 * v3 + 6 * u2 * v2
  skew_f <- m3 / m2^1.5
  # Fourth central moment of U + V + U V minus 3 m2^2, the canceling terms
  # removed algebraically.
  e4 <- c4u + c4v + c4u * c4v + 3 * c4u * v2^2 + 3 * u2^2 * c4v +
    6 * u2^2 * v2^2 + 4 * u4 * v3 + 4 * u3 * v4 + 6 * c4u * v2 +
    12 * u2^2 * v2 + 6 * u2 * c4v + 12 * u2 * v2^2 + 12 * u3 * v2 +
    12 * u2 * v3 + 12 * u3 * v3
  exkurt_f <- e4 / m2^2

  .as_dmar_tbl(data.frame(
    term  = c("mean", "variance", "sd", "skewness", "excess_kurtosis",
              "df_1", "df_2", "ncp"),
    value = c(mean_f, var_f, sd_f, skew_f, exkurt_f, df_1, df_2, ncp),
    stringsAsFactors = FALSE
  ))
}
