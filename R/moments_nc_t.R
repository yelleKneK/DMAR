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
#'   number). May be \code{Inf}.
#' @param ncp Noncentrality parameter \eqn{\delta}, a single number (may be
#'   negative, which mirrors the distribution about 0). Must be finite.
#'   Defaults to 0, the central \emph{t}.
#'
#' @details
#' Writing the noncentral \emph{t} as \eqn{T = (Z + \delta)/\sqrt{W/\nu}} with
#' \eqn{Z \sim N(0, 1)} and \eqn{W \sim \chi^2_\nu} independent, the raw moments
#' are
#' \deqn{\mathrm{E}[T^k] = \mathrm{E}[(Z + \delta)^k]\,
#'       \Bigl(\tfrac{\nu}{2}\Bigr)^{k/2}\,
#'       \frac{\Gamma\!\bigl((\nu - k)/2\bigr)}{\Gamma(\nu/2)},
#'       \qquad \nu > k.}
#' The mean is formed from the logarithm of the chi mean factor that
#' \code{\link{sd_unbiased}} uses. The variance, skewness, and excess
#' kurtosis are formed from the central moments of \eqn{\sqrt{\nu/W}}, so
#' that no raw moments are subtracted; from 100 degrees of freedom on, the
#' two smallest of those are their expansions in \eqn{1/\nu}. Against values
#' computed to several hundred digits, the largest relative error found is
#' about \eqn{10^{-12}} for degrees of freedom from just above each boundary
#' to \eqn{10^{300}} and any finite \code{ncp}. A mean, variance, or
#' standard deviation above the largest number R holds is returned as
#' \code{Inf}. With \code{df = Inf} the result is the normal distribution
#' with mean \code{ncp} and variance 1.
#'
#' The mean exists for
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
#'   \doi{10.1080/00401706.1968.10490594}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#'   size and related estimators. \emph{Journal of Educational Statistics,
#'   6}(2), 107--128.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{moments_nc_F}} for the noncentral \emph{F};
#'   \code{\link{expected_smd}} and \code{\link{var_smd}} for the same moments
#'   specialized to Cohen's \emph{d}; \code{\link[stats]{dt}} for the density.
#'
#' @family noncentral distribution moments
#'
#' @keywords distribution
#'
#' @examples
#' # A noncentral t with 20 df and noncentrality 2.5.
#' moments_nc_t(df = 20, ncp = 2.5)
#'
#' # ncp = 0 is the central t: mean 0, variance df / (df - 2), no skew.
#' moments_nc_t(df = 10)
#'
#' # The mean is the noncentrality times the Hedges bias factor's reciprocal,
#' # which is why Cohen's d (a scaled noncentral t) is upward biased.
#' m <- moments_nc_t(df = 18, ncp = 1.2)
#' m$value[m$term == "mean"]
#'
#' @export
moments_nc_t <- function(df, ncp = 0) {
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    .dmar_refuse_argument("'df' must be a single positive number.", arg = "df")
  }
  if (!is.numeric(ncp) || length(ncp) != 1L || is.na(ncp)) {
    .dmar_refuse_argument("'ncp' must be a single number.", arg = "ncp")
  }

  if (!is.finite(ncp)) .dmar_refuse_argument(
    "'ncp' must be finite.", arg = "ncp")

  m <- .nc_t_moments(df, ncp)

  .as_dmar_tbl(data.frame(
    term  = c("mean", "variance", "sd", "skewness", "excess_kurtosis",
              "df", "ncp"),
    value = c(unname(m), df, ncp),
    stringsAsFactors = FALSE
  ))
}


# .nc_t_moments(df, ncp): the mean, variance, standard deviation, skewness,
# and excess kurtosis of the noncentral t, in that order, for a single
# df > 0 (Inf allowed) and a single finite ncp. A moment that does not
# exist is NA: the mean needs df > 1, the variance and standard deviation
# df > 2, the skewness df > 3, the excess kurtosis df > 4.
#
# Representation. T = (Z + ncp) A with A = sqrt(df / W), Z standard normal
# and W chi square on df degrees of freedom, independent. Write
#   a  = E(A)   = sqrt(df / (df - 1)) / c(df - 1),
#   K2 = E(A^2) = df / (df - 2),
#   vA = Var(A) = K2 {1 - a^2 / K2},
# where c(m) = exp{L(m)} is the chi mean factor of R/c4_internals.R and
# a^2 / K2 = {1 - 1 / (df - 1)} exp{-2 L(df - 1)}. Then
#   E(T)      = ncp a,
#   Var(T)    = K2 + ncp^2 vA,
#   mu3(T)    = ncp^3 mu3(A) + 3 ncp a K2 / (df - 3),
#   kappa4(T) = ncp^4 kappa4(A)
#               + 12 ncp^2 K2 {K2 + vA (df - 4)} / {(df - 4)(df - 3)}
#               + 6 K2^2 / (df - 4),
# with mu3(A) = a {K2 / (df - 3) - 2 vA} and
# kappa4(A) = 2 K2^2 / (df - 4) - 4 a^2 K2 / (df - 3) + 4 K2 vA - 6 vA^2.
#
# Why this form. The raw moments E(T^k) grow like ncp^k while the central
# moments built from them are differences that grow far more slowly, so
# m2 - m1^2 and its higher-order relatives lose digits as df or ncp grows
# (the excess kurtosis at df = 1e4, ncp = 100 came out negative, and the
# variance at df = 1e16 came out negative). Here every quantity is a sum of
# terms of one sign or a central moment of A formed without a difference
# of large numbers: a comes from L(df - 1), which .chi_log_mean_factor()
# gives without cancellation, and 1 - a^2 / K2 from expm1() and log1p().
#
# The standardized ratios r = ncp / sd(T) and q = K2 / Var(T) are formed
# first, and the skewness and excess kurtosis are written in them, so that
# ncp^3 and ncp^4 are never formed (no overflow at a large ncp) and powers
# of 1 / df are multiplied by r before they can underflow. sd(T) itself is
# the scaled root of a sum of two squares, which does not overflow unless
# the standard deviation does.
#
# From df = 100 on, mu3(A) and kappa4(A) are of order df^-2 and df^-3 and
# their closed forms above lose about log10(df) and 2 log10(df) digits, so
# they are taken from their expansions in x = 1 / df,
#   mu3(A)    = a x^2 Q(x),   kappa4(A) = x^3 R(x),
# with Q through x^10 and R through x^9. The coefficients are the exact
# expansion, checked by rational power-series algebra and against a
# 900-digit evaluation; the first omitted terms are below 2e-13 of the
# leading one at df = 100. Below df = 100 the closed forms are used.
#
# Accuracy. Against references computed to several hundred digits over df
# from just above each boundary to 1e300 and abs(ncp) from 0 to 1e6, the
# largest relative errors found were about 3e-16 (mean), 9e-16 (variance),
# 4e-14 (skewness), and 1.3e-12 (excess kurtosis), the worst of them just
# below df = 100. Above df of about 4.5e307 the series of
# .chi_log_mean_factor() returns 0 (4 df overflows), which bounds the
# statement at 1e300.
.nc_t_moments <- function(df, ncp) {
  out <- c(mean = NA_real_, variance = NA_real_, sd = NA_real_,
           skewness = NA_real_, excess_kurtosis = NA_real_)
  if (is.infinite(df)) {
    out[] <- c(ncp, 1, 1, 0, 0)
    return(out)
  }
  if (df <= 1) return(out)

  log_c <- .chi_log_mean_factor(df - 1)
  a <- exp(0.5 * log1p(1 / (df - 1)) - log_c)
  out[["mean"]] <- ncp * a
  if (df <= 2) return(out)

  K2 <- df / (df - 2)
  # 1 - a^2 / K2, in (0, 1).
  w <- -expm1(log1p(-1 / (df - 1)) - 2 * log_c)
  vA <- K2 * w
  sK <- sqrt(K2)
  sA <- sqrt(vA)
  big <- max(sK, abs(ncp) * sA)
  # When abs(ncp) * sA is above the largest double the standard deviation
  # is too, and it is reported as Inf, as an overflowing variance is; the
  # ratio below would be Inf / Inf.
  overflow <- !is.finite(big)
  sd_t <- if (overflow) Inf else
    big * sqrt((sK / big)^2 + (abs(ncp) * sA / big)^2)
  out[["variance"]] <- sd_t^2
  out[["sd"]] <- sd_t
  if (df <= 3) return(out)

  # ncp / sd and (sK / sd)^2, at their limits when the standard deviation
  # overflows: there sd is abs(ncp) * sA to every digit.
  r <- if (overflow) sign(ncp) / sA else ncp / sd_t
  q <- if (overflow) 0 else (sK / sd_t)^2
  x <- 1 / df
  horner <- function(coef, x) {
    ans <- 0
    for (v in rev(coef)) ans <- ans * x + v
    ans
  }
  series <- df >= 100
  # r^3 mu3(A).
  skew_A <- if (series) {
    a * (r * x)^2 * r *
      horner(c(5 / 4, 69 / 8, 2555 / 64, 19845 / 128, 279685 / 512,
               1857597 / 1024, 95118035 / 16384, 594990465 / 32768,
               7338401455 / 131072, 44838263547 / 262144,
               1088860758895 / 2097152), x)
  } else {
    r^3 * a * (K2 / (df - 3) - 2 * vA)
  }
  out[["skewness"]] <- skew_A + 3 * r * a * q / (df - 3)
  if (df <= 4) return(out)

  # r^4 kappa4(A).
  kurt_A <- if (series) {
    (r^2 * x)^2 * x *
      horner(c(6, 1179 / 16, 8859 / 16, 422775 / 128, 69003 / 4,
               339404283 / 4096, 1542657093 / 4096, 54032237145 / 32768,
               115270685427 / 16384, 15457890206331 / 524288), x)
  } else {
    r^4 * (2 * K2^2 / (df - 4) - 4 * a^2 * K2 / (df - 3) + 4 * K2 * vA -
             6 * vA^2)
  }
  out[["excess_kurtosis"]] <- kurt_A +
    12 * (r^2 / (df - 4)) * q * ((K2 + vA * (df - 4)) / (df - 3)) +
    6 * q^2 / (df - 4)
  out
}
