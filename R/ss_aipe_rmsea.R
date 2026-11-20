#' Sample Size Planning for RMSEA in SEM
#'
#' @description
#' Sample size planning for the population root mean square error of approximation (RMSEA) from the accuracy
#' in parameter estimation (AIPE) perspective. The sample size is planned so that the expected width of a
#' confidence interval for the population RMSEA is no larger than desired.
#'
#' @param true_rmsea The population RMSEA planned for
#' @param df Degrees of freedom of the model
#' @param width Desired confidence interval width
#' @param conf_level Desired confidence level of the planned interval. The
#'   default, 0.90, plans the conventional interval for the RMSEA, the one
#'   \code{\link{ci_rmsea}} forms by default (see Details).
#'
#' @details
#' The width planned is that of the interval \code{\link{ci_rmsea}} forms
#' when the sample RMSEA equals \code{true_rmsea}, with the sample size
#' multiplier \eqn{N - 1} of Kelley and Lai (2011). The RMSEA is
#' conventionally reported with a 90 percent interval, because a
#' \eqn{100(1 - 2\alpha)} percent interval decides the one-sided tests of
#' close and not close fit at level \eqn{\alpha} (MacCallum, Browne, &
#' Sugawara, 1996, p. 136), and every RMSEA interval in the package
#' defaults to it, the interval planned here included. Kelley and Lai
#' (2011) planned 95 percent intervals in their examples; a higher level
#' needs a larger sample for the same width. For an RMSEA of .045 on 5
#' degrees of freedom and a width of .060, the plan is 802 at the default
#' and 1068 with \code{conf_level = 0.95}. The printed table states the
#' level beneath it.
#'
#' @return
#' Returns the necessary total sample size in order to achieve the desired degree of accuracy
#' (i.e., the sufficiently narrow confidence interval).
#'
#' @references
#' Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
#'   root mean square error of approximation: Sample size planning for
#'   narrow confidence intervals.
#'   \emph{Multivariate Behavioral Research, 46}, 1--32.
#'   \doi{10.1080/00273171.2011.543027}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_rmsea}}
#'
#' @examples
#' # The expected width of the 90 percent interval, the default, is .05.
#' ss_aipe_rmsea(true_rmsea = .035, df = 50, width = .05)
#'
#' # The same width for a 95 percent interval needs a larger sample.
#' ss_aipe_rmsea(true_rmsea = .035, df = 50, width = .05, conf_level = .95)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_rmsea <- function(true_rmsea, df, width, conf_level = 0.90) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_numbers("ss_aipe_rmsea", "true_rmsea", required = "true_rmsea")
  if (conf_level > 50 && conf_level < 100) {
    conf_level <- conf_level / 100
  }
  if (conf_level < 0.5 || conf_level > 0.9999) {
    .dmar_refuse_argument(
      "The value of 'conf_level' must be between .5 and .9999",
      arg = "conf_level")
  }

  omega <- width

  width_n <- function(true_rmsea, df, omega, conf_level, N) {
    ci_i <- ci_rmsea(
      rmsea = true_rmsea, df = df, N = N,
      conf_level = conf_level
    )
    w_i <- ci_i[3, 2] - ci_i[1, 2]
    return(abs(omega - w_i))
  }

  ci_2000 <- ci_rmsea(rmsea = true_rmsea, df = df, N = 2000, conf_level = conf_level)
  w_2000 <- ci_2000[3, 2] - ci_2000[1, 2]
  small_N <- ifelse(omega > w_2000, TRUE, FALSE)
  med_N <- lg_N <- FALSE

  if (small_N) {
    ss <- optimize(width_n, c(2, 2000), true_rmsea = true_rmsea, df = df, omega = omega, conf_level = conf_level)
    N <- ceiling(ss$minimum)
  }

  if (!small_N) {
    ci_5000 <- ci_rmsea(rmsea = true_rmsea, df = df, N = 5000, conf_level = conf_level)
    w_5000 <- ci_5000[3, 2] - ci_5000[1, 2]
    med_N <- ifelse(omega > w_5000, TRUE, FALSE)
    lg_N <- ifelse(omega < w_5000, TRUE, FALSE)
  }

  if (med_N) {
    ss <- optimize(width_n, c(2000, 5000), true_rmsea = true_rmsea, df = df, omega = omega, conf_level = conf_level)
    N <- ceiling(ss$minimum)
  }

  if (lg_N) {
    ss <- optimize(width_n, c(5000, 100000), true_rmsea = true_rmsea, df = df, omega = omega, conf_level = conf_level)
    N <- ceiling(ss$minimum)
  }

  ci_v <- ci_rmsea(rmsea = true_rmsea, df = df, N = N, conf_level = conf_level)
  w_v <- ci_v[3, 2] - ci_v[1, 2]
  warn_msg <- paste(
    "The sample size calculated is correct but is probably at the function break point.",
    "\n", "It is recommended to use function 'ci_rmsea' to verify this sample size."
  )

  if (!isTRUE(all.equal(w_v, omega, tolerance = .001))) {
    ci_v1 <- ci_rmsea(rmsea = true_rmsea, df = df, N = N - 1, conf_level = conf_level)
    w_v1 <- ci_v1[3, 2] - ci_v1[1, 2]
    ci_v2 <- ci_rmsea(rmsea = true_rmsea, df = df, N = N + 1, conf_level = conf_level)
    w_v2 <- ci_v2[3, 2] - ci_v2[1, 2]
    if (w_v1 < omega || w_v2 > omega) {
      ss <- optimize(width_n, c(2, 20000), true_rmsea = true_rmsea, df = df, omega = omega, conf_level = conf_level)
      N <- ceiling(ss$minimum)
    } # else warning(warn_msg,immediate. = TRUE)
  }

  .as_dmar_tbl(data.frame(term = "necessary_N", value = N), conf_level = conf_level, subclass = "dmar_ss_aipe")
}
