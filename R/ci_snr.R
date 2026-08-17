#' Confidence Interval for the Signal-to-Noise Ratio
#'
#' Computes the exact confidence interval for the signal-to-noise ratio in a
#' fixed effects analysis of variance, the variance due to the factor of
#' interest divided by the error variance, expressing the magnitude of an
#' effect relative to the unexplained variability.
#'
#' @param F_value Observed \emph{F}-value from the analysis of variance
#' @param df_1 Numerator degrees of freedom
#' @param df_2 Denominator degrees of freedom
#' @param N Sample size
#' @param conf_level Confidence interval coverage (i.e., 1 - Type I error rate), default is .95
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' The confidence level must be specified in one of following two ways: using confidence interval
#' coverage (\code{conf_level}), or lower and upper confidence limits (\code{alpha_lower} and \code{alpha_upper}).
#'
#' This function uses the confidence interval transformation principle (Steiger, 2004) to transform
#' the confidence limits for the noncentrality parameter to the confidence limits for the population's
#' signal-to-noise ratio. The confidence interval for noncentral \emph{F} parameter can be obtained
#' from the \code{conf_limits_ncf} function in DMAR, which is used internally within this function.
#'
#' @return
#' A 2-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"lower_limit"} and \code{"upper_limit"},
#' giving the lower and upper confidence limits on the signal-to-noise ratio.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis.  \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note The signal to noise ratio is defined as the variance due to the particular factor over the error variance (i.e., the mean square error).
#'
#' @seealso
#' \code{\link{ci_srsnr}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' ## Bargman (1970) gave an example in which a 5-group ANOVA with 11 subjects in each
#' ## group is conducted and the observed F value is 11.221. This example was
#' ## used in Venables (1975),  Fleishman (1980), and Steiger (2004). If one wants to calculate
#' ## the exact confidence interval for the signal-to-noise ratio of that example, this
#' ## function can be used.
#'
#' ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#'
#' ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#'
#' ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_snr <- function(F_value = NULL, df_1 = NULL, df_2 = NULL, N = NULL, conf_level = .95,
                   alpha_lower = NULL, alpha_upper = NULL, ...) {

  ####################################################################### 33
  # Preliminary information and function set-up.

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_upper - alpha_lower
  }

  if (!is.null(alpha_lower) && is.null(alpha_upper)) stop("This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').")
  if (is.null(alpha_lower) && !is.null(alpha_upper)) stop("This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').")

  if (is.null(df_1) || is.null(df_2) || is.null(N)) stop("You need to specify 'df_1', 'df_2', and 'N'.")
  if (is.null(F_value)) stop("You must specify the observed F-value ('F_value') from the analysis of variance.")

  if (alpha_lower > .5 || alpha_lower < 0) stop(" 'alpha_lower' must be smaller than .5 and nonnegative.")
  if (conf_level > 1 || conf_level < 0) stop(" 'conf_level' must be larger than 0 and smaller than 1. ")
  if (F_value <= 0) stop(" 'F_value' must be larger than 0. ")
  if (N <= 0 || N <= df_1 + df_2) stop("N must be larger than df_1+df_2")
  ##########################################################################

  Lims <- .conf_limits_ncf_for(
    caller = "ci_snr", quantity = "the signal-to-noise ratio",
    F_value = F_value, conf_level = NULL, df_1 = df_1,
    df_2 = df_2, alpha_lower = alpha_lower, alpha_upper = alpha_upper, ...
  )

  term <- c("lower_limit", "upper_limit")
  value <- c(Lims$value[Lims$term == "lower_limit"] / N, Lims$value[Lims$term == "upper_limit"] / N)
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level))
}
