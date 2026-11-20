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
#' @param alpha_upper Type I error for the upper confidence limit.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound.
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' The confidence level must be specified in one of following two ways: using confidence interval
#' coverage (\code{conf_level}), or lower and upper confidence limits (\code{alpha_lower} and \code{alpha_upper}).
#'
#' This function uses the confidence interval transformation principle (Steiger, 2004) to transform
#' the confidence limits for the noncentrality parameter to the confidence limits for the population's
#' signal-to-noise ratio. The confidence interval for noncentral \emph{F} parameter can be obtained
#' from the \code{ci_nc_F} function in DMAR, which is used internally within this function.
#'
#' @return
#' A 2-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"ci_lower"} and \code{"ci_upper"},
#' giving the lower and upper confidence limits on the signal-to-noise ratio.
#'
#' @references
#' Bargmann, R. E. (1970). Interpretation and use of a generalized
#'   discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
#'   Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), \emph{Essays in
#'   probability and statistics} (pp. 35--60). University of North
#'   Carolina Press.
#'
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'   \doi{10.1177/001316448004000309}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis.  \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' Venables, W. (1975). Calculation of confidence intervals for
#'   noncentrality parameters. \emph{Journal of the Royal Statistical
#'   Society, Series B, 37}(3), 406--412.
#'   \doi{10.1111/j.2517-6161.1975.tb01554.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note The signal to noise ratio is defined as the variance due to the particular factor over the error variance (i.e., the mean square error).
#'
#' @seealso
#' \code{\link{ci_srsnr}}, \code{\link{ci_nc_F}}
#'
#' @examples
#' ## Bargmann (1970) gave an example of a one-way ANOVA with an observed
#' ## F value of 11.2213 on 4 and 50 degrees of freedom, which Venables (1975),
#' ## Fleishman (1980), and Steiger (2004) used in turn. Steiger describes it
#' ## as five groups of 11 and reports F as 11.221, the value used here.
#' ## The exact confidence interval for the signal-to-noise ratio of that example:
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
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_snr", c("df_1", "df_2", "N"))

  ####################################################################### 33
  # Preliminary information and function set-up.

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_upper - alpha_lower
  }

  if (!is.null(alpha_lower) && is.null(alpha_upper)) .dmar_refuse_argument(
    "This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').",
    arg = c("alpha_lower", "alpha_upper"))
  if (is.null(alpha_lower) && !is.null(alpha_upper)) .dmar_refuse_argument(
    "This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').",
    arg = c("alpha_upper", "alpha_lower"))

  if (is.null(df_1) || is.null(df_2) || is.null(N)) .dmar_refuse_argument(
    "You need to specify 'df_1', 'df_2', and 'N'.",
    arg = c("df_1", "df_2", "N"))
  if (is.null(F_value)) .dmar_refuse_argument(
    "You must specify the observed F-value ('F_value') from the analysis of variance.",
    arg = "F_value")

  if (alpha_lower > .5 || alpha_lower < 0) .dmar_refuse_argument(
    "'alpha_lower' must be between 0 and .5.",
    arg = "alpha_lower")
  if (conf_level > 1 || conf_level < 0) .dmar_refuse_argument(
    "'conf_level' must be between 0 and 1.",
    arg = "conf_level")
  if (F_value <= 0) .dmar_refuse_argument(
    "'F_value' must be larger than 0.", arg = "F_value")
  if (N <= 0 || N <= df_1 + df_2) .dmar_refuse_argument(
    "'N' must be larger than 'df_1' + 'df_2'.", arg = c("N", "df_1", "df_2"))
  ##########################################################################

  Lims <- .ci_nc_F_for(
    caller = "ci_snr", quantity = "the signal-to-noise ratio",
    F_value = F_value, conf_level = NULL, df_1 = df_1,
    df_2 = df_2, alpha_lower = alpha_lower, alpha_upper = alpha_upper, ...
  )

  term <- c("ci_lower", "ci_upper")
  value <- c(Lims$value[Lims$term == "ci_lower"] / N, Lims$value[Lims$term == "ci_upper"] / N)
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level,
                      tail_areas = .tail_areas(alpha_lower, alpha_upper)))
}
