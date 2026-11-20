#' Confidence Interval for the Standardized Mean
#'
#' Computes the exact confidence interval for the standardized mean,
#' the mean divided by the standard deviation, by inverting the
#' noncentral \emph{t} distribution. The standardized mean is the
#' one-sample analog of the standardized mean difference and shares
#' its noncentral interval theory.
#'
#' @param sm Standardized mean
#' @param mean Mean
#' @param sd Standard deviation
#' @param ncp Noncentral parameter
#' @param N Sample size
#' @param conf_level Confidence interval coverage (i.e., 1 - Type I error rate); default is .95. A level above
#'   .999999998, whose tail areas are below 1e-9, is refused
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound. A nonzero tail
#'   area below 1e-9 is refused (see "Small tail areas" in
#'   \code{\link{ci_nc_t}}).
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' The user must specify the standardized mean in one and only one of the three ways:
#' a) mean and standard deviation (\code{mean} and \code{sd}),
#' b) standardized mean (\code{sm}), and
#' c) noncentral parameter (\code{ncp}).
#' The confidence level must be specified in one of following two ways: using confidence interval coverage (\code{conf_level}),
#' or lower and upper confidence limits (\code{alpha_lower} and \code{alpha_upper}). This function uses the exact confidence
#' interval method based on noncentral \emph{t}-distributions. The confidence interval for noncentral \emph{t}-parameter can
#' be obtained from the \code{ci_nc_t} function in DMAR.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"ci_lower"} (the lower confidence limit
#' on the standardized mean), \code{"std_mean"} (the standardized mean),
#' and \code{"ci_upper"} (the upper confidence limit).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note The standardized mean is the mean divided by the standard deviation.
#'
#' @seealso \code{\link{ci_nc_t}}
#'
#' @examples
#' ci_sm(sm = 2.037905, N = 13, conf_level = .95)
#' ci_sm(mean = 30, sd = 14.721, N = 13, conf_level = .95)
#' ci_sm(ncp = 7.347771, N = 13, conf_level = .95)
#' ci_sm(sm = 2.037905, N = 13, alpha_lower = .05, alpha_upper = 0)
#' ci_sm(mean = 50, sd = 10, N = 25, conf_level = .95)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_sm <- function(sm = NULL, mean = NULL, sd = NULL, ncp = NULL, N = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, ...) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_sm", c("sm", "N"))
  SM <- NULL
  if (!is.null(mean) || !is.null(sd)) {
    if (is.null(mean) || is.null(sd)) .dmar_refuse_argument(
      "Since either 'mean' or 'sd' was specified, they must both be specified.",
      arg = c("mean", "sd"))
    if (!is.null(sm) || !is.null(ncp)) .dmar_refuse_argument(
      "since you specified 'mean' and 'sd', you should not specify 'sm' or 'ncp'.",
      arg = c("mean", "sd", "sm"))


    if (sd < 0) .dmar_refuse_argument(
      "The estimated standard deviation is less than zero, yet it should be positive.",
      arg = NULL)
    SM <- mean / sd
  }

  if (!is.null(sm)) {
    if (!is.null(mean) || !is.null(sd) || !is.null(ncp)) .dmar_refuse_argument(
      "Give the standardized mean in one way only: 'sm', 'ncp', or 'mean' with 'sd'.", arg = c("sm", "ncp", "mean", "sd"))
    SM <- sm
  }

  if (!is.null(ncp)) {
    if (!is.null(mean) || !is.null(sd) || !is.null(sm)) .dmar_refuse_argument(
      "Give the standardized mean in one way only: 'sm', 'ncp', or 'mean' with 'sd'.", arg = c("sm", "ncp", "mean", "sd"))
  }

  if (is.null(N)) .dmar_refuse_argument(
    "'N', the sample size, must be given.", arg = "N")

  # Resolve conf_level versus explicit tail areas, mirroring ci_smd(). A supplied
  # pair of alphas defines the coverage as 1 - alpha_lower - alpha_upper, so the
  # default conf_level must not also be attached to the result (its footer would
  # otherwise label, say, a 90% interval as 95%); the result records that
  # coverage instead. Reject mixing an explicit conf_level with the alphas.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      .dmar_refuse_argument(
        "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    conf_level <- NULL
  }

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (!is.null(alpha_lower) && is.null(alpha_upper)) .dmar_refuse_argument(
    "This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').",
    arg = c("alpha_lower", "alpha_upper"))
  if (is.null(alpha_lower) && !is.null(alpha_upper)) .dmar_refuse_argument(
    "This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').",
    arg = c("alpha_upper", "alpha_lower"))
  # The floor on a tail area (R/dmar_conditions.R), here so that a level
  # that leaves too small a tail is refused in its own name; ci_nc_t() is
  # given the tails and checks them against the same floor.
  .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = conf_level)

  ###################################################################################################

  if (!is.null(ncp)) {
    SM <- ncp / sqrt(N)
  }

  if (!is.null(SM)) {
    ncp <- SM * sqrt(N)
  }

  Conf_Limits <- ci_nc_t(ncp = ncp, df = (N - 1), conf_level = NULL, alpha_lower = alpha_lower, alpha_upper = alpha_upper, ...)

  ll <- Conf_Limits[which(Conf_Limits$term == "ci_lower"), 2] / sqrt(N)
  ul <- Conf_Limits[which(Conf_Limits$term == "ci_upper"), 2] / sqrt(N)
  term <- c("ci_lower", "std_mean", "ci_upper")
  value <- c(ll, SM, ul)

  # print(paste("The", 1 - (alpha_lower + alpha_upper), "confidence limits for the standardized mean are given as:"))
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level %||% (1 - alpha_lower - alpha_upper),
                      tail_areas = .tail_areas(alpha_lower, alpha_upper)))
}
