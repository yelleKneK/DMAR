#' Confidence Interval for the Square Root of the Signal-to-Noise Ratio
#'
#' Computes the exact confidence interval for the square root of the
#' signal-to-noise ratio, the standard deviation of the group means
#' relative to the error standard deviation. On this root scale the
#' quantity is an effect size in standard deviation units, the
#' multi-group analog of the standardized mean difference.
#'
#' @param F_value Observed \emph{F}-value from the analysis of variance.
#'   Use this argument when re-analyzing existing data.
#' @param df_1 Numerator degrees of freedom
#' @param df_2 Denominator degrees of freedom
#' @param N Sample size
#' @param means Numeric vector of population or hypothesized group means.
#'   Supply together with \code{sigma_squared} and \code{n_per_group} as a
#'   design-stage alternative to \code{F_value}: the function then computes
#'   the \emph{F}-value implied by these design parameters and proceeds with
#'   the same noncentral F machinery.
#' @param sigma_squared The within-group variance. Used with \code{means}.
#' @param n_per_group A single per-group sample size, or a vector of per-group
#'   sample sizes the same length as \code{means}. Used with \code{means}.
#' @param conf_level Confidence interval coverage (i.e., 1 - Type I error rate); default is .95
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' The confidence level must be specified in one of following two ways: using confidence interval coverage (\code{conf_level}),
#' or lower and upper confidence limits (\code{alpha_lower} and \code{alpha_upper}).
#'
#' The square root of the signal-to-noise ratio is defined as the standard deviation due to the particular factor over the
#' standard deviation of the error (i.e., the square root of the mean square error). This function uses the confidence
#' interval transformation principle (Steiger, 2004) to transform the confidence limits for the noncentrality parameter to
#' the confidence limits for square root of signal-to-noise ratio. The confidence interval for noncentral \emph{F} parameter
#' can be obtained from function \code{conf_limits_ncf} in DMAR.
#'
#' @return
#' A 2-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"lower_limit"} and \code{"upper_limit"},
#' giving the square roots of the corresponding signal-to-noise-ratio
#' confidence limits.
#'
#' @references
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_snr}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' ## To illustrate the calculation of the confidence interval for noncentral
#' ## F parameter,Bargman (1970) gave an example in which a 5-group ANOVA with
#' ## 11 subjects in each group is conducted and the observed F value is 11.221.
#' ## This example continued to be used in Venables (1975),  Fleishman (1980),
#' ## and Steiger (2004). If one wants to calculate the exact confidence interval
#' ## for square root of the signal-to-noise ratio of that example, this
#' ## function can be used.
#'
#' ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#'
#' ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#'
#' ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#'
#' # Design-stage call with population means + within-group variance + n.
#' # Useful when planning a study, before data are observed: derives the
#' # implied F-value internally and returns the resulting CI on the square
#' # root of the signal-to-noise ratio.
#' ci_srsnr(means = c(94, 91, 92, 83), sigma_squared = 67.375, n_per_group = 6)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_srsnr <- function(F_value = NULL, df_1 = NULL, df_2 = NULL, N = NULL,
                     means = NULL, sigma_squared = NULL, n_per_group = NULL,
                     conf_level = .95,
                     alpha_lower = NULL, alpha_upper = NULL, ...) {

  # Design-stage entry: derive F_value, df_1, df_2, N from means + sigma^2 + n.
  if (!is.null(means)) {
    if (is.null(sigma_squared) || is.null(n_per_group)) {
      stop("When 'means' is supplied, also supply 'sigma_squared' and 'n_per_group'.", call. = FALSE)
    }
    if (!is.null(F_value)) {
      stop("Specify either 'F_value' (analysis-of-data mode) or 'means' (design mode), not both.", call. = FALSE)
    }
    if (length(n_per_group) == 1L) n_per_group <- rep(n_per_group, length(means))
    if (length(n_per_group) != length(means)) {
      stop("'n_per_group' must be a scalar or a vector the same length as 'means'.", call. = FALSE)
    }
    if (sigma_squared <= 0) stop("'sigma_squared' must be positive.", call. = FALSE)
    if (any(n_per_group <= 1)) stop("Each per-group sample size must be at least 2.", call. = FALSE)
    J     <- length(means)
    N     <- sum(n_per_group)
    df_1  <- J - 1L
    df_2  <- N - J
    mu_bar     <- sum(n_per_group * means) / N
    SS_between <- sum(n_per_group * (means - mu_bar)^2)
    F_value    <- (SS_between / df_1) / sigma_squared
  }

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_upper - alpha_lower
  }

  if (!is.null(alpha_lower) && is.null(alpha_upper)) stop("This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').")
  if (is.null(alpha_lower) && !is.null(alpha_upper)) stop("This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').")

  if (is.null(F_value)) stop("Specify either 'F_value' (analysis-of-data mode) or 'means' + 'sigma_squared' + 'n_per_group' (design-stage mode).", call. = FALSE)
  if (is.null(df_1) || is.null(df_2) || is.null(N)) stop("You need to specify 'df_1', 'df_2', and 'N'.")

  if (alpha_lower > .5 || alpha_lower < 0) stop(" 'alpha_lower' must be smaller than .5 and nonnegative.")
  if (conf_level > 1 || conf_level < 0) stop(" 'conf_level' must be larger than 0 and smaller than 1. ")
  if (F_value <= 0) stop(" 'F_value' must be larger than 0. ")
  if (N <= 0 || N <= df_1 + df_2) stop("N must be larger than df_1+df_2")

  Lims <- .conf_limits_ncf_for(
    caller = "ci_srsnr",
    quantity = "the square root of the signal-to-noise ratio",
    F_value = F_value, conf_level = NULL, df_1 = df_1,
    df_2 = df_2, alpha_lower = alpha_lower, alpha_upper = alpha_upper, ...
  )

  # print(paste("The", 1-(alpha_lower + alpha_upper), "confidennce limits for the signal to noise ratio are given as:"))
  term <- c("lower_limit", "upper_limit")
  value <- c(sqrt(max(0, Lims[which(Lims$term == "lower_limit"), 2] / N)), sqrt(Lims[which(Lims$term == "upper_limit"), 2] / N))
  length(value) <- length(term)
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level))
}
