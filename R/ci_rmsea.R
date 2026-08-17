#' Confidence Interval for the Population Root Mean Square Error of Approximation
#'
#' Constructs a confidence interval for the population root mean square
#' error of approximation (RMSEA), a population badness-of-fit index for
#' structural equation models. The interval is obtained by inverting the
#' noncentral chi square distribution of the sample fit function
#' \eqn{T = (N - 1) \hat F_{ML}} under the model implied covariance
#' structure, mapping the resulting noncentrality limits to the RMSEA
#' metric (Steiger & Lind, 1980; Browne & Cudeck, 1993).
#'
#' @param rmsea Observed root mean square error of approximation
#' @param df Degrees of freedom of the model
#' @param N Sample size
#' @param conf_level Desired confidence level (e.g., .90, .95, .99)
#' @param alpha_lower The Type I error rate for the lower tail
#' @param alpha_upper The Type I error rate for the upper tail
#'
#' @details
#' The RMSEA expresses the badness of model fit per degree of freedom on
#' the noncentrality scale. Under the noncentral chi square model for the
#' sample fit statistic, the sample \eqn{T = (N - 1) \hat F_{ML}} has
#' approximate noncentral chi square distribution with \eqn{df}
#' degrees of freedom and noncentrality parameter
#' \eqn{\lambda = (N - 1) df \cdot \mathrm{RMSEA}^2}. The CI on
#' \eqn{\mathrm{RMSEA}^2} is obtained by inverting the noncentral chi
#' square distribution at the requested confidence level
#' (\code{\link{conf_limits_nc_chisq}} does the inversion); the bounds are
#' then mapped back to the RMSEA scale via the square root. When the
#' lower noncentrality limit hits zero (\emph{i.e.}, the data are
#' compatible with a well-fitting model), the lower RMSEA limit is
#' truncated at zero because RMSEA is non-negative by construction.
#'
#' The 90 percent CI (rather than the usual 95 percent) is the
#' conventional reporting choice for RMSEA (Browne & Cudeck, 1993) because
#' the upper limit of the 90 percent CI plays a one-sided role in the
#' test of close fit (\eqn{H_0: \mathrm{RMSEA} \le 0.05}). \code{ci_rmsea}
#' defaults to \code{conf_level = 0.95} in line with the rest of the
#' package; pass \code{conf_level = 0.90} when the close fit test is the
#' intended use.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} values are \code{"lower_limit"} (the lower bound of
#' the confidence interval on the population RMSEA, truncated at zero by
#' definition), \code{"rmsea"} (the observed point estimate), and
#' \code{"upper_limit"} (the upper bound).
#'
#' @examples
#' # 1. A typical 95 percent CI on RMSEA.
#' ci_rmsea(rmsea = .055, df = 40, N = 425, conf_level = .95)
#'
#' # 2. The 90 percent CI is the conventional choice when interpretation
#' #    will follow the Browne and Cudeck (1993) close fit decision rule
#' #    (the test of H_0: RMSEA <= 0.05 vs. the upper CI limit). Here the
#' #    upper limit is 0.052, just above the close fit threshold of 0.05
#' #    that Browne and Cudeck recommend, so close fit is not established
#' #    even though the point estimate sits comfortably below it.
#' ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = .90)
#'
#' # 3. Wider model with smaller N: more uncertainty, wider CI.
#' ci_rmsea(rmsea = .055, df = 10, N = 100, conf_level = .90)
#'
#' @references
#' Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing
#'   model fit. In K. A. Bollen & J. S. Long (Eds.), \emph{Testing
#'   structural equation models} (pp. 136--162). Sage.
#'
#' Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
#'   root mean square error of approximation: Sample size planning for
#'   narrow confidence intervals.
#'   \emph{Multivariate Behavioral Research, 46}, 1--32.
#'   \doi{10.1080/00273171.2011.543027}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Steiger, J. H., & Lind, J. C. (1980). \emph{Statistically-based tests for the number of common
#' factors}. Paper presented at the annual Spring meeting of the Psychometric Society, Iowa City, IA.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_rmsea <- function(rmsea, df, N, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL) {
  if (!is.null(alpha_lower) || !is.null(alpha_upper)) {
    if (is.null(alpha_upper) || is.null(alpha_lower)) stop("Both 'alpha_lower' and 'alpha_upper' must be specified")
    if (sum(alpha_lower, alpha_upper) >= 1) stop("There is a problem with the specified confidence limits (their sum should not be greater than 1).", call. = FALSE)
    if (alpha_lower >= 1 || alpha_lower < 0) stop("'alpha_lower' is not correctly specified.", call. = FALSE)
    if (alpha_upper >= 1 || alpha_upper < 0) stop("'alpha_upper' is not correctly specified.", call. = FALSE)
    conf_level <- NULL
  }

  if (!is.null(conf_level)) {
    # Reached only when the alphas were left NULL: the alpha branch above
    # sets conf_level to NULL whenever the alphas are supplied, so the
    # explicit alphas take precedence over the (default) conf_level.
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (rmsea < 0) stop("Your RMSEA cannot be less than zero, as it is defined as the minimum of 0 or the typical point estimate.", call. = FALSE)
  if (df < 0) stop("Your degrees of freedom (i.e., 'df') cannot be less than zero.", call. = FALSE)
  if (N < 0) stop("Your sample size (i.e., 'N') cannot be less than zero.", call. = FALSE)

  chi_sq_statistic <- rmsea^2 * df * (N - 1) + df
  # print("The Chi Square Statistic is:")
  # print(chi_sq_statistic)

  # The clamp-to-zero warning from conf_limits_nc_chisq fires routinely for
  # well-fitting models where the noncentrality lower limit hits the floor;
  # surfacing it on every RMSEA call would be noise. Suppress it here only.
  chi_sq_conf_limits <- suppressWarnings(conf_limits_nc_chisq(
    chi_square = chi_sq_statistic, conf_level = NULL, df = df,
    alpha_lower = alpha_lower, alpha_upper = alpha_upper
  ))

  Low_Limit <- chi_sq_conf_limits$value[chi_sq_conf_limits$term == "lower_limit"]
  Up_Limit  <- chi_sq_conf_limits$value[chi_sq_conf_limits$term == "upper_limit"]
  if (is.na(Low_Limit)) Low_Limit <- 0
  Low_Limit <- max(0, Low_Limit)

  # Verify probabilities.

  # Prob.Greater.Upper <- pchisq(q=chi_sq_statistic, df=df, ncp=chi_sq_conf_limits[2,2])

  ###########################################################################################

  if (Low_Limit == 0) Low_RMSEA <- 0
  if (Low_Limit > 0) Low_RMSEA <- (Low_Limit / (df * (N - 1)))^.5
  if (Low_RMSEA <= 0) message("Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.")
  if (is.na(Low_RMSEA)) {
    Low_RMSEA <- 0
  }
  Low_RMSEA <- max(0, Low_RMSEA)

  # chi.sq.statistic.Low <- Low_RMSEA^2*df*(N-1)+df
  # if(chi.sq.statistic.Low==0) Prob.Less.Lower <- 0
  # if(chi.sq.statistic.Low>0) Prob.Less.Lower <- 1-pchisq(q=chi.sq.statistic.Low, df=df, ncp=Low_Limit)

  # if(round((Prob.Less.Lower + Prob.Greater.Upper), 3) != round((alpha_lower + alpha_upper), 3))
  # {
  # warning("The computed confidence interval does not have the same coverage as the specified confidence interval.", call. = FALSE)
  # if(alpha_lower!=0 && Prob.Less.Lower==0) warning("This is the case in this situation because the confidence interval was truncated at zero (the RMSEA cannot be #negative), and so as to not have the lower limit be negative (which would be statistically absurd), the lower limit was set to zero.", call. = FALSE)
  # }

  Up_RMSEA <- if (is.infinite(Up_Limit)) Inf else (Up_Limit / (df * (N - 1)))^.5

  output <- data.frame(
    term = c("lower_limit", "rmsea", "upper_limit"),
    value = c(Low_RMSEA, rmsea, Up_RMSEA)
  )
  return(.as_dmar_tbl(output, conf_level = conf_level))
}
