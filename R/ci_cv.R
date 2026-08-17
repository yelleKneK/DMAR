#' Confidence Interval for the Coefficient of Variation
#'
#' Computes the noncentral \emph{t}-based confidence interval for the
#' population coefficient of variation, the standard deviation relative to the
#' mean, so variability can be reported on a scale that is free of the
#' measurement units.
#'
#' @param cv Coefficient of variation
#' @param mean Sample mean
#' @param sd Sample standard deviation (square root of the unbiased estimate of the variance
#' @param n Sample size
#' @param data Vector of data for which the confidence interval for the coefficient of variation is to be calculated
#' @param conf_level Desired confidence level (1-Type I error rate)
#' @param alpha_lower The proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf_level})
#' @param alpha_upper The proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf_level})
#' @param ... Allows one to potentially include parameter values for inner functions
#'
#' @details Uses the noncentral \emph{t}-distribution to calculate the confidence interval for the population coefficient of variation.
#'
#' @return
#' A 4-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The rows are ordered so the two
#' point estimates sit between the confidence limits:
#' \code{"lower_limit"} (lower confidence limit on the coefficient of
#' variation), \code{"c_of_v"} (the sample coefficient of variation),
#' \code{"c_of_v_unbiased"} (the unbiased estimator), and
#' \code{"upper_limit"} (upper confidence limit). The \code{prob_less}
#' and \code{prob_greater} columns report the achieved tail probabilities
#' of the noncentral t search at the limit values; they are \code{NA} for
#' the point-estimate rows.
#'
#' @references
#' Johnson, N. L., & Welch, B. L. (1940). Applications of the non-central
#'   \emph{t}-distribution. \emph{Biometrika, 31}(3--4), 362--389.
#'   \doi{10.1093/biomet/31.3-4.362}
#'
#' Kelley, K. (2007). Sample size planning for the coefficient of variation from the accuracy in parameter estimation approach. \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3.)
#'
#' McKay, A. T. (1932). Distribution of the coefficient of variation and
#'   the extended \emph{t} distribution. \emph{Journal of the Royal
#'   Statistical Society, 95}(4), 695--698.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv}}
#'
#' @examples
#' set.seed(113)
#' N <- 15
#' X <- rnorm(N, 5, 1)
#' mean.X <- mean(X)
#' sd.X <- var(X)^.5
#'
#' ci_cv(mean = mean.X, sd = sd.X, n = N, alpha_lower = .025,
#' alpha_upper = .025, conf_level = NULL)
#' ci_cv(data = X, conf_level = .95)
#' ci_cv(cv = sd.X / mean.X, n = N, conf_level = .95)
#'
#' @keywords models htest
#'
#' @import stats
#' @family confidence intervals for effect sizes
#'
#' @export

ci_cv <- function(cv = NULL, mean = NULL, sd = NULL, n = NULL, data = NULL,
                  conf_level = 0.95, alpha_lower = NULL, alpha_upper = NULL,
                  ...) {
  if (is.null(conf_level)) {
    if (is.null(alpha_lower) || is.null(alpha_upper)) {
      stop("When 'conf_level' is NULL, both 'alpha_lower' and 'alpha_upper' must be specified.")
    }
    if (alpha_lower >= 1 || alpha_lower < 0) {
      stop("'alpha_lower' is not correctly specified.")
    }
    if (alpha_upper >= 1 || alpha_upper < 0) {
      stop("'alpha_upper' is not correctly specified.")
    }
  }
  if (!is.null(conf_level)) {
    if (!is.null(alpha_lower) || !is.null(alpha_upper)) {
      stop("Since 'conf_level' is specified, 'alpha_lower' and 'alpha_upper' should be 'NULL'.")
    }
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  if (is.null(data) && is.null(cv)) {
    if (is.null(mean)) {
      stop("Either input the whole data set using 'data' or specify the sample mean.")
    }
    if (mean <= 0) {
      stop("The sample mean must be some non-zero positive value (Does taking the absolute value of the mean make sense?).")
    }
    if (is.null(sd)) {
      stop("Either input the whole data set using 'data' or specify the sample standard deviation (using 'n'-1 in the denominator).")
    }
    if (is.null(n)) {
      stop("Either input the whole data set using 'data' or specify the sample size")
    }
    k <- sd / mean
    ncp_estimate <- sqrt(n) / k
    CI_NCP <- conf_limits_nct(
      ncp = ncp_estimate, df = n -
        1, alpha_lower = alpha_upper, alpha_upper = alpha_lower,
      conf_level = NULL
    )
    Low_lim <- CI_NCP[which(CI_NCP$term == "upper_limit"), 2]
    Up_lim <- CI_NCP[which(CI_NCP$term == "lower_limit"), 2]
    Low_Lim <- sqrt(n) / Low_lim
    Up_Lim <- sqrt(n) / Up_lim
    if (Up_lim <= 0) {
      Up_Lim <- Inf
    }

    Result <- data.frame(
      term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
      value = c(Low_Lim, k, k * (1 + 1 / (4 * n)), Up_Lim),
      prob_less = c(alpha_lower, NA, NA, 1 - alpha_upper),
      prob_greater = c(1 - alpha_lower, NA, NA, alpha_upper)
    )
    if (alpha_lower == 0) {
      Result <- data.frame(
        term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
        value = c(-Inf, k, k * (1 + 1 / (4 * n)), Up_Lim),
        prob_less = c(0, NA, NA, 1 - alpha_upper),
        prob_greater = c(1, NA, NA, alpha_upper)
      )
    }
    if (alpha_upper == 0) {
      Result <- data.frame(
        term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
        value = c(Low_Lim, k, k * (1 + 1 / (4 * n)), Inf),
        prob_less = c(alpha_lower, NA, NA, 1),
        prob_greater = c(1 - alpha_lower, NA, NA, 0)
      )
    }
    if (Up_Lim == Inf) {
      Result[4, 3] <- 1
      Result[4, 4] <- 0
    }
    if (round(
      (Result[1, 3] + Result[4, 4]),
      3
    ) != round((alpha_lower + alpha_upper), 3)) {
      warning("The computed confidence interval does not have the same coverage as the specified confidence interval.",
        call. = FALSE
      )
    }
    return(.as_dmar_tbl(Result, conf_level = conf_level))
  }
  if (!is.null(data) && is.null(cv)) {
    if (!is.null(mean)) {
      stop("Since 'data' is specified, do not specify the 'mean'.")
    }
    if (!is.null(sd)) {
      stop("Since 'data' is specified, do not specify the 'sd'.")
    }
    if (!is.null(n)) {
      stop("Since 'data' is specified, do not specify the 'n'.")
    }
    n <- length(data)
    sd_data <- (var(data))^0.5
    mean_data <- mean(data)
    k <- sd_data / mean_data
    ncp_estimate <- sqrt(n) / k
    CI_NCP <- conf_limits_nct(
      ncp = ncp_estimate, df = n -
        1, alpha_lower = alpha_upper, alpha_upper = alpha_lower,
      conf_level = NULL
    )
    Low_lim <- CI_NCP[2, 2]
    Up_lim <- CI_NCP[1, 2]
    Low_Lim <- sqrt(n) / Low_lim
    Up_Lim <- sqrt(n) / Up_lim
    if (Up_lim <= 0) {
      Up_Lim <- Inf
    }
    Result <- data.frame(
      term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
      value = c(Low_Lim, k, k * (1 + 1 / (4 * n)), Up_Lim),
      prob_less = c(alpha_lower, NA, NA, 1 - alpha_upper),
      prob_greater = c(1 - alpha_lower, NA, NA, alpha_upper)
    )
    if (alpha_lower == 0) {
      Result <- data.frame(
        term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
        value = c(-Inf, k, k * (1 + 1 / (4 * n)), Up_Lim),
        prob_less = c(0, NA, NA, 1 - alpha_upper),
        prob_greater = c(1, NA, NA, alpha_upper)
      )
    }
    if (alpha_upper == 0) {
      Result <- data.frame(
        term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
        value = c(Low_Lim, k, k * (1 + 1 / (4 * n)), Inf),
        prob_less = c(alpha_lower, NA, NA, 1),
        prob_greater = c(1 - alpha_lower, NA, NA, 0)
      )
    }
    if (Up_Lim == Inf) {
      Result[4, 3] <- 1
      Result[4, 4] <- 0
    }
    if (round(
      (Result[1, 3] + Result[4, 4]),
      3
    ) != round((alpha_lower + alpha_upper), 3)) {
      warning("The computed confidence interval does not have the same coverage as the specified confidence interval.",
        call. = FALSE
      )
    }
    return(.as_dmar_tbl(Result, conf_level = conf_level))
  }
  if (!is.null(cv)) {
    k <- cv
    if (is.null(n)) {
      stop("Since you specified the coefficient of variation directly ('cv'), you must specify the sample size.")
    }
    if (!is.null(data)) {
      stop("Since you specified the coefficient of variation ('cv') directly, do not include the raw data.")
    }
    if (!is.null(mean)) {
      stop("Since you specified the coefficient of variation ('cv') directly, do not specify the mean ('mean').")
    }
    if (!is.null(sd)) {
      stop("Since you specified the coefficient of variation ('cv') directly, do not specify the standard deviation ('sd').")
    }
    ncp_estimate <- sqrt(n) / k
    CI_NCP <- conf_limits_nct(
      ncp = ncp_estimate, df = n -
        1, alpha_lower = alpha_upper, alpha_upper = alpha_lower,
      conf_level = NULL
    )
    Low_lim <- CI_NCP[2, 2]
    Up_lim <- CI_NCP[1, 2]
    Low_Lim <- sqrt(n) / Low_lim
    Up_Lim <- sqrt(n) / Up_lim
    if (Up_lim <= 0) {
      Up_Lim <- Inf
    }
    Result <- data.frame(
      term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
      value = c(Low_Lim, k, k * (1 + 1 / (4 * n)), Up_Lim),
      prob_less = c(alpha_lower, NA, NA, 1 - alpha_upper),
      prob_greater = c(1 - alpha_lower, NA, NA, alpha_upper)
    )
    if (alpha_lower == 0) {
      Result <- data.frame(
        term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
        value = c(-Inf, k, k * (1 + 1 / (4 * n)), Up_Lim),
        prob_less = c(0, NA, NA, 1 - alpha_upper),
        prob_greater = c(1, NA, NA, alpha_upper)
      )
    }
    if (alpha_upper == 0) {
      Result <- data.frame(
        term = c("lower_limit", "c_of_v", "c_of_v_unbiased", "upper_limit"),
        value = c(Low_Lim, k, k * (1 + 1 / (4 * n)), Inf),
        prob_less = c(alpha_lower, NA, NA, 1),
        prob_greater = c(1 - alpha_lower, NA, NA, 0)
      )
    }
    if (Up_Lim == Inf) {
      Result[4, 3] <- 1
      Result[4, 4] <- 0
    }
    if (round(
      (Result[1, 3] + Result[4, 4]),
      3
    ) != round((alpha_lower + alpha_upper), 3)) {
      warning("The computed confidence interval does not have the same coverage as the specified confidence interval.",
        call. = FALSE
      )
    }
    return(.as_dmar_tbl(Result, conf_level = conf_level))
  }
}
