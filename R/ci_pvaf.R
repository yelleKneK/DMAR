#' Confidence Interval for the Proportion of Variance Accounted for (in the Dependent Variable by Knowing the Levels of the Factor)
#'
#' Computes the exact confidence limits for the proportion of variance in the
#' dependent variable accounted for by knowing the levels of the factor (group
#' status in a single factor design) in a fixed effects analysis of variance,
#' so an omnibus \emph{F}-test is accompanied by an effect size with a
#' statement of its precision.
#'
#' @param F_value Observed \emph{F}-value from fixed effects analysis of variance
#' @param df_1 Numerator degrees of freedom
#' @param df_2 Denominator degrees of freedom
#' @param N Sample size
#' @param conf_level Confidence interval coverage (i.e., 1-Type I error rate); default is .95
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' The confidence level must be specified in one of following two ways: using confidence interval coverage (\code{conf_level}),
#' or lower and upper confidence limits (\code{alpha_lower} and \code{alpha_upper}).
#'
#' This function uses the confidence interval transformation principle (Steiger, 2004) to transform the confidence limits for
#' the noncentrality parameter to the confidence limits for the population proportion of variance accounted for by knowing the group status.
#' The confidence interval for the noncentral \emph{F} parameter can be obtained from the function \code{\link{conf_limits_ncf}}, which is used within this function.
#'
#' @return
#' A 4-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The \code{term} values are
#' \code{"lower_limit"} (the lower confidence limit on the proportion of
#' variance accounted for, on the [0, 1] scale), \code{"pvaf"} (the sample
#' proportion of variance accounted for,
#' \code{df_1 * F_value / (df_1 * F_value + df_2)}, the same value that eta
#' squared reports, so the point estimate sits between its confidence
#' limits), \code{"upper_limit"} (the upper confidence limit), and
#' \code{"actual_coverage"} (the achieved coverage probability, which equals
#' \code{conf_level} when both tail targets are met). The \code{prob_less}
#' and \code{prob_greater} columns report the achieved tail-error
#' probabilities at the two limits; \code{NA} on the \code{"pvaf"} and
#' \code{"actual_coverage"} rows.
#'
#' @references
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43},
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note This function can be used for single or factorial ANOVA designs.
#'
#' @seealso \code{\link{conf_limits_ncf}}
#'
#' @examples
#' ## Bargman (1970) gave an example in which a 5-group ANOVA with 11 subjects in each
#' ## group is conducted and the observed F value is 11.221. This example was used
#' ## in Venables (1975),  Fleishman (1980), and Steiger (2004). If one wants to calculate the
#' ## exact confidence interval for the proportion of variance accounted for in that example,
#' ## this function can be used.
#' ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#'
#' ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#'
#' ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = 0, alpha_upper = .05)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_pvaf <- function(F_value = NULL, df_1 = NULL, df_2 = NULL, N = NULL,
                    conf_level = 0.95, alpha_lower = NULL, alpha_upper = NULL,
                    ...) {
  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_upper - alpha_lower
  }
  if (!is.null(alpha_lower) && is.null(alpha_upper)) {
    stop("This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').")
  }
  if (is.null(alpha_lower) && !is.null(alpha_upper)) {
    stop("This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').")
  }
  if (is.null(df_1) || is.null(df_2) || is.null(N)) {
    stop("You need to specify 'df_1', 'df_2', and 'N'.")
  }
  if (is.null(F_value)) {
    stop("You must specify the observed F-value ('F_value') from the analysis of variance.")
  }
  if (alpha_lower > 0.5 || alpha_lower < 0) {
    stop(" 'alpha_lower' must be smaller than .5 and nonnegative.")
  }
  if (conf_level > 1 || conf_level < 0) {
    stop(" 'conf_level' must be larger than 0 and smaller than 1. ")
  }
  if (F_value <= 0) {
    stop(" 'F_value' must be larger than 0. ")
  }
  if (N <= 0 || N <= df_1 + df_2) {
    stop("N must be larger than df_1+df_2")
  }

  Lims <- .conf_limits_ncf_for(
    caller = "ci_pvaf",
    quantity = "the proportion of variance accounted for",
    F_value = F_value, conf_level = NULL,
    df_1 = df_1, df_2 = df_2, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, ...
  )

  lower_lambda <- Lims$value[Lims$term == "lower_limit"]
  upper_lambda <- Lims$value[Lims$term == "upper_limit"]
  achieved_alpha_lower <- Lims$prob_greater[Lims$term == "lower_limit"]
  achieved_alpha_upper <- Lims$prob_less[Lims$term == "upper_limit"]

  Actual_Coverage <- 1 - (achieved_alpha_lower + achieved_alpha_upper)

  Lower_lim <- lower_lambda / (lower_lambda + N)
  Upper_lim <- if (is.infinite(upper_lambda)) 1 else upper_lambda / (upper_lambda + N)

  # The sample proportion of variance accounted for, in the closed form
  # implied by the F statistic; the same point estimate eta squared reports.
  Pvaf <- df_1 * F_value / (df_1 * F_value + df_2)

  term <- c("lower_limit", "pvaf", "upper_limit", "actual_coverage")
  value <- c(Lower_lim, Pvaf, Upper_lim, Actual_Coverage)
  prob_less <- c(achieved_alpha_lower, NA_real_, 1 - achieved_alpha_upper,
                 NA_real_)
  prob_greater <- c(1 - achieved_alpha_lower, NA_real_, achieved_alpha_upper,
                    NA_real_)

  out <- data.frame(term, value, prob_less, prob_greater)
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ci_long")
}
