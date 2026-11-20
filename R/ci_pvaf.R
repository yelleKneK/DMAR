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
#' The confidence level must be specified in one of following two ways: using confidence interval coverage (\code{conf_level}),
#' or lower and upper confidence limits (\code{alpha_lower} and \code{alpha_upper}).
#'
#' This function uses the confidence interval transformation principle (Steiger, 2004) to transform the confidence limits for
#' the noncentrality parameter to the confidence limits for the population proportion of variance accounted for by knowing the group status.
#' The confidence interval for the noncentral \emph{F} parameter can be obtained from the function \code{\link{ci_nc_F}}, which is used within this function.
#'
#' \strong{A lower limit of 0.} When the observed \emph{F} is below the
#' \code{alpha_lower} critical value of the central \emph{F} distribution
#' (an \emph{F}-test that is not significant at \code{alpha_lower}), no
#' positive noncentrality parameter leaves as little as \code{alpha_lower}
#' above the observed \emph{F}, and the lower limit is 0, with a warning.
#' The interval from 0 to the upper limit has the stated coverage,
#' \code{1 - alpha_lower - alpha_upper}. A lower limit misses only when it
#' lies above the population value, which happens with probability
#' \code{alpha_lower} whatever that value is, and a lower limit of 0 never
#' lies above it, so a limit of 0 gives up no coverage. In 20,000 simulated
#' analyses of variance of 5 groups of 11 at each of five population values
#' from 0 to .20, the 95\% interval contained the population value in
#' 95.0\% to 95.1\% of them (Monte Carlo standard error 0.15 percentage
#' points), although its lower limit was 0 in 27\% to 98\% of them. An
#' analysis whose \emph{F} was so small that the central \emph{F}
#' distribution leaves \code{alpha_upper} or less below it has no upper
#' limit (the \code{"ci_upper"} and \code{"actual_coverage"} rows are
#' \code{NA}) and was counted as not containing the population value:
#' 2.5\% of the analyses at a population value of 0, and fewer at each
#' larger value.
#'
#' @return
#' A 4-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The \code{term} values are
#' \code{"ci_lower"} (the lower confidence limit on the proportion of
#' variance accounted for, on the [0, 1] scale), \code{"pvaf"} (the sample
#' proportion of variance accounted for,
#' \code{df_1 * F_value / (df_1 * F_value + df_2)}, the same value that eta
#' squared reports, so the point estimate sits between its confidence
#' limits), \code{"ci_upper"} (the upper confidence limit), and
#' \code{"actual_coverage"} (the coverage probability of the interval,
#' \code{1 - alpha_lower - alpha_upper}; \code{NA} when the upper limit is
#' undefined). The \code{prob_less} and \code{prob_greater} columns hold, on
#' the row of each limit, the two tail areas at the observed \code{F_value}
#' of the noncentral \emph{F} distribution whose noncentrality parameter
#' corresponds to that limit: \code{prob_less} is the area at or above
#' \code{F_value} and \code{prob_greater} the area at or below it (the
#' reverse of the columns of the same names in \code{\link{ci_nc_F}}, where
#' \code{prob_less} is the area at or below \code{F_value}), so
#' \code{prob_less} is \code{alpha_lower} on the \code{"ci_lower"} row and
#' \code{prob_greater} is \code{alpha_upper} on the \code{"ci_upper"} row.
#' With a tail area of 0 no limit is computed on that side (the lower limit
#' is 0, or the upper limit is 1), and the row of that limit holds the
#' requested area, 0, and its complement, 1. When \code{alpha_lower} is
#' larger than 0 and the observed \emph{F} is below the \code{alpha_lower}
#' critical value, so that the lower limit is 0, the \code{"ci_lower"} row
#' holds the tail areas of the central \emph{F} distribution, whose area at
#' or above \code{F_value} is the \emph{p}-value of the \emph{F}-test and
#' is larger than \code{alpha_lower}; those areas describe the central
#' distribution and not the coverage of the interval. Both columns are
#' \code{NA} on the \code{"pvaf"} and \code{"actual_coverage"} rows, and
#' on the \code{"ci_upper"} row when the upper limit is undefined.
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
#' Venables, W. (1975). Calculation of confidence intervals for
#'   noncentrality parameters. \emph{Journal of the Royal Statistical
#'   Society, Series B, 37}(3), 406--412.
#'   \doi{10.1111/j.2517-6161.1975.tb01554.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note This function can be used for single or factorial ANOVA designs.
#'
#' @seealso \code{\link{ci_nc_F}}
#'
#' @examples
#' ## Bargmann (1970) gave an example of a one-way ANOVA with an observed
#' ## F value of 11.2213 on 4 and 50 degrees of freedom, which Venables (1975),
#' ## Fleishman (1980), and Steiger (2004) used in turn. Steiger describes it
#' ## as five groups of 11 and reports F as 11.221, the value used here.
#' ## The exact confidence interval for the proportion of variance accounted for of that example:
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
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_pvaf", c("df_1", "df_2", "N"))
  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_upper - alpha_lower
  }
  if (!is.null(alpha_lower) && is.null(alpha_upper)) {
    .dmar_refuse_argument(
      "This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').",
      arg = c("alpha_lower", "alpha_upper"))
  }
  if (is.null(alpha_lower) && !is.null(alpha_upper)) {
    .dmar_refuse_argument(
      "This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').",
      arg = c("alpha_upper", "alpha_lower"))
  }
  if (is.null(df_1) || is.null(df_2) || is.null(N)) {
    .dmar_refuse_argument(
      "You need to specify 'df_1', 'df_2', and 'N'.",
      arg = c("df_1", "df_2", "N"))
  }
  if (is.null(F_value)) {
    .dmar_refuse_argument(
      "You must specify the observed F-value ('F_value') from the analysis of variance.",
      arg = "F_value")
  }
  if (alpha_lower > 0.5 || alpha_lower < 0) {
    .dmar_refuse_argument(
      "'alpha_lower' must be between 0 and .5.",
      arg = "alpha_lower")
  }
  if (conf_level > 1 || conf_level < 0) {
    .dmar_refuse_argument(
      "'conf_level' must be between 0 and 1.",
      arg = "conf_level")
  }
  if (F_value <= 0) {
    .dmar_refuse_argument("'F_value' must be larger than 0.", arg = "F_value")
  }
  if (N <= 0 || N <= df_1 + df_2) {
    .dmar_refuse_argument("'N' must be larger than 'df_1' + 'df_2'.", arg = c("N", "df_1", "df_2"))
  }

  Lims <- .ci_nc_F_for(
    caller = "ci_pvaf",
    quantity = "the proportion of variance accounted for",
    F_value = F_value, conf_level = NULL,
    df_1 = df_1, df_2 = df_2, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, ...
  )

  lower_lambda <- Lims$value[Lims$term == "ci_lower"]
  upper_lambda <- Lims$value[Lims$term == "ci_upper"]
  achieved_alpha_lower <- Lims$prob_greater[Lims$term == "ci_lower"]
  achieved_alpha_upper <- Lims$prob_less[Lims$term == "ci_upper"]

  # A lower limit of 0 costs no coverage: the limit misses only when it
  # lies above the population value, which a limit of 0 never does, and the
  # observed F exceeds the alpha_lower critical value of the F distribution
  # at the population value with probability alpha_lower. The tail area of
  # the central F at the observed F (achieved_alpha_lower there) is kept in
  # the tail columns and is not the interval's chance of missing, so the
  # coverage row uses the requested tail area.
  miss_lower <- if (lower_lambda == 0) alpha_lower else achieved_alpha_lower
  Actual_Coverage <- 1 - (miss_lower + achieved_alpha_upper)

  Lower_lim <- lower_lambda / (lower_lambda + N)
  Upper_lim <- if (is.infinite(upper_lambda)) 1 else upper_lambda / (upper_lambda + N)

  # The sample proportion of variance accounted for, in the closed form
  # implied by the F statistic; the same point estimate eta squared reports.
  Pvaf <- df_1 * F_value / (df_1 * F_value + df_2)

  term <- c("ci_lower", "pvaf", "ci_upper", "actual_coverage")
  value <- c(Lower_lim, Pvaf, Upper_lim, Actual_Coverage)
  prob_less <- c(achieved_alpha_lower, NA_real_, 1 - achieved_alpha_upper,
                 NA_real_)
  prob_greater <- c(1 - achieved_alpha_lower, NA_real_, achieved_alpha_upper,
                    NA_real_)

  out <- data.frame(term, value, prob_less, prob_greater)
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ci_long",
               tail_areas = .tail_areas(alpha_lower, alpha_upper))
}
