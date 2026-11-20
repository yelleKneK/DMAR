#' Confidence Interval for a Standardized Regression Coefficient
#'
#' Computes a confidence interval for a population standardized
#' regression coefficient, by the standard \emph{t}-based approach or
#' the noncentral \emph{t} approach. A thin convenience wrapper around
#' \code{\link{ci_reg_coef}}, which is the general engine; for the
#' coefficient in its raw metric use \code{\link{ci_rc}}.
#'
#' @param std_reg_coef_j The standardized regression coefficient, or \code{NULL}
#'   (the default) when it is given in the units of the data as \code{reg_coef_j}
#'   or when the other inputs determine it (see Details)
#' @param se_std_reg_coef_j The standard error of the standardized regression coefficient
#' @param N Sample size
#' @param p The number of predictors
#' @param R2_Y_X The squared multiple correlation coefficient predicting \emph{Y} from the \emph{p} predictor variables
#' @param R2_j_X_without_j The squared multiple correlation coefficient predicting the \emph{j}th predictor variable (i.e., the predictor of interest) from the remaining \emph{p}-1 predictor variables
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate). A level
#'   above .999999998, whose tail areas are below 1e-9, is refused
#' @param R2_Y_X_without_j The squared multiple correlation coefficient predicting \emph{Y} from the \emph{p}-1 predictor variable with the \emph{j}th predictor of interest excluded
#' @param t_value The \emph{t}-value evaluating the null hypothesis that the population regression coefficient for the \emph{j}th predictor equals zero
#' @param reg_coef_j The unstandardized regression coefficient
#' @param se_reg_coef_j The standard error of the unstandardized regression coefficient
#' @param sd_Y Standard deviation of \emph{Y}, the dependent variable
#' @param sd_X Standard deviation of \emph{X}, the predictor variable of interest
#' @param alpha_lower The Type I error rate for the lower confidence interval limit
#' @param alpha_upper The Type I error rate for the upper confidence interval limit.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound. A nonzero tail
#'   area below 1e-9 is refused (see "Small tail areas" in
#'   \code{\link{ci_nc_t}}).
#'
#' @details
#' For standardized variables, do not specify the standard deviation of the variables and input the
#' standardized regression coefficient for \code{std_reg_coef_j}.
#'
#' The interval needs \code{N}, \code{p}, and one of these: the coefficient
#' (\code{std_reg_coef_j}, or \code{reg_coef_j} with \code{sd_Y} and
#' \code{sd_X}) with its standard error (\code{se_std_reg_coef_j} or
#' \code{se_reg_coef_j}), with \code{t_value}, or with \code{R2_Y_X} and
#' \code{R2_j_X_without_j}; \code{t_value} with \code{se_std_reg_coef_j},
#' or with \code{R2_Y_X} and \code{R2_j_X_without_j}; or the three squared
#' multiple correlations \code{R2_Y_X}, \code{R2_Y_X_without_j}, and
#' \code{R2_j_X_without_j}. The squared multiple correlations fix the
#' magnitude of the coefficient but not its sign: the positive root is
#' returned, with a warning, as in \code{\link{ci_reg_coef}}. A call that
#' gives none of these is refused with a condition of class
#' \code{dmar_argument_invalid}. The interval is computed by
#' \code{\link{ci_reg_coef}}, and its conditions name the arguments of this
#' function (\code{std_reg_coef_j}, not \code{reg_coef_j}). A standardized
#' coefficient larger than 1.1 in absolute value is reported with a warning
#' of class \code{dmar_result_qualified}, since a coefficient that large is
#' more often one in the units of the data.
#'
#' \strong{Fixed and random predictors.} The standard error of the
#' standardized coefficient, whether given as \code{se_std_reg_coef_j},
#' rescaled from \code{se_reg_coef_j} by \code{sd_X / sd_Y}, or computed
#' from the squared multiple correlations, treats the sample standard
#' deviations of \emph{Y} and of \eqn{X_j} as fixed constants (Kelley &
#' Maxwell, 2003, pp. 315--316). The standard deviation of \emph{Y} varies
#' from sample to sample in every design, and those of the predictors vary
#' as well when the predictors are random, as they are in most
#' observational studies, so this standard error is approximate either way
#' and the interval does not have exactly its nominal coverage. When the
#' values of the predictors are set by the design, the standard error is,
#' in large samples, too large whenever the population coefficient is not
#' 0, so the interval covers above its nominal level. When the predictors
#' are random, the standard error can be too large or too small, and which
#' it is depends on the correlations among \emph{Y} and the predictors,
#' not on \eqn{R^2} alone. In 2,193 two-predictor structures drawn with
#' each of the three correlations uniform between -.9 and .9 (structures
#' whose correlation matrix had an eigenvalue below .02 were set aside),
#' the large-sample standard error of the first coefficient was too small
#' in a quarter of them: in none of the 743 in which the predictors
#' explained at most 30 percent of the variance of \emph{Y}, in 8 percent
#' of those explaining 30 to 60 percent, 56 percent of those explaining 60
#' to 80, and 79 percent of those explaining 80 to 90, and in all 70 that
#' explained more than 90 percent. A structure outside that range can
#' reverse the last: with one predictor correlated .96 with \emph{Y}, a
#' second correlated .1 with \emph{Y}, and the two uncorrelated
#' (\eqn{R^2 = .93}), the large-sample standard error of the first
#' coefficient is 2.2 times its sampling standard deviation, and at
#' \emph{N} = 200 both intervals covered it in every one of 4,000 samples.
#'
#' In a simulation of multivariate normal data with random predictors (20
#' conditions with 1 to 4 predictors and \emph{N} from 50 to 1,000; 4,000
#' or 20,000 samples per condition, with Monte Carlo standard errors of
#' .001 to .005), the 95 percent interval of this function covered the
#' population standardized coefficient .961 to .998 of the time with a
#' single predictor correlated .3 to .7 with \emph{Y}; .964 to .979 for
#' coefficients of .34 to .52 among two to four predictors correlated
#' positively with \emph{Y} and with one another and explaining less than
#' 60 percent of its variance; .952 to .956 in the five conditions whose
#' coefficient was 0 or near it (at most .17 in absolute value); and .955
#' to .965 for suppressors (a predictor
#' correlated 0 to .1 with \emph{Y} and .6 to .8 with another). With two
#' uncorrelated predictors that explained 72 and 90 percent of the variance
#' of \emph{Y}, it covered .982 and .937. The central \emph{t} interval
#' with the same standard error (that of \code{\link{ci_reg_coef}} with
#' \code{noncentral = FALSE}) covered less: .882 to .952 for the
#' suppressors, and .938 and .881 in the last two conditions. With the
#' values of the predictors held fixed from sample to sample instead (six
#' conditions), this function's interval covered .971 to .998 whenever the
#' coefficient was not 0.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The \code{term} values are
#' \code{"ci_lower"}, \code{"src"} (the standardized regression
#' coefficient point estimate), and \code{"ci_upper"}, so the estimate
#' sits between its confidence limits. The \code{prob_less} and
#' \code{prob_greater} columns report the achieved tail probabilities at
#' each limit when the noncentral t method is used (\code{NA} for the
#' estimate row).
#'
#' @examples
#' ci_src(std_reg_coef_j = .6707, .1761, N = 30, p = 6, conf_level = .95)
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 4 on individual
#'   comparisons of means and Chapter 6 on trend analysis.)
#'
#' Smithson, M. (2003). \emph{Confidence intervals}. Thousand Oaks, CA: Sage Publications.
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close
#' fit in the Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' This function calls upon \code{ci_reg_coef} in DMAR, but has a different naming scheme.
#' See  \code{\link{ci_reg_coef}} for more details.
#'
#' To form a confidence interval for the unstandardized regression coefficient, use \code{ci_rc}.
#' This function is used to form a confidence interval for the standardized regression coefficient.
#'
#' Not all of the values need to be specified, only those that contain all of the necessary information
#' in order to compute the confidence interval (options are thus given for the values that need to be specified).
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{ci_nc_t}}, \code{\link{ci_reg_coef}}, \code{\link{ci_rc}}
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_src <- function(std_reg_coef_j = NULL, se_std_reg_coef_j = NULL, N = NULL, p = NULL, R2_Y_X = NULL, R2_j_X_without_j = NULL,
                   conf_level = .95, R2_Y_X_without_j = NULL, t_value = NULL, reg_coef_j = NULL, se_reg_coef_j = NULL, sd_Y = NULL, sd_X = NULL,
                   alpha_lower = NULL, alpha_upper = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_src", c("std_reg_coef_j", "se_std_reg_coef_j", "N", "p"))
  if (!is.null(reg_coef_j)) {
    if (is.null(sd_Y)) .dmar_refuse_argument(
      "Since you have specified the unstandardized regression coefficient, you must also specify the standard deviation of Y (so that the function can compute the standardized regression coefficient).",
      arg = "sd_Y")
    if (is.null(sd_X)) .dmar_refuse_argument(
      "Since you have specified the unstandardized regression coefficient, you must also specify the standard deviation of X (so that the function can compute the standardized regression coefficient).",
      arg = "sd_X")

    std_reg_coef_j <- reg_coef_j * (sd_X / sd_Y)

    if (!is.null(se_reg_coef_j)) se_std_reg_coef_j <- se_reg_coef_j * (sd_X / sd_Y)
  }

  if (is.null(N) || is.null(p)) .dmar_refuse_argument(
    "'N' and 'p' must both be given: the interval is on N - p - 1 degrees of freedom.",
    arg = c("N", "p"))

  # Without the coefficient, the routes the help names: the t statistic with
  # its standard error, the t statistic with R2_Y_X and R2_j_X_without_j
  # (which give the standard error), or the three squared multiple
  # correlations, which give the magnitude of the coefficient but not its
  # sign. Anything less is refused here rather than reaching R's own
  # "argument is of length zero" (the 1.1.0 development builds compared a
  # missing coefficient with 1.1 and stopped there).
  if (is.null(std_reg_coef_j)) {
    if (!is.null(t_value) && !is.null(se_std_reg_coef_j)) {
      std_reg_coef_j <- t_value * se_std_reg_coef_j
    } else {
      by_t <- !is.null(t_value) && !is.null(R2_Y_X) && !is.null(R2_j_X_without_j)
      by_R2 <- is.null(t_value) && !is.null(R2_Y_X) && !is.null(R2_Y_X_without_j) &&
        !is.null(R2_j_X_without_j)
      if (!by_t && !by_R2) .dmar_refuse_argument(
        "Without 'std_reg_coef_j' (or 'reg_coef_j' with 'sd_Y' and 'sd_X'), give 't_value' with 'se_std_reg_coef_j', 't_value' with 'R2_Y_X' and 'R2_j_X_without_j', or all three of 'R2_Y_X', 'R2_Y_X_without_j', and 'R2_j_X_without_j'.",
        arg = "std_reg_coef_j")
    }
  } else if (is.null(se_std_reg_coef_j) && is.null(t_value) &&
             (is.null(R2_Y_X) || is.null(R2_j_X_without_j))) {
    .dmar_refuse_argument(
      "With 'std_reg_coef_j' given, its standard error comes from 'se_std_reg_coef_j', from 't_value', or from 'R2_Y_X' and 'R2_j_X_without_j'; give one of them.",
      arg = c("se_std_reg_coef_j", "t_value", "R2_Y_X", "R2_j_X_without_j"))
  }

  # The engine's reg_coef_j and se_reg_coef_j are this function's
  # std_reg_coef_j and se_std_reg_coef_j; its conditions are stated in this
  # function's names.
  result <- .reg_coef_in_wrapper_names(
    ci_reg_coef(
      reg_coef_j = std_reg_coef_j, se_reg_coef_j = se_std_reg_coef_j, sd_Y = 1, sd_X = 1, N = N, p = p,
      R2_Y_X = R2_Y_X, R2_j_X_without_j = R2_j_X_without_j, conf_level = conf_level,
      R2_Y_X_without_j = R2_Y_X_without_j, t_value = t_value, alpha_lower = alpha_lower,
      alpha_upper = alpha_upper, noncentral = TRUE
    ),
    c(reg_coef_j = "std_reg_coef_j", se_reg_coef_j = "se_std_reg_coef_j"))

  ll_i <- which(result$term == 'ci_lower')
  ul_i <- which(result$term == 'ci_upper')
  # The estimate is the engine's, which is the coefficient given or the one
  # the other inputs determine.
  estimate <- result$value[result$term == 'reg_coef']

  # A standardized coefficient can exceed 1 in absolute value, but rarely
  # does; a large one in either direction is more often a coefficient in
  # the units of the data.
  if (abs(estimate) > 1.1) .dmar_warn(
    c("The standardized coefficient is ", format(estimate, digits = 4),
      ". A standardized regression coefficient can be larger than 1 in ",
      "absolute value, but rarely is; check that 'std_reg_coef_j' is in ",
      "standardized units (ci_rc() is the interval for a coefficient in ",
      "the units of the data)."),
    "result_qualified")

  # Report the standardized regression coefficient (the point estimate) between
  # its limits, so the rows read ci_lower, src, ci_upper.
  term <- c('ci_lower', 'src', 'ci_upper')
  value <- c(result[ll_i, 2], estimate, result[ul_i, 2])
  prob_less <- c(result[ll_i, 3], NA_real_, 1 - result[ul_i, 4])
  prob_greater <- c(1 - result[ll_i, 3], NA_real_, result[ul_i, 4])

  result_new <- data.frame(term, value, prob_less, prob_greater)
  return(.as_dmar_tbl(result_new, conf_level = attr(result, "conf_level", exact = TRUE),
                      tail_areas = attr(result, "tail_areas", exact = TRUE)))
}
