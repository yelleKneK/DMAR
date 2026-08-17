#' Confidence Interval for an Unstandardized Regression Coefficient
#'
#' @description
#' Computes a confidence interval for a population regression
#' coefficient in its raw (unstandardized) metric, by the standard
#' \emph{t}-based approach or the noncentral \emph{t} approach. A thin
#' convenience wrapper around \code{\link{ci_reg_coef}}, which is the
#' general engine; for the standardized coefficient use
#' \code{\link{ci_src}}.
#'
#' @param b_j Value of the regression coefficient for the \emph{j}th predictor variable
#' @param SE_b_j Standard error for the \emph{j}th predictor variable
#' @param s_Y Standard deviation of \emph{Y}, the dependent variable
#' @param s_X Standard deviation of \emph{X}, the predictor variable of interest
#' @param N Sample size
#' @param p The number of predictors
#' @param R2_Y_X The squared multiple correlation coefficient predicting \emph{Y} from the \emph{p} predictor variables
#' @param R2_j_X_without_j The squared multiple correlation coefficient predicting the \emph{j}th predictor variable (i.e., the predictor of interest) from the remaining \emph{p}-1 predictor variables
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param R2_Y_X_without_j The squared multiple correlation coefficient predicting \emph{Y} from the \emph{p}-1 predictor variable with the \emph{j}th predictor of interest excluded
#' @param t_value The \emph{t}-value evaluating the null hypothesis that the population regression coefficient for the \emph{j}th predictor equals zero
#' @param alpha_lower The Type I error rate for the lower confidence interval limit
#' @param alpha_upper The Type I error rate for the upper confidence interval limit
#' @param noncentral \code{TRUE} or \code{FALSE} statement specifying whether or not the noncentral approach to confidence intervals should be used
#' @param \dots Optional additional specifications for nested functions
#'
#' @details
#' Returns the confidence limits for the regression coefficient of interest from the standard
#' approach to confidence interval formation or from the noncentral approach to confidence interval
#' formation using the noncentral \emph{t}-distribution.
#'
#' @return
#' A 2-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The \code{term} values are
#' \code{"lower_limit"} and \code{"upper_limit"}, and \code{value} holds
#' the confidence limits on the regression coefficient in its raw metric.
#' The \code{prob_less} and \code{prob_greater} columns report the tail
#' probabilities below and above each limit; when the noncentral \emph{t}
#' approach is used they are the achieved tail probabilities. Unlike
#' \code{\link{ci_src}} and \code{\link{ci_reg_coef}}, which place the
#' point estimate between its limits as a third row, \code{ci_rc} returns
#' the two limits only.
#'
#' @examples
#' ci_rc(b_j = 0.61319, SE_b_j = 0.16098, N = 30, p = 6, conf_level = 0.95)
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application,
#' and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
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
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the
#' Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Not all of the values need to be specified, only those that contain all of the necessary information
#' in order to compute the confidence interval (options are thus given for the values that need to be specified).
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{conf_limits_nct}}, \code{\link{ci_reg_coef}}, \code{\link{ci_src}}
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_rc <- function(b_j, SE_b_j = NULL, s_Y = NULL, s_X = NULL, N, p, R2_Y_X = NULL, R2_j_X_without_j = NULL,
                  conf_level = .95, R2_Y_X_without_j = NULL, t_value = NULL, alpha_lower = NULL, alpha_upper = NULL,
                  noncentral = FALSE, ...) {
  result <- ci_reg_coef(
    b_j = b_j, SE_b_j = SE_b_j, s_Y = s_Y, s_X = s_X, N = N, p = p, R2_Y_X = R2_Y_X,
    R2_j_X_without_j = R2_j_X_without_j, conf_level = conf_level,
    R2_Y_X_without_j = R2_Y_X_without_j, t_value = t_value, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, noncentral = noncentral, ...
  )
  term <- c("lower_limit", "upper_limit")
  value <- c(result[which(result$term == "lower_limit"), 2], result[which(result$term == "upper_limit"), 2])
  prob_less <- c(result[which(result$term == "lower_limit"), 3], 1 - result[which(result$term == "upper_limit"), 4])
  prob_greater <- c(1 - result[which(result$term == "lower_limit"), 3], result[which(result$term == "upper_limit"), 4])

  result_new <- data.frame(term, value, prob_less, prob_greater)
  return(.as_dmar_tbl(result_new, conf_level = conf_level))
}
