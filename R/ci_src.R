#' Confidence Interval for a Standardized Regression Coefficient
#'
#' Computes a confidence interval for a population standardized
#' regression coefficient, by the standard \emph{t}-based approach or
#' the noncentral \emph{t} approach. A thin convenience wrapper around
#' \code{\link{ci_reg_coef}}, which is the general engine; for the
#' coefficient in its raw metric use \code{\link{ci_rc}}.
#'
#' @param beta_j The standardized regression coefficient
#' @param SE_beta_j The standard error of the standardized regression coefficient
#' @param N Sample size
#' @param p The number of predictors
#' @param R2_Y_X The squared multiple correlation coefficient predicting \emph{Y} from the \emph{p} predictor variables
#' @param R2_j_X_without_j The squared multiple correlation coefficient predicting the \emph{j}th predictor variable (i.e., the predictor of interest) from the remaining \emph{p}-1 predictor variables
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param R2_Y_X_without_j The squared multiple correlation coefficient predicting \emph{Y} from the \emph{p}-1 predictor variable with the \emph{j}th predictor of interest excluded
#' @param t_value The \emph{t}-value evaluating the null hypothesis that the population regression coefficient for the \emph{j}th predictor equals zero
#' @param b_j The unstandardized regression coefficient
#' @param SE_b_j The standard error of the unstandardized regression coefficient
#' @param s_Y Standard deviation of \emph{Y}, the dependent variable
#' @param s_X Standard deviation of \emph{X}, the predictor variable of interest
#' @param alpha_lower The Type I error rate for the lower confidence interval limit
#' @param alpha_upper The Type I error rate for the upper confidence interval limit
#' @param \dots Optional additional specifications for nested functions
#'
#' @details
#' For standardized variables, do not specify the standard deviation of the variables and input the
#' standardized regression coefficient for \code{b_j}.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The \code{term} values are
#' \code{"lower_limit"}, \code{"src"} (the standardized regression
#' coefficient point estimate), and \code{"upper_limit"}, so the estimate
#' sits between its confidence limits. The \code{prob_less} and
#' \code{prob_greater} columns report the achieved tail probabilities at
#' each limit when the noncentral t method is used (\code{NA} for the
#' estimate row).
#'
#' @examples
#' ci_src(beta_j = .6707, .1761, N = 30, p = 6, conf_level = .95)
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
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{conf_limits_nct}}, \code{\link{ci_reg_coef}}, \code{\link{ci_rc}}
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_src <- function(beta_j = NULL, SE_beta_j = NULL, N = NULL, p = NULL, R2_Y_X = NULL, R2_j_X_without_j = NULL,
                   conf_level = .95, R2_Y_X_without_j = NULL, t_value = NULL, b_j = NULL, SE_b_j = NULL, s_Y = NULL, s_X = NULL,
                   alpha_lower = NULL, alpha_upper = NULL, ...) {
  if (!is.null(b_j)) {
    if (is.null(s_Y)) stop("Since you have specified the unstandardized regression coefficient, you must also specify the standard deviation of Y (so that the function can compute the standardized regression coefficient).")
    if (is.null(s_X)) stop("Since you have specified the unstandardized regression coefficient, you must also specify the standard deviation of X (so that the function can compute the standardized regression coefficient).")

    beta_j <- b_j * (s_X / s_Y)

    if (!is.null(SE_b_j)) SE_beta_j <- SE_b_j * (s_X / s_Y)
  }


  if (beta_j > 1.1) warning("This function is only for standardized regression coefficients. Is your 'b_j' in standardized units (the observed value), although possible, seems quite large?", call. = FALSE)

  result <- ci_reg_coef(
    b_j = beta_j, SE_b_j = SE_beta_j, s_Y = 1, s_X = 1, N = N, p = p,
    R2_Y_X = R2_Y_X, R2_j_X_without_j = R2_j_X_without_j, conf_level = conf_level,
    R2_Y_X_without_j = R2_Y_X_without_j, t_value = t_value, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, noncentral = TRUE, ...
  )

  ll_i <- which(result$term == 'lower_limit')
  ul_i <- which(result$term == 'upper_limit')
  # Report the standardized regression coefficient (the point estimate) between
  # its limits, so the rows read lower_limit, src, upper_limit.
  term <- c('lower_limit', 'src', 'upper_limit')
  value <- c(result[ll_i, 2], beta_j, result[ul_i, 2])
  prob_less <- c(result[ll_i, 3], NA_real_, 1 - result[ul_i, 4])
  prob_greater <- c(1 - result[ll_i, 3], NA_real_, result[ul_i, 4])

  result_new <- data.frame(term, value, prob_less, prob_greater)
  return(.as_dmar_tbl(result_new, conf_level = conf_level))
}
