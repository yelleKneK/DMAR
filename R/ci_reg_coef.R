#' Confidence Interval for a Regression Coefficient, Raw or Standardized
#'
#' @description
#' The general engine behind \code{\link{ci_rc}} (unstandardized) and
#' \code{\link{ci_src}} (standardized): computes a confidence interval
#' for a population regression coefficient by the standard
#' \emph{t}-based approach or the noncentral \emph{t} approach, in
#' whichever metric the inputs are supplied.
#'
#' @param b_j Value of the regression coefficient for the \emph{j}th predictor variable
#' @param SE_b_j Standard error for the \emph{j}th predictor variable
#' @param s_Y Standard deviation of \emph{Y}, the dependent variable
#' @param s_X Standard deviation of \eqn{X_j}, the predictor variable of interest
#' @param N Sample size
#' @param p The number of predictors
#' @param R2_Y_X The squared multiple correlation coefficient predicting \code{Y} from the \code{p} predictor variables
#' @param R2_j_X_without_j The squared multiple correlation coefficient predicting the \emph{\code{j}}th predictor variable (i.e., the predictor of interest) from the remaining \code{p}-1 predictor variables
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param R2_Y_X_without_j The squared multiple correlation coefficient predicting \code{Y} from the \emph{\code{p}}-1 predictor variable with the \code{j}th predictor of interest excluded
#' @param t_value The \emph{t}-value evaluating the null hypothesis that the population regression coefficient for the \code{j}th predictor equals zero
#' @param alpha_lower The Type I error rate for the lower confidence interval limit
#' @param alpha_upper The Type I error rate for the upper confidence interval limit
#' @param noncentral \code{TRUE} or \code{FALSE}, specifying whether or not the noncentral approach to confidence intervals should be used
#' @param ... Optional additional specifications for nested functions
#'
#' @details
#' For standardized variables, do not specify the standard deviation of the variables and input the standardized
#' regression coefficient for \code{b_j}.
#'
#' When \code{b_j} is reconstructed from squared multiple correlations (that is,
#' from \code{R2_Y_X}, \code{R2_Y_X_without_j}, and \code{R2_j_X_without_j}
#' rather than a supplied \code{b_j}, \code{SE_b_j}, or \code{t_value}), only the
#' magnitude of the coefficient is identifiable; its sign is not. The positive
#' root is returned and a warning is issued. If the coefficient is negative,
#' negate the point estimate and swap and negate the confidence limits, or
#' supply \code{b_j} directly.
#'
#' @examples
#' ci_reg_coef(b_j = 0.61319, SE_b_j = 0.16098, N = 30, p = 6)
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The rows are ordered
#' \code{"lower_limit"}, \code{"reg_coef"} (the regression coefficient point
#' estimate), and \code{"upper_limit"}, so the point estimate sits between
#' its confidence limits. The lower and upper rows give the confidence limits
#' on the regression coefficient. The \code{prob_less} and \code{prob_greater}
#' columns report the achieved tail probabilities at each limit when the
#' noncentral t method is used (they are \code{NA} for the \code{"reg_coef"}
#' estimate row).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical
#' Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
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
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Not all of the values need to be specified, only those that
#' contain all of the necessary information in order to compute the
#' confidence interval (options are thus given for the values that
#' need to be specified).
#'
#' The function \code{ci_rc} in DMAR also calculates the confidence interval
#' for the population (unstandardized) regression coefficient. The
#' function \code{ci_src} also calculates the confidence interval
#' for the population (standardized) regression coefficient. These two
#' functions perform the same tasks as \code{ci_reg_coef} does and
#' are preferred to it because of simpler arguments.
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{conf_limits_nct}}, \code{\link{ci_rc}}, \code{\link{ci_src}}
#'
#' @keywords htest
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_reg_coef <- function(b_j, SE_b_j = NULL, s_Y = NULL, s_X = NULL, N, p, R2_Y_X = NULL, R2_j_X_without_j = NULL,
                        conf_level = .95, R2_Y_X_without_j = NULL, t_value = NULL, alpha_lower = NULL,
                        alpha_upper = NULL, noncentral = FALSE, ...) {

  # Determine if NC was used instead of noncentral
  # tmp <- try(is.null(NC), silent=TRUE)
  # if(tmp==TRUE || tmp==FALSE) noncentral <- NC

  if (!is.null(t_value)) {
    obs_t <- t_value
  }

  if (!is.null(b_j) && !is.null(SE_b_j)) {
    obs_t <- b_j / SE_b_j
  }

  if (is.null(SE_b_j)) {
    {      if (is.null(t_value)) {
      if (is.null(s_Y) && is.null(s_X) && noncentral == TRUE) {
        s_Y <- 1
        s_X <- 1
      }

      if (is.null(s_Y) || is.null(s_X)) stop("You need to specify 's_Y' and 's_X'.")

      SE_b_j <- sqrt(((1 - R2_Y_X) / ((1 - R2_j_X_without_j) * (N - p - 1)))) * (s_Y / s_X)
    }
    if (!is.null(t_value) && !is.null(b_j)) SE_b_j <- b_j / t_value
    if (!is.null(t_value) && is.null(b_j)) {
      SE_b_j <- sqrt(((1 - R2_Y_X) / ((1 - R2_j_X_without_j) * (N - p - 1)))) * (s_Y / s_X)
      b_j <- obs_t * SE_b_j
    }    }
    if (is.null(t_value)) obs_t <- b_j / SE_b_j
  }


  if (is.null(b_j)) {
    b_j <- (((R2_Y_X - R2_Y_X_without_j) / (1 - R2_j_X_without_j))^.5) * (s_Y / s_X)
    # The R^2 inputs are all squared quantities, so they fix the magnitude of the
    # coefficient but not its sign. The positive root is returned; warn loudly
    # rather than silently place the interval on the wrong side of zero.
    warning("The sign of 'b_j' is not identifiable from the R-squared inputs; the positive root is returned. If the coefficient is negative, negate the point estimate and swap (and negate) the confidence limits, or supply 'b_j' directly (optionally with 'SE_b_j' or 't_value').", call. = FALSE)
    if (is.null(SE_b_j)) {
      SE_b_j <- sqrt(((1 - R2_Y_X) / ((1 - R2_j_X_without_j) * (N - p - 1)))) * (s_Y / s_X)
    }
    obs_t <- b_j / SE_b_j
  }


  # A supplied pair of alphas sets the coverage to 1 - alpha_lower - alpha_upper.
  # Drop the default conf_level so it is neither used to overwrite those alphas
  # nor attached as a mislabeled footer; reject an explicitly supplied conf_level
  # mixed with the alphas.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      stop("Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.", call. = FALSE)
    }
    conf_level <- NULL
  }

  if (!is.null(conf_level)) {
    if (conf_level >= 1 || conf_level <= 0) stop("You have not properly specified \'conf_level\'", call. = FALSE)

    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }
  if (is.null(conf_level)) {
    if (is.null(alpha_lower) || is.null(alpha_upper)) stop("You need to specify either \'conf_level\', or both \'alpha_lower\' and \'alpha_upper\'.", call. = FALSE)
    if (alpha_lower > .5 || alpha_lower < 0) stop("You have not properly specified \'alpha_lower\' correctly.", call. = FALSE)
    if (alpha_upper > .5 || alpha_upper < 0) stop("You have not properly specified \'alpha_upper\' correctly.", call. = FALSE)
  }

  if (noncentral == FALSE) {
    term <- c("lower_limit", "reg_coef", "upper_limit")
    value <- c(b_j + qt(alpha_lower, df = N - p - 1) * SE_b_j, b_j, b_j + qt(1 - alpha_upper, df = N - p - 1) * SE_b_j)
    prob_less <- c(alpha_lower, NA_real_, 1 - alpha_upper)
    prob_greater <- c(1 - alpha_lower, NA_real_, alpha_upper)
    out <- data.frame(term, value, prob_less, prob_greater); attr(out, "conf_level") <- conf_level; class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }

  if (noncentral == TRUE) {
    NC_t_values <- conf_limits_nct(ncp = obs_t, df = N - p - 1, conf_level = NULL, alpha_lower = alpha_lower, alpha_upper = alpha_upper)
    term <- c("lower_limit", "reg_coef", "upper_limit")
    value <- c(NC_t_values[which(NC_t_values$term == "lower_limit"), 2] * SE_b_j, b_j, NC_t_values[which(NC_t_values$term == "upper_limit"), 2] * SE_b_j)
    # conf_limits_nct reports prob_less / prob_greater as the mass of the
    # distribution below / above each limit; convert to the achieved-tail-error
    # convention used by the central branch, ci_R2, and MBESS so that
    # prob_less on the lower row is the lower-tail error, not 1 - it.
    achieved_alpha_lower <- NC_t_values[which(NC_t_values$term == "lower_limit"), 4]
    achieved_alpha_upper <- NC_t_values[which(NC_t_values$term == "upper_limit"), 3]
    prob_less <- c(achieved_alpha_lower, NA_real_, 1 - achieved_alpha_upper)
    prob_greater <- c(1 - achieved_alpha_lower, NA_real_, achieved_alpha_upper)
    out <- data.frame(term, value, prob_less, prob_greater); attr(out, "conf_level") <- conf_level; class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }
}
