#' Confidence Interval for a Regression Coefficient, Raw or Standardized
#'
#' @description
#' The general engine behind \code{\link{ci_rc}} (unstandardized) and
#' \code{\link{ci_src}} (standardized): computes a confidence interval
#' for a population regression coefficient by the standard
#' \emph{t}-based approach or the noncentral \emph{t} approach, in
#' whichever metric the inputs are supplied.
#'
#' @param reg_coef_j Value of the regression coefficient for the \emph{j}th predictor variable, or
#'   \code{NULL} (the default) when the other inputs determine it: \code{t_value}
#'   with \code{R2_Y_X} and \code{R2_j_X_without_j}, or the three squared
#'   multiple correlations, which fix its magnitude but not its sign (see Details)
#' @param se_reg_coef_j Standard error for the \emph{j}th predictor variable
#' @param sd_Y Standard deviation of \emph{Y}, the dependent variable
#' @param sd_X Standard deviation of \eqn{X_j}, the predictor variable of interest
#' @param N Sample size
#' @param p The number of predictors
#' @param R2_Y_X The squared multiple correlation coefficient predicting \code{Y} from the \code{p} predictor variables
#' @param R2_j_X_without_j The squared multiple correlation coefficient predicting the \emph{\code{j}}th predictor variable (i.e., the predictor of interest) from the remaining \code{p}-1 predictor variables
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate). A level
#'   above .999999998, whose tail areas are below 1e-9, is refused
#' @param R2_Y_X_without_j The squared multiple correlation coefficient predicting \code{Y} from the \emph{\code{p}}-1 predictor variable with the \code{j}th predictor of interest excluded
#' @param t_value The \emph{t}-value evaluating the null hypothesis that the population regression coefficient for the \code{j}th predictor equals zero
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
#' @param noncentral \code{TRUE} or \code{FALSE}, specifying whether or not the noncentral approach to confidence intervals should be used
#'
#' @details
#' For standardized variables, do not specify the standard deviation of the variables and input the standardized
#' regression coefficient for \code{reg_coef_j}.
#'
#' \strong{Exactness.} For a coefficient in the units of the data, the
#' central interval (\code{noncentral = FALSE}, the default) is exact under
#' normal errors given the values of the predictors, and so for random
#' predictors as well. The noncentral interval multiplies the confidence
#' limits for the noncentrality parameter by the estimated standard error,
#' while the noncentrality parameter is the coefficient divided by its
#' population standard error, so that interval is approximate on either
#' scale. In a simulation with one predictor held fixed and normal errors
#' (4,000 samples per condition, \emph{N} of 15, 40, and 200, Monte Carlo
#' standard errors of .001 to .004), the 95 percent noncentral interval for
#' a coefficient in the units of the data covered .944 to .955 when the
#' noncentrality parameter was 1, and from .961 (\emph{N} = 200) to .993
#' (\emph{N} = 15) when it was 7; the central interval covered .943 to .953
#' in every condition.
#'
#' \strong{Standardized data.} For a coefficient of standardized data, or
#' one rescaled by the sample standard deviations, the standard error
#' treats those standard deviations as fixed constants (Kelley & Maxwell,
#' 2003, pp. 315--316). The standard deviation of \emph{Y} varies from
#' sample to sample in every design, and those of the predictors vary as
#' well when the predictors are random, so this standard error is
#' approximate either way and the interval does not have exactly its
#' nominal coverage. In the simulation of multivariate normal data
#' described in \code{\link{ci_src}}, the central 95 percent interval (the
#' default here) covered the population standardized coefficient .957 to
#' .993 of the time with a single random predictor correlated .3 to .7 with
#' \emph{Y}, .959 to .967 for coefficients of .34 to .52 among two to four
#' random predictors correlated positively with \emph{Y} and with one
#' another and explaining less than 60 percent of its variance, and .952
#' to .956 in the conditions whose coefficient was 0 or near it, but only
#' .882 to .952 for suppressors, and .938 and .881 with two uncorrelated
#' predictors that explained 72 and 90 percent of the variance of
#' \emph{Y}. With the values of the predictors held
#' fixed it covered .955 to .989 whenever the coefficient was not 0. The
#' noncentral interval, the one \code{\link{ci_src}} reports, covered more
#' in every condition except those whose coefficient was 0 or near it,
#' where the two agreed within Monte Carlo error.
#'
#' When \code{reg_coef_j} is reconstructed from squared multiple correlations (that is,
#' from \code{R2_Y_X}, \code{R2_Y_X_without_j}, and \code{R2_j_X_without_j}
#' rather than a supplied \code{reg_coef_j}, \code{se_reg_coef_j}, or \code{t_value}), only the
#' magnitude of the coefficient is identifiable; its sign is not. The positive
#' root is returned and a warning is issued. If the coefficient is negative,
#' negate the point estimate and swap and negate the confidence limits, or
#' supply \code{reg_coef_j} directly.
#'
#' @examples
#' ci_reg_coef(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098, N = 30, p = 6)
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The rows are ordered
#' \code{"ci_lower"}, \code{"reg_coef"} (the regression coefficient point
#' estimate), and \code{"ci_upper"}, so the point estimate sits between
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
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{ci_nc_t}}, \code{\link{ci_rc}}, \code{\link{ci_src}}
#'
#' @keywords htest
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_reg_coef <- function(reg_coef_j = NULL, se_reg_coef_j = NULL, sd_Y = NULL, sd_X = NULL, N, p, R2_Y_X = NULL, R2_j_X_without_j = NULL,
                        conf_level = .95, R2_Y_X_without_j = NULL, t_value = NULL, alpha_lower = NULL,
                        alpha_upper = NULL, noncentral = FALSE) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)

  # Determine if NC was used instead of noncentral
  # tmp <- try(is.null(NC), silent=TRUE)
  # if(tmp==TRUE || tmp==FALSE) noncentral <- NC

  if (!is.null(t_value)) {
    obs_t <- t_value
  }

  if (!is.null(reg_coef_j) && !is.null(se_reg_coef_j)) {
    obs_t <- reg_coef_j / se_reg_coef_j
  }

  if (is.null(se_reg_coef_j)) {
    {      if (is.null(t_value)) {
      if (is.null(sd_Y) && is.null(sd_X) && noncentral == TRUE) {
        sd_Y <- 1
        sd_X <- 1
      }

      if (is.null(sd_Y) || is.null(sd_X)) .dmar_refuse_argument(
        "You need to specify 'sd_Y' and 'sd_X'.", arg = c("sd_Y", "sd_X"))

      se_reg_coef_j <- sqrt(((1 - R2_Y_X) / ((1 - R2_j_X_without_j) * (N - p - 1)))) * (sd_Y / sd_X)
    }
    if (!is.null(t_value) && !is.null(reg_coef_j)) se_reg_coef_j <- reg_coef_j / t_value
    if (!is.null(t_value) && is.null(reg_coef_j)) {
      se_reg_coef_j <- sqrt(((1 - R2_Y_X) / ((1 - R2_j_X_without_j) * (N - p - 1)))) * (sd_Y / sd_X)
      reg_coef_j <- obs_t * se_reg_coef_j
    }    }
    if (is.null(t_value)) obs_t <- reg_coef_j / se_reg_coef_j
  }


  if (is.null(reg_coef_j)) {
    reg_coef_j <- (((R2_Y_X - R2_Y_X_without_j) / (1 - R2_j_X_without_j))^.5) * (sd_Y / sd_X)
    # The R^2 inputs are all squared quantities, so they fix the magnitude of the
    # coefficient but not its sign. The positive root is returned; warn loudly
    # rather than silently place the interval on the wrong side of zero.
    .dmar_warn(
      "The sign of 'reg_coef_j' is not identifiable from the R-squared inputs; the positive root is returned. If the coefficient is negative, negate the point estimate and swap (and negate) the confidence limits, or supply 'reg_coef_j' directly (optionally with 'se_reg_coef_j' or 't_value').",
      "result_qualified")
    if (is.null(se_reg_coef_j)) {
      se_reg_coef_j <- sqrt(((1 - R2_Y_X) / ((1 - R2_j_X_without_j) * (N - p - 1)))) * (sd_Y / sd_X)
    }
    obs_t <- reg_coef_j / se_reg_coef_j
  }


  # A supplied pair of alphas sets the coverage to 1 - alpha_lower - alpha_upper.
  # Drop the default conf_level so it is neither used to overwrite those alphas
  # nor attached as a mislabeled footer (the table records the coverage the
  # alphas define); reject an explicitly supplied conf_level mixed with the
  # alphas.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      .dmar_refuse_argument(
        "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    conf_level <- NULL
  }

  if (!is.null(conf_level)) {
    if (conf_level >= 1 || conf_level <= 0) .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")

    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }
  if (is.null(conf_level)) {
    if (is.null(alpha_lower) || is.null(alpha_upper)) .dmar_refuse_argument(
      "You need to specify either \'conf_level\', or both \'alpha_lower\' and \'alpha_upper\'.",
      arg = c("conf_level", "alpha_lower", "alpha_upper"))
    if (alpha_lower > .5 || alpha_lower < 0) .dmar_refuse_argument(
      "'alpha_lower' must be between 0 and .5.",
      arg = "alpha_lower")
    if (alpha_upper > .5 || alpha_upper < 0) .dmar_refuse_argument(
      "'alpha_upper' must be between 0 and .5.",
      arg = "alpha_upper")
  }
  # The floor on a tail area (R/dmar_conditions.R), for the central interval
  # as well as the noncentral one, whose limits come from ci_nc_t().
  .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = conf_level)

  if (noncentral == FALSE) {
    term <- c("ci_lower", "reg_coef", "ci_upper")
    # The upper quantile is taken in the upper tail (lower.tail = FALSE):
    # 1 - alpha_upper rounds a small tail area away.
    value <- c(reg_coef_j + qt(alpha_lower, df = N - p - 1) * se_reg_coef_j, reg_coef_j, reg_coef_j + qt(alpha_upper, df = N - p - 1, lower.tail = FALSE) * se_reg_coef_j)
    prob_less <- c(alpha_lower, NA_real_, 1 - alpha_upper)
    prob_greater <- c(1 - alpha_lower, NA_real_, alpha_upper)
    out <- data.frame(term, value, prob_less, prob_greater); attr(out, "conf_level") <- conf_level %||% (1 - alpha_lower - alpha_upper); attr(out, "tail_areas") <- .tail_areas(alpha_lower, alpha_upper); class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }

  if (noncentral == TRUE) {
    NC_t_values <- ci_nc_t(ncp = obs_t, df = N - p - 1, conf_level = NULL, alpha_lower = alpha_lower, alpha_upper = alpha_upper)
    term <- c("ci_lower", "reg_coef", "ci_upper")
    value <- c(NC_t_values[which(NC_t_values$term == "ci_lower"), 2] * se_reg_coef_j, reg_coef_j, NC_t_values[which(NC_t_values$term == "ci_upper"), 2] * se_reg_coef_j)
    # ci_nc_t reports prob_less / prob_greater as the mass of the
    # distribution below / above each limit; convert to the achieved-tail-error
    # convention used by the central branch, ci_R2, and MBESS so that
    # prob_less on the lower row is the lower-tail error, not 1 - it.
    achieved_alpha_lower <- NC_t_values[which(NC_t_values$term == "ci_lower"), 4]
    achieved_alpha_upper <- NC_t_values[which(NC_t_values$term == "ci_upper"), 3]
    prob_less <- c(achieved_alpha_lower, NA_real_, 1 - achieved_alpha_upper)
    prob_greater <- c(1 - achieved_alpha_lower, NA_real_, achieved_alpha_upper)
    out <- data.frame(term, value, prob_less, prob_greater); attr(out, "conf_level") <- conf_level %||% (1 - alpha_lower - alpha_upper); attr(out, "tail_areas") <- .tail_areas(alpha_lower, alpha_upper); class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }
}
