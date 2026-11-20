#' Sample Size Necessary for the Accuracy in Parameter Estimation Approach for an Unstandardized Regression Coefficient of Interest
#'
#' @description
#' A function used to plan sample size from the accuracy in parameter estimation perspective for an
#' unstandardized regression coefficient of interest given the input specification.
#'
#' @param true_R2_Y_X Population value of the squared multiple correlation coefficient
#' @param true_R2_j_X_without_j Population value of the squared multiple correlation coefficient predicting the \emph{j}th predictor variable from the remaining \emph{p}-1 predictor variables
#' @param p The number of predictor variables
#' @param true_reg_coef_j The population regression coefficient for the \emph{j}th predictor variable (i.e., the predictor of interest), in the
#'   units \code{true_sd_Y} and \code{true_sd_X} describe: the standardized coefficient times \code{true_sd_Y / true_sd_X}. On the
#'   covariance structure path it may be left out; when it is given it must equal that value (see \code{\link{ss_aipe_reg_coef}})
#' @param width The desired width of the confidence interval
#' @param which_width Which Width (\code{"Full"}, \code{"Lower"}, or \code{"Upper"}) the width refers to (at present, only \code{"Full"} can be specified)
#' @param true_sd_Y The population standard deviation of \emph{Y} (i.e., the dependent variable)
#' @param true_sd_X The population standard deviation of the \emph{j}th \emph{X} variable (i.e., the predictor variable of interest)
#' @param true_cor_XX Population correlation matrix for the \emph{p} predictor variables
#' @param true_cor_YX Population \emph{p} length vector of correlation between the dependent variable (\emph{Y}) and the \emph{p} independent variables
#' @param which_predictor Identifies which of the \emph{p} predictors is of interest
#' @param alpha_lower Type I error rate for the lower confidence interval limit, at least 1e-9 and at most .5; a tail area
#'   below 1e-9 is refused, as in \code{\link{ss_aipe_reg_coef}}
#' @param alpha_upper Type I error rate for the upper confidence interval limit, at least 1e-9 and at most .5
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate); a level above
#'   .999999998, whose tail areas are below 1e-9, is refused
#' @param assurance Optional probability with which the realized confidence interval is to be no wider than \code{width}; when
#'   \code{NULL} (the default) the plan is for the width of the interval at the expected values of the sample squared
#'   multiple correlations (see \code{\link{ss_aipe_reg_coef}}), not for the mean width over samples
#' @param assurance_method How the probability that the realized interval is no wider than \code{width} is found, used
#'   only with \code{assurance}. Only \code{"analytic"} (the default) is offered here: the probability is computed exactly
#'   (an \emph{F} probability) and reported as \code{assurance_attained}, and the plan is the smallest sample size, at or
#'   above the plan without \code{assurance}, at which it is at least \code{assurance}. \code{"simulation"}, which
#'   \code{\link{ss_aipe_sc_ancova}}, \code{\link{ss_aipe_reg_coef}} with \code{noncentral = TRUE}, and
#'   \code{\link{ss_aipe_src}} offer, where the probability has no exact form, is refused with a condition of class
#'   \code{dmar_argument_invalid}. See "The assurance step" in \code{\link{ss_aipe_reg_coef}}.
#'
#' @details
#' Not all of the arguments need to be specified, only those that provide all of the necessary information
#' so that the sample size can be determined for the conditions specified.
#'
#' The interval planned for is the central \emph{t} interval for the coefficient in the units of the data, the
#' interval \code{\link{ci_rc}} reports, for random predictors from a multivariate normal population, with the
#' \emph{t} quantiles of its two tails. The
#' plan without \code{assurance} is the smallest sample size, of at least \eqn{p + 2}, at which that interval, at the
#' expected squared multiple correlations, is no wider than \code{width}; it can be a case or more smaller than that of
#' \pkg{MBESS} (see \code{\link{ss_aipe_reg_coef}}). With
#' \code{assurance} the plan is the smallest sample size at which that interval is no wider than \code{width} with at
#' least that probability, computed exactly (\code{assurance_method = "analytic"}). The published chi square
#' approximation of Kelley and Maxwell (2003) can give a smaller or a larger sample size than the exact plan.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}: \code{necessary_N}, the necessary total sample
#' size; with \code{assurance}, the rows \code{assurance} and \code{assurance_attained}, the exact probability at
#' \code{necessary_N} that the realized interval is no wider than \code{width}. With \code{assurance} the printed footer
#' names the assurance method, which the attribute \code{assurance_method} records, with \code{assurance_probability}
#' \code{"exact"}.
#'
#' @references
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 4 on individual
#'   comparisons of means and Chapter 6 on trend analysis.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' This function calls upon \code{ss_aipe_reg_coef} in DMAR but has a different naming scheme.
#' See \code{ss_aipe_reg_coef} for more details.
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef_sensitivity}}, \code{\link{ci_nc_t}}, \code{\link{ss_aipe_reg_coef}}, \code{\link{ss_aipe_src}}
#'
#' @examples
#' # Exchangeable correlation structure
#' true_cor_YX <- c(.3, .3, .3, .3, .3)
#' true_cor_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
#' c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))
#'
#' ss_aipe_rc(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
#'            true_cor_XX = true_cor_XX,
#'            true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05)
#'
#' # The smallest N at which the interval is no wider than .1 with
#' # probability at least .85; the last row is that probability.
#' ss_aipe_rc(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
#'            true_cor_XX = true_cor_XX,
#'            true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05, assurance = .85)
#'
#' # In the units of the data (Y with standard deviation 10, the predictor
#' # with standard deviation 2), a width of .1 * 10 / 2 gives the same plan.
#' ss_aipe_rc(width = .5, true_sd_Y = 10, true_sd_X = 2, true_cor_XX = true_cor_XX,
#'            true_cor_YX = true_cor_YX, which_predictor = 1, assurance = .85)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_rc <- function(true_R2_Y_X = NULL, true_R2_j_X_without_j = NULL, p = NULL, true_reg_coef_j = NULL, width, which_width = "Full",
                       true_sd_Y = 1, true_sd_X = 1, true_cor_XX = NULL, true_cor_YX = NULL, which_predictor = NULL, alpha_lower = NULL,
                       alpha_upper = NULL, conf_level = .95, assurance = NULL, assurance_method = "analytic") {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_rc", c("true_sd_Y", "true_sd_X", "which_predictor"))
  # The probability is exact here, so "simulation" is not offered
  # (R/ss_aipe_assurance_internals.R).
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_rc", simulation = FALSE, both = FALSE)
  # The engine's conditions are stated in this function's name.
  .reg_coef_in_wrapper_names(
    ss_aipe_reg_coef(
      true_R2_Y_X = true_R2_Y_X, true_R2_j_X_without_j = true_R2_j_X_without_j, p = p, true_reg_coef_j = true_reg_coef_j, width = width,
      which_width = which_width, true_sd_Y = true_sd_Y, true_sd_X = true_sd_X, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
      which_predictor = which_predictor, noncentral = FALSE, alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      conf_level = conf_level, assurance = assurance, assurance_method = assurance_method
    ),
    c(ss_aipe_reg_coef = "ss_aipe_rc"))
}
