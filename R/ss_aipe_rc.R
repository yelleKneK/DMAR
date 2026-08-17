#' Sample Size Necessary for the Accuracy in Parameter Estimation Approach for an Unstandardized Regression Coefficient of Interest
#'
#' @description
#' A function used to plan sample size from the accuracy in parameter estimation perspective for an
#' unstandardized regression coefficient of interest given the input specification.
#'
#' @param rho2_Y_X Population value of the squared multiple correlation coefficient
#' @param Rho2_j_X_without_j Population value of the squared multiple correlation coefficient predicting the \emph{j}th predictor variable from the remaining \emph{p}-1 predictor variables
#' @param p The number of predictor variables
#' @param b_j The regression coefficient for the \emph{j}th predictor variable (i.e., the predictor of interest)
#' @param width The desired width of the confidence interval
#' @param which_width Which Width (\code{"Full"}, \code{"Lower"}, or \code{"Upper"}) the width refers to (at present, only \code{"Full"} can be specified)
#' @param sigma_Y The population standard deviation of \emph{Y} (i.e., the dependent variables)
#' @param sigma_X_j The population standard deviation of the \emph{j}th \emph{X} variable (i.e., the predictor variable of interest)
#' @param rho_XX Population correlation matrix for the \emph{p} predictor variables
#' @param rho_YX Population \emph{p} length vector of correlation between the dependent variable (\emph{Y}) and the \emph{p} independent variables
#' @param which_predictor Identifies which of the \emph{p} predictors is of interest
#' @param alpha_lower Type I error rate for the lower confidence interval limit
#' @param alpha_upper Type I error rate for the upper confidence interval limit
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param assurance Degree of certainty that the obtained confidence interval will be sufficiently narrow
#'
#' @details
#' Not all of the arguments need to be specified, only those that provide all of the necessary information
#' so that the sample size can be determined for the conditions specified.
#'
#' @return
#' Returns the necessary sample size in order for the goals of accuracy in parameter estimation to be
#' satisfied for the confidence interval for a particular regression coefficient given the input specifications.
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
#' \code{\link{ss_aipe_reg_coef_sensitivity}}, \code{\link{conf_limits_nct}}, \code{\link{ss_aipe_reg_coef}}, \code{\link{ss_aipe_src}}
#'
#' @examples
#' # Exchangeable correlation structure
#' rho_YX <- c(.3, .3, .3, .3, .3)
#' rho_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
#' c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))
#'
#' ss_aipe_rc(width = .1, which_width = "Full", sigma_Y = 1, sigma_X = 1, rho_XX = rho_XX,
#'            rho_YX = rho_YX, which_predictor = 1, conf_level = 1 - .05)
#'
#' ss_aipe_rc(width = .1, which_width = "Full", sigma_Y = 1, sigma_X = 1, rho_XX = rho_XX,
#'            rho_YX = rho_YX, which_predictor = 1, conf_level = 1 - .05, assurance = .85)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_rc <- function(rho2_Y_X = NULL, Rho2_j_X_without_j = NULL, p = NULL, b_j = NULL, width, which_width = "Full",
                       sigma_Y = 1, sigma_X_j = 1, rho_XX = NULL, rho_YX = NULL, which_predictor = NULL, alpha_lower = NULL,
                       alpha_upper = NULL, conf_level = .95, assurance = NULL) {
  return(ss_aipe_reg_coef(
    rho2_Y_X = rho2_Y_X, rho2_j_X_without_j = Rho2_j_X_without_j, p = p, b_j = b_j, width = width,
    which_width = which_width, sigma_Y = sigma_Y, sigma_X = sigma_X_j, rho_XX = rho_XX, rho_YX = rho_YX,
    which_predictor = which_predictor, noncentral = FALSE, alpha_lower = alpha_lower, alpha_upper = alpha_upper,
    conf_level = conf_level, assurance = assurance
  ))
}
