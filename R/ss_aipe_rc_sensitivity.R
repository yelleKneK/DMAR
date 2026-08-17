#' @title Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter
#' Estimation Perspective for the Unstandardized Regression Coefficient
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter
#' Estimation Perspective for the unstandardized regression coefficient.
#'
#' @param true_var_Y Population variance of the dependent variable (\emph{Y})
#' @param true_cov_YX Population covariances vector between the \emph{p} predictor variables and the dependent variable (\emph{Y})
#' @param true_cov_XX Population covariance matrix of the \emph{p} predictor variables
#' @param estimated_var_Y Estimated variance of the dependent variable (\emph{Y})
#' @param estimated_cov_YX Estimated covariances vector between the \emph{p} predictor variables and the dependent variable (\emph{Y})
#' @param estimated_cov_XX Estimated Population covariance matrix of the \emph{p} predictor variables
#' @param specified_N Directly specified sample size (instead of planning one from the estimated covariance structure)
#' @param which_predictor identifies which of the \emph{p} predictors is of interest
#' @param w desired confidence interval width for the regression coefficient of interest
#' @param noncentral specify with a \code{TRUE} or \code{FALSE} statement whether or not the noncentral approach to sample size planning should be used
#' @param standardize specify with a \code{TRUE} or \code{FALSE} statement whether or not the regression coefficient will be standardized; default is \code{TRUE}
#' @param conf_level desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param assurance degree of certainty that the obtained confidence interval will be sufficiently narrow (i.e., the probability that the observed interval will be no larger than desired)
#' @param G the number of generations (i.e., replications) of the simulation within the function
#' @param print_iter specify with a \code{TRUE/FALSE} statement if the iteration number should be printed as the simulation within the function runs
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#'
#' @details Direct specification of \code{true_cov_YX} and \code{true_cov_XX} is necessary, even if one is interested in
#' a single regression coefficient, so that the covariance/correlation structure can be specified when the simulation
#' within the function runs.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis. This function
#' delegates to \code{\link{ss_aipe_reg_coef_sensitivity}} and inherits
#' its return structure: mean / median / SD summaries of the realized
#' unstandardized regression coefficient, the realized interval widths,
#' and the realized squared multiple correlation coefficient; the
#' proportion of intervals at or below the planning target
#' (\code{pct_ci_less_w}); the tail-specific and overall empirical
#' non-coverage rates (\code{pct_ci_miss_low}, \code{pct_ci_miss_high},
#' \code{total_type_I_error}), all proportions on the 0 to 1 scale; and
#' the input echoes (\code{total_N}, \code{p}, \code{which_predictor},
#' \code{true_b_j}, \code{estimated_b_j}, \code{width},
#' \code{conf_level}, and, when one was supplied, \code{assurance}). See
#' \code{\link{ss_aipe_reg_coef_sensitivity}} for the full row list.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
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
#' @note
#' Note that when the true and estimated covariance structures agree (\code{true_cov_YX} equals \code{estimated_cov_YX} and \code{true_cov_XX} equals \code{estimated_cov_XX}),
#' the results are not literally from a sensitivity analysis, rather the function performs a standard simulation
#' study. A simulation study can be helpful in order to determine if the sample size procedure
#' under or overestimates necessary sample size. See \code{ss_aipe_reg_coef_sensitivity} in DMAR for more details.
#'
#' @seealso \code{\link{ss_aipe_reg_coef_sensitivity}}, \code{\link{ss_aipe_src_sensitivity}}, \code{\link{ss_aipe_reg_coef}}, \code{\link{ci_reg_coef}}
#'
#' @examples
#' # Sensitivity analysis for an unstandardized regression coefficient
#' # with two correlated predictors. G is kept small here so the example
#' # runs quickly; raise G for a stable Monte Carlo summary.
#' set.seed(113)
#' Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' rho_YX  <- c(0.4, 0.3)
#' cov_YX  <- rho_YX
#' ss_aipe_rc_sensitivity(
#'   true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
#'   estimated_var_Y = 1, estimated_cov_YX = cov_YX, estimated_cov_XX = Sigma_X,
#'   which_predictor = 1, w = 0.20, conf_level = 0.95,
#'   G = 20, print_iter = FALSE
#' )
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_rc_sensitivity <- function(true_var_Y = NULL, true_cov_YX = NULL, true_cov_XX = NULL,
                                   estimated_var_Y = NULL, estimated_cov_YX = NULL, estimated_cov_XX = NULL, specified_N = NULL,
                                   which_predictor = 1, w = NULL, noncentral = FALSE, standardize = FALSE, conf_level = .95,
                                   assurance = NULL, G = 1000, print_iter = TRUE, save = FALSE,
                                   filename = "ss_aipe_rc_sensitivity_result.csv") {
  return(ss_aipe_reg_coef_sensitivity(
    true_var_Y = true_var_Y, true_cov_YX = true_cov_YX, true_cov_XX = true_cov_XX,
    estimated_var_Y = estimated_var_Y, estimated_cov_YX = estimated_cov_YX, estimated_cov_XX = estimated_cov_XX, specified_N = specified_N,
    which_predictor = which_predictor, w = w, noncentral = noncentral, standardize = standardize, conf_level = conf_level,
    assurance = assurance, G = G, print_iter = print_iter, save = save, filename = filename
  ))
}
