#' @title Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter
#' Estimation Perspective for the Unstandardized Regression Coefficient
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter
#' Estimation Perspective for the unstandardized regression coefficient. It
#' is \code{\link{ss_aipe_reg_coef_sensitivity}} under the name that pairs
#' with \code{\link{ss_aipe_rc}}; \code{\link{ss_aipe_src_sensitivity}} is
#' the counterpart for the standardized coefficient.
#'
#' @param true_var_Y Population variance of the dependent variable (\emph{Y})
#' @param true_cov_YX Population covariances vector between the \emph{p} predictor variables and the dependent variable (\emph{Y})
#' @param true_cov_XX Population covariance matrix of the \emph{p} predictor variables
#' @param supposed_var_Y Supposed variance of the dependent variable (\emph{Y}), the planning value from which the sample size is planned
#' @param supposed_cov_YX Supposed covariances vector (planning values) between the \emph{p} predictor variables and the dependent variable (\emph{Y})
#' @param supposed_cov_XX Supposed covariance matrix (planning values) of the \emph{p} predictor variables
#' @param specified_N Directly specified sample size (instead of planning one from the supposed covariance structure)
#' @param which_predictor identifies which of the \emph{p} predictors is of interest
#' @param width Desired confidence interval width for the regression coefficient of interest, in the units of the data
#' @param noncentral specify with a \code{TRUE} or \code{FALSE} statement whether or not the noncentral approach to sample size planning should be used; the noncentral plan is for a standardized coefficient, so a plan from the supposed values with \code{noncentral = TRUE} needs supposed variances of 1 for \emph{Y} and the targeted predictor and is refused otherwise (see \code{\link{ss_aipe_reg_coef_sensitivity}})
#' @param conf_level desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param assurance degree of certainty that the obtained confidence interval will be sufficiently narrow (i.e., the probability that the observed interval will be no larger than desired)
#' @param assurance_method Passed to \code{\link{ss_aipe_reg_coef}}, which
#'   plans the sample size from the supposed values with \code{assurance}:
#'   \code{"analytic"} (the default) for the smallest sample size at which
#'   the exact probability of an interval no wider than \code{width} is at
#'   least \code{assurance} with the central interval, and for the
#'   published approximation of \pkg{MBESS}, which does not compute that
#'   probability, with \code{noncentral = TRUE}; \code{"simulation"}, offered with
#'   \code{noncentral = TRUE} only, for the smallest sample size at which a
#'   simulated estimate of the probability (made by the planner from its
#'   default \code{G} of 20000 samples with the random number stream that
#'   \code{seed} sets) is at least \code{assurance}.
#'   With \code{noncentral = FALSE} \code{"simulation"} is refused with a
#'   condition of class \code{dmar_argument_invalid}. Not used with
#'   \code{specified_N} or without \code{assurance}.
#' @param G the number of generations (i.e., replications) of the simulation within the function
#' @param print_iter specify with a \code{TRUE/FALSE} statement if the iteration number should be printed as the simulation within the function runs
#' @param filename Optional path of a CSV file to receive the per-replication
#'   results (the unstandardized coefficient estimate \code{reg_coef_j}, its
#'   confidence limits \code{ci_lower} and \code{ci_upper}, the observed
#'   \eqn{R^2} \code{R2}, the standard error \code{se_reg_coef_j}, and the
#'   \emph{t} statistic \code{t_value}), appended when the file already
#'   exists and created otherwise; the default \code{NULL} writes nothing,
#'   and a throwaway run that wants the file should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details Direct specification of \code{true_cov_YX} and \code{true_cov_XX} is necessary, even if one is interested in
#' a single regression coefficient, so that the covariance/correlation structure can be specified when the simulation
#' within the function runs.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis, the table of
#' \code{\link{ss_aipe_reg_coef_sensitivity}}: mean / median / SD summaries
#' of the realized unstandardized regression coefficient
#' (\code{mean_reg_coef_j}, \code{median_reg_coef_j},
#' \code{sd_reg_coef_j}), the realized interval widths, and the realized
#' squared multiple correlation coefficient; the proportion of intervals at
#' or below the planning target (\code{pct_ci_less_width}); the
#' tail-specific and overall empirical non-coverage rates
#' (\code{pct_ci_miss_low}, \code{pct_ci_miss_high},
#' \code{total_type_I_error}), all proportions on the 0 to 1 scale; and
#' the input echoes (\code{total_N}, \code{p}, \code{which_predictor},
#' \code{true_reg_coef_j}, \code{supposed_reg_coef_j}, \code{width},
#' \code{conf_level}, and, when one was supplied, \code{assurance}). See
#' \code{\link{ss_aipe_reg_coef_sensitivity}} for the full row list and the
#' scale of each value.
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
#' Note that when the true and supposed covariance structures agree (\code{true_cov_YX} equals \code{supposed_cov_YX} and \code{true_cov_XX} equals \code{supposed_cov_XX}),
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
#' cov_XX <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' cor_YX <- c(0.4, 0.3)
#' cov_YX <- cor_YX
#' ss_aipe_rc_sensitivity(
#'   true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
#'   supposed_var_Y = 1, supposed_cov_YX = cov_YX, supposed_cov_XX = cov_XX,
#'   which_predictor = 1, width = 0.20, conf_level = 0.95,
#'   G = 20, print_iter = FALSE
#' )
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export


ss_aipe_rc_sensitivity <- function(true_var_Y = NULL, true_cov_YX = NULL, true_cov_XX = NULL,
                                   supposed_var_Y = NULL, supposed_cov_YX = NULL, supposed_cov_XX = NULL, specified_N = NULL,
                                   which_predictor = 1, width = NULL, noncentral = FALSE, conf_level = .95,
                                   assurance = NULL, G = 1000, print_iter = TRUE, filename = NULL,
                                   assurance_method = c("analytic", "simulation"), seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!is.null(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_rc_sensitivity", "which_predictor")
  if (!is.null(G)) G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)
  # The unstandardized coefficient, as in ss_aipe_reg_coef_sensitivity(); the
  # engine is in R/reg_coef_sensitivity_internals.R.
  .ss_aipe_reg_coef_sens(
    "ss_aipe_rc_sensitivity", scale = "raw", assurance_method = assurance_method,
    true_var_Y = true_var_Y, true_cov_YX = true_cov_YX, true_cov_XX = true_cov_XX,
    supposed_var_Y = supposed_var_Y, supposed_cov_YX = supposed_cov_YX,
    supposed_cov_XX = supposed_cov_XX, specified_N = specified_N,
    which_predictor = which_predictor, width = width, noncentral = noncentral,
    conf_level = conf_level, assurance = assurance, G = G,
    print_iter = print_iter, filename = filename)
}
