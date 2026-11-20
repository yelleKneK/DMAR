#' @title Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter
#' Estimation Perspective for the Standardized Regression Coefficient
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation
#' Perspective for the standardized regression coefficient.
#' \code{\link{ss_aipe_reg_coef_sensitivity}} and
#' \code{\link{ss_aipe_rc_sensitivity}} are the counterparts for the
#' unstandardized coefficient.
#'
#' @param true_var_Y Population variance of the dependent variable (\emph{Y})
#' @param true_cov_YX Population covariances vector between the \emph{p} predictor variables and the dependent variable (\emph{Y})
#' @param true_cov_XX Population covariance matrix of the \emph{p} predictor variables
#' @param supposed_var_Y Supposed variance of the dependent variable (\emph{Y}), the planning value from which the sample size is planned
#' @param supposed_cov_YX Supposed covariances vector (planning values) between the \emph{p} predictor variables and the dependent variable (\emph{Y})
#' @param supposed_cov_XX Supposed covariance matrix (planning values) of the \emph{p} predictor variables
#' @param specified_N Directly specified sample size (instead of planning one from the supposed covariance structure)
#' @param which_predictor identifies which of the \emph{p} predictors is of interest
#' @param width Desired confidence interval width for the standardized regression coefficient of interest
#' @param noncentral \code{TRUE} (the default) for the noncentral \emph{t}
#'   interval for the standardized coefficient (that of \code{\link{ci_src}}),
#'   \code{FALSE} for the central \emph{t} interval of Kelley and Maxwell
#'   (2003, Equation 1). The interval of each replication is that interval,
#'   and a sample size planned from the supposed values is planned for it by
#'   \code{\link{ss_aipe_src}} with the same \code{noncentral}.
#' @param conf_level desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param assurance degree of certainty that the obtained confidence interval will be sufficiently narrow
#' @param assurance_method Passed to \code{\link{ss_aipe_src}}, which plans
#'   the sample size from the supposed values with \code{assurance}:
#'   \code{"analytic"} (the default) for the published approximation, which
#'   does not compute the probability of an interval no wider than
#'   \code{width} (for the noncentral interval that of \pkg{MBESS}, and
#'   for the central one the chi square approximation of Kelley and
#'   Maxwell, 2003, Equation 3), at unit standard deviations;
#'   \code{"simulation"} for the smallest sample size at which a
#'   simulated estimate of that probability, made by the planner from its
#'   default \code{G} of 20000 samples with the random number stream that
#'   \code{seed} sets, is at least \code{assurance}. Not used with
#'   \code{specified_N} or without \code{assurance}.
#' @param G the number of generations/replication of the simulation study within the function
#' @param print_iter specify with a \code{TRUE/FALSE} statement if the iteration number should be printed as the simulation within the function runs
#' @param filename Optional path of a CSV file to receive the per-replication
#'   results (the standardized coefficient estimate \code{std_reg_coef_j},
#'   its confidence limits \code{ci_lower} and \code{ci_upper}, the observed
#'   \eqn{R^2} \code{R2}, the standard error \code{se_std_reg_coef_j}, and
#'   the \emph{t} statistic \code{t_value}), appended when the file already
#'   exists and created otherwise; the default \code{NULL} writes nothing,
#'   and a throwaway run that wants the file should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' Direct specification of \code{true_cov_YX} and \code{true_cov_XX} is necessary, even if one is interested in
#' a single regression coefficient, so that the covariance/correlation structure can be specified when
#' the simulation study within the function runs.
#'
#' \strong{The coefficient.} The data of every replication are standardized
#' with their own sample standard deviations before the regression is fit,
#' so each estimate and each interval is for the standardized coefficient,
#' the quantity \code{true_std_reg_coef_j} of \code{\link{ss_aipe_src}}
#' stands for. The true value the intervals are judged against is the
#' population standardized coefficient
#' \deqn{\beta_j = b_j \sigma_{X_j} / \sigma_Y,}
#' the unstandardized coefficient \eqn{b_j} (the \emph{j}th element of
#' \eqn{\Sigma_{XX}^{-1}\sigma_{XY}}, from \code{true_cov_XX} and
#' \code{true_cov_YX}) multiplied by the ratio of the population standard
#' deviations (Kelley & Maxwell, 2003, p. 316), where \eqn{\sigma_Y} is the
#' square root of \code{true_var_Y} and \eqn{\sigma_{X_j}} the square root
#' of the \emph{j}th diagonal element of \code{true_cov_XX}. The supposed
#' coefficient is put on the same scale with the supposed standard
#' deviations (\code{supposed_var_Y} and the diagonal of
#' \code{supposed_cov_XX}), and a sample size planned from the supposed
#' values is planned by \code{\link{ss_aipe_src}} for a coefficient on
#' that scale, with the supposed correlations and the same
#' \code{noncentral}, so that \code{width} is a width in the units of the
#' intervals simulated and the plan is for the interval simulated.
#'
#' \strong{The plan with assurance.} With \code{noncentral = FALSE} and
#' \code{assurance}, the plan is that of \code{\link{ss_aipe_src}} for the
#' central interval for the standardized coefficient: by default the chi
#' square approximation of Kelley and Maxwell (2003, Equation 3), and with
#' \code{assurance_method = "simulation"} the smallest sample size at which
#' that interval is no wider than \code{width} with at least that
#' probability. For two predictors correlated .3 with correlations of .4
#' and .3 with \emph{Y}, a width of .3, and an assurance of .90, the
#' approximation is 181, at which .986 of the standardized intervals are no
#' wider than .3, and the simulated plan is 171 (another seed can give one
#' a case away; see \code{\link{ss_aipe_src}}).
#'
#' The interval of each replication is that of \code{\link{ci_reg_coef}}
#' with the standard error of the regression fitted to the standardized
#' data. That standard error treats the sample standard deviations as
#' constants, so it is approximate, and the predictors simulated here are
#' random, so those of the predictors vary along with that of \emph{Y}
#' (Kelley & Maxwell, 2003, pp. 315--316; see \code{\link{ci_src}}); the
#' coverage rows show how good the approximation is for the structure given.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis: mean / median / SD
#' summaries of the realized standardized regression coefficient
#' (\code{mean_std_reg_coef_j}, \code{median_std_reg_coef_j},
#' \code{sd_std_reg_coef_j}), the realized interval widths
#' (\code{mean_ci_width}, \code{median_ci_width}, \code{sd_ci_width}), and
#' the realized squared multiple correlation coefficient (\code{mean_R2},
#' \code{median_R2}, \code{sd_R2}); the proportion of intervals at or below
#' the planning target (\code{pct_ci_less_width}); the tail-specific and
#' overall empirical non-coverage of the population standardized
#' coefficient (\code{pct_ci_miss_low}, \code{pct_ci_miss_high},
#' \code{total_type_I_error}), all proportions on the 0 to 1 scale; and the
#' input echoes (\code{total_N}, \code{p}, \code{which_predictor},
#' \code{true_std_reg_coef_j} and \code{supposed_std_reg_coef_j}, the
#' population and planning values of the standardized coefficient implied
#' by the covariance structures, \code{width}, \code{conf_level}, and, when
#' one was supplied, \code{assurance}). The per-replication file written
#' with \code{filename} has the columns \code{std_reg_coef_j},
#' \code{ci_lower}, \code{ci_upper}, \code{R2}, \code{se_std_reg_coef_j},
#' and \code{t_value}.
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
#' @note Note that when the true and supposed covariance structures agree (\code{true_cov_YX} equals \code{supposed_cov_YX} and \code{true_cov_XX} equals \code{supposed_cov_XX}), the results are not
#' literally from a sensitivity analysis, rather the function performs a standard simulation study. A simulation study can
#' be helpful in order to determine if the sample size procedure under or overestimates necessary sample size.
#' See \code{ss_aipe_reg_coef_sensitivity} in DMAR for more details.
#'
#' @seealso \code{\link{ss_aipe_src}}, \code{\link{ss_aipe_reg_coef_sensitivity}}, \code{\link{ss_aipe_rc_sensitivity}}, \code{\link{ss_aipe_reg_coef}}, \code{\link{ci_reg_coef}}
#'
#' @examples
#' # Sensitivity analysis for a standardized regression coefficient
#' # with two correlated predictors. A production run uses many more
#' # generations (G = 1000 is typical); G is reduced here so the
#' # example runs quickly.
#' set.seed(113)
#' cov_XX  <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' cov_YX  <- c(0.4, 0.3)
#' ss_aipe_src_sensitivity(
#'   true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
#'   supposed_var_Y = 1, supposed_cov_YX = cov_YX, supposed_cov_XX = cov_XX,
#'   which_predictor = 1, width = 0.20, conf_level = 0.95,
#'   G = 50, print_iter = FALSE
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


ss_aipe_src_sensitivity <- function(true_var_Y = NULL, true_cov_YX = NULL, true_cov_XX = NULL,
                                    supposed_var_Y = NULL, supposed_cov_YX = NULL, supposed_cov_XX = NULL, specified_N = NULL,
                                    which_predictor = 1, width = NULL, noncentral = TRUE, conf_level = .95,
                                    assurance = NULL, G = 1000, print_iter = TRUE, filename = NULL,
                                    assurance_method = c("analytic", "simulation"), seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!is.null(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_src_sensitivity", "which_predictor")
  if (!is.null(G)) G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)
  # The standardized coefficient: the engine in
  # R/reg_coef_sensitivity_internals.R standardizes every replication and
  # names the rows std_reg_coef_j.
  .ss_aipe_reg_coef_sens(
    "ss_aipe_src_sensitivity", scale = "std", assurance_method = assurance_method,
    true_var_Y = true_var_Y, true_cov_YX = true_cov_YX, true_cov_XX = true_cov_XX,
    supposed_var_Y = supposed_var_Y, supposed_cov_YX = supposed_cov_YX,
    supposed_cov_XX = supposed_cov_XX, specified_N = specified_N,
    which_predictor = which_predictor, width = width, noncentral = noncentral,
    conf_level = conf_level, assurance = assurance, G = G,
    print_iter = print_iter, filename = filename)
}
