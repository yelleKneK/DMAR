#' Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Unstandardized Regression Coefficient
#'
#' @description
#' This function performs a sensitivity analysis when planning sample size from the Accuracy in Parameter
#' Estimation Perspective for the unstandardized regression coefficient.
#' \code{\link{ss_aipe_rc_sensitivity}} is the same analysis under the name
#' that pairs with \code{\link{ss_aipe_rc}}, and
#' \code{\link{ss_aipe_src_sensitivity}} is its counterpart for the
#' standardized regression coefficient.
#'
#' @param true_var_Y Population variance of the dependent variable (\emph{Y})
#' @param true_cov_YX Population covariances vector between the \code{p} predictor variables and the dependent variable (\emph{Y})
#' @param true_cov_XX Population covariance matrix of the \code{p} predictor variables
#' @param supposed_var_Y Supposed variance of the dependent variable (\emph{Y}), the planning value from which the sample size is planned
#' @param supposed_cov_YX Supposed covariances vector (planning values) between the \code{p} predictor variables and the dependent variable (\code{Y})
#' @param supposed_cov_XX Supposed covariance matrix (planning values) of the \code{p} predictor variables
#' @param specified_N Directly specified sample size (instead of planning one from the supposed covariance structure)
#' @param which_predictor Identifies which of the \emph{p} predictors is of interest
#' @param width Desired confidence interval width for the regression coefficient of interest, in the units of the data
#' @param noncentral Specify with a \code{TRUE} or \code{FALSE} statement whether or not the noncentral approach to sample size planning should be used.
#'   The noncentral plan is for a standardized coefficient, so when the
#'   sample size is planned from the supposed values, \code{noncentral = TRUE}
#'   needs supposed variances of 1 for \emph{Y} and the targeted predictor,
#'   where the two scales coincide; otherwise the call is refused (with a
#'   condition of class \code{dmar_argument_invalid}), and
#'   \code{\link{ss_aipe_src_sensitivity}} is the function for the
#'   standardized coefficient.
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param assurance Degree of certainty that the obtained confidence interval will be sufficiently narrow
#' @param assurance_method Passed to \code{\link{ss_aipe_reg_coef}} when
#'   the sample size is planned from the supposed values with
#'   \code{assurance}. \code{"analytic"} (the default) plans, for the central
#'   interval, the smallest sample size at which the exact probability of an
#'   interval no wider than \code{width} is at least \code{assurance}, and,
#'   for the noncentral interval, the approximation that \pkg{MBESS} uses,
#'   which does not compute that probability. \code{"simulation"}, offered with
#'   \code{noncentral = TRUE} only, plans the smallest sample size at which a
#'   simulated estimate of the probability (made by the planner from its
#'   default \code{G} of 20000 samples with the random number stream that
#'   \code{seed} sets) is at least \code{assurance}; with
#'   \code{noncentral = FALSE} it is refused with a condition of class
#'   \code{dmar_argument_invalid}, as \code{\link{ss_aipe_reg_coef}} refuses
#'   it. Not used with \code{specified_N} or without \code{assurance}.
#' @param G The number of generations (i.e., replications) of the simulation within the function
#' @param print_iter Specify with a \code{TRUE}/\code{FALSE} statement if the iteration number should be printed as the simulation within the function runs
#' @param filename Optional path of a CSV file to receive the per-replication
#'   results (the coefficient estimate \code{reg_coef_j}, its confidence limits
#'   \code{ci_lower} and \code{ci_upper}, the observed \eqn{R^2} \code{R2}, the
#'   standard error \code{se_reg_coef_j}, and the \emph{t} statistic
#'   \code{t_value}), appended when
#'   the file already exists and created otherwise; the default \code{NULL}
#'   writes nothing, and a throwaway run that wants the file should point it
#'   at \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' Direct specification of \code{true_cov_YX} and \code{true_cov_XX} is necessary, even if one is
#' interested in a single regression coefficient, so that the covariance/correlation structure can be
#' specified when the simulation within the function runs.
#'
#' \strong{The coefficient.} The coefficient is the unstandardized one, the
#' \emph{j}th element of \eqn{\Sigma_{XX}^{-1}\sigma_{XY}}, where
#' \eqn{\Sigma_{XX}} is \code{true_cov_XX} and \eqn{\sigma_{XY}} is
#' \code{true_cov_YX}; call it \eqn{b_j}. It is the true value the intervals
#' are judged against, and the supposed coefficient echoed beside it is the
#' same element computed from \code{supposed_cov_XX} and
#' \code{supposed_cov_YX}. Every replication is fit in the units of the data,
#' and a sample size planned from the supposed values is planned by
#' \code{\link{ss_aipe_reg_coef}} with the supposed standard deviations of
#' \emph{Y} and of the targeted predictor (the square roots of
#' \code{supposed_var_Y} and of the \emph{j}th diagonal element of
#' \code{supposed_cov_XX}), so that \code{width} is a width in the units of
#' the data. The standardized coefficient,
#' \eqn{\beta_j = b_j \sigma_{X_j} / \sigma_Y} (Kelley & Maxwell, 2003,
#' p. 316), is the subject of \code{\link{ss_aipe_src_sensitivity}}.
#'
#' The interval of each replication is that of \code{\link{ci_reg_coef}}
#' with the standard error of the fitted regression, central or, with
#' \code{noncentral = TRUE}, noncentral.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries
#' are: \code{mean_reg_coef_j}, \code{median_reg_coef_j},
#' \code{sd_reg_coef_j} (summaries of the realized regression coefficient
#' point estimates); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the realized
#' interval widths); \code{pct_ci_less_width} (proportion of intervals at or
#' below the planning target \code{width}); \code{pct_ci_miss_low} and
#' \code{pct_ci_miss_high} (tail-specific empirical non-coverage of the
#' population coefficient); \code{total_type_I_error} (overall empirical
#' non-coverage, the sum of the two tails); \code{mean_R2},
#' \code{median_R2}, \code{sd_R2} (summaries of the realized squared
#' multiple correlation coefficient); and the input echoes
#' \code{total_N} (the sample size evaluated), \code{p},
#' \code{which_predictor}, \code{true_reg_coef_j} and \code{supposed_reg_coef_j} (the
#' population and planning values of the targeted coefficient implied by
#' the supplied covariance structures), \code{width}, \code{conf_level},
#' and \code{assurance} (present only when an assurance was supplied). The
#' proportion and Type I error rows are proportions on the 0 to 1
#' scale, not percentages, so \code{total_type_I_error} is the sum of
#' \code{pct_ci_miss_low} and \code{pct_ci_miss_high}.
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
#' @keywords design
#'
#' @note
#' Note that when the true and supposed covariance structures agree
#' (\code{true_cov_YX} equals \code{supposed_cov_YX} and \code{true_cov_XX}
#' equals \code{supposed_cov_XX}), the results are not literally from a
#' sensitivity analysis, rather the function performs a standard
#' simulation study. A simulation study can be helpful in order to determine if the sample size procedure
#' under or overestimates necessary sample size.
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{ci_reg_coef}}
#'
#' @examples
#' # Sensitivity analysis for an unstandardized regression coefficient
#' # with two predictors at a modest R squared. The Monte Carlo loop is
#' # run with a small number of generations (G) here so the example is
#' # fast; use a larger G (for example G = 1000) in real applications.
#' set.seed(113)
#' cov_XX <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' cov_YX <- c(0.4, 0.3)
#' ss_aipe_reg_coef_sensitivity(
#'   true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
#'   supposed_var_Y = 1, supposed_cov_YX = cov_YX, supposed_cov_XX = cov_XX,
#'   which_predictor = 1, width = 0.20, conf_level = 0.95,
#'   G = 100, print_iter = FALSE
#' )
#'
#' # The outcome in units whose variance is 4. The unstandardized
#' # coefficient of the first predictor is about .68 units of Y per unit of
#' # the predictor, and the 0.55 width is in those units. For the same
#' # structure ss_aipe_src_sensitivity() reports the standardized
#' # coefficient, .68 times 1 / 2, about .34, with a width in standardized
#' # units.
#' ss_aipe_reg_coef_sensitivity(
#'   true_var_Y = 4, true_cov_YX = c(0.8, 0.6), true_cov_XX = cov_XX,
#'   specified_N = 200, which_predictor = 1, width = 0.55,
#'   G = 100, print_iter = FALSE, seed = 113
#' )
#' ss_aipe_src_sensitivity(
#'   true_var_Y = 4, true_cov_YX = c(0.8, 0.6), true_cov_XX = cov_XX,
#'   specified_N = 200, which_predictor = 1, width = 0.28,
#'   G = 100, print_iter = FALSE, seed = 113
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


ss_aipe_reg_coef_sensitivity <- function(true_var_Y = NULL, true_cov_YX = NULL, true_cov_XX = NULL,
                                         supposed_var_Y = NULL, supposed_cov_YX = NULL, supposed_cov_XX = NULL, specified_N = NULL,
                                         which_predictor = 1, width = NULL, noncentral = FALSE, conf_level = .95,
                                         assurance = NULL, G = 1000, print_iter = TRUE, filename = NULL,
                                         assurance_method = c("analytic", "simulation"), seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!is.null(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_reg_coef_sensitivity", "which_predictor")
  if (!is.null(G)) G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)
  # The unstandardized coefficient: the engine in
  # R/reg_coef_sensitivity_internals.R names its rows reg_coef_j.
  .ss_aipe_reg_coef_sens(
    "ss_aipe_reg_coef_sensitivity", scale = "raw", assurance_method = assurance_method,
    true_var_Y = true_var_Y, true_cov_YX = true_cov_YX, true_cov_XX = true_cov_XX,
    supposed_var_Y = supposed_var_Y, supposed_cov_YX = supposed_cov_YX,
    supposed_cov_XX = supposed_cov_XX, specified_N = specified_N,
    which_predictor = which_predictor, width = width, noncentral = noncentral,
    conf_level = conf_level, assurance = assurance, G = G,
    print_iter = print_iter, filename = filename)
}
