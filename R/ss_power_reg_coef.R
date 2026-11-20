#' Sample Size for a Targeted Regression Coefficient
#'
#' Determine the necessary sample size for a targeted regression coefficient or determine the degree of power given a specified sample size
#'
#' @param true_R2_Y_X Population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p} predictor variables (i.e., the \emph{X} variables)
#' @param true_R2_Y_X_without_j Population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p}-1 predictor variables, where the one not used is the predictor of interest
#' @param p Number of predictor variables in the full model, a single whole number of at least 1. It is required when the effect is given through \code{true_R2_Y_X}, \code{true_reg_coef_j}, or \code{true_f2}, and is not given with \code{true_cor_XX} and \code{true_cor_YX}, whose dimension sets it
#' @param desired_power Desired degree of statistical power for the test of targeted regression coefficient
#' @param alpha_level Type I error rate
#' @param directional Whether or not a direction or a nondirectional test is to be used (usually \code{directional=FALSE})
#' @param true_reg_coef_j Population value of the regression coefficient for the predictor of interest, the name it has in \code{\link{ss_aipe_reg_coef}} and \code{\link{ss_aipe_rc}}. It is on the scale that \code{true_sd_X} and \code{true_sd_Y} describe: with both standard deviations 1 it is the standardized coefficient
#' @param true_sd_X Population standard deviation for the predictor variable of interest
#' @param true_sd_Y Population standard deviation for the outcome variable
#' @param true_R2_j_X_without_j Population squared multiple correlation coefficient predicting the predictor variable of interest from the remaining p-1 predictor variables
#' @param true_cor_XX Population correlation matrix for the \code{p} predictor variables: symmetric, with 1 on the diagonal, and positive definite; another matrix is refused with a condition of class \code{dmar_argument_invalid}
#' @param true_cor_YX Population vector of correlation coefficient between the \code{p} predictor variables and the criterion variable; with \code{true_cor_XX} it must imply a squared multiple correlation below 1 (the correlation matrix of the outcome and the predictors is then positive definite)
#' @param true_f2 Cohen's (1988) definition for an effect size for a targeted regression coefficient: \code{(true_R2_Y_X-true_R2_Y_X_without_j)/(1-true_R2_Y_X)}
#' @param which_predictor Identifies the predictor of interest when \code{true_cor_XX} and \code{true_cor_YX} are specified: a single whole number from 1 to \code{p}, its position in \code{true_cor_XX}
#' @param specified_N Sample size for which power should be evaluated. This is the \emph{total} sample size.
#' @param print_iter If the progress of the iterative procedure is printed to the screen as the iterations are occurring
#'
#' @details
#' Determines the necessary sample size given a desired level of statistical power. Alternatively,
#' determines the statistical power for a given a specified sample size. There are a number of ways
#' that the specification regarding the size of the regression coefficient can be entered. The most
#' basic, and often the simplest, is to specify \code{true_R2_Y_X} and \code{true_R2_Y_X_without_j}.
#' See the examples section for several options.
#'
#' Power is computed from a noncentral \emph{t} distribution with noncentrality
#' \eqn{\sqrt{N}\,f}, which treats the predictors as fixed (their values held
#' constant across hypothetical replications). This is the standard
#' fixed-predictor power analysis; under random predictors, where the predictor
#' values themselves vary from sample to sample, the sample size required for a
#' given level of power is somewhat larger. The test is the case of
#' \code{\link{ss_power_R2_change}} with one predictor added and fixed
#' predictors (\code{true_R2_Y_X} is its \code{true_R2_full} and
#' \code{true_R2_Y_X_without_j} its \code{true_R2_restricted}); that function
#' also plans the test with random predictors and the test of a set of added
#' predictors.
#'
#' With \code{true_cor_XX} and \code{true_cor_YX}, the squared multiple
#' correlation of the full model is
#' \eqn{\rho_{YX}' \mathbf{R}_{XX}^{-1} \rho_{YX}}, where \eqn{\rho_{YX}} is
#' \code{true_cor_YX} and \eqn{\mathbf{R}_{XX}} is \code{true_cor_XX}, and that of the model
#' without the predictor of interest is the same quantity with the
#' predictor's row and column removed from \eqn{\mathbf{R}_{XX}} and its
#' element from \eqn{\rho_{YX}}. With a single predictor (\code{p} of 1)
#' the model without
#' it has no predictors and its squared multiple correlation is 0, so the
#' effect is that of \code{true_f2} \eqn{= \rho^2 / (1 - \rho^2)}, where
#' \eqn{\rho} is \code{true_cor_YX}.
#'
#' When the planning values give the predictor an effect of exactly 0
#' (\code{true_f2} of 0, equal squared multiple correlations, or a
#' \code{true_reg_coef_j} of 0), the null hypothesis is true and the \emph{t}
#' statistic has the central \emph{t} distribution at every sample size,
#' so the probability that the test rejects is \code{alpha_level}, for the
#' two-sided test and for the directional one alike. A sample size is then
#' defined only for a \code{desired_power} no greater than
#' \code{alpha_level}, and it is the smallest admissible one,
#' \code{p + 2}; a larger \code{desired_power} is refused at once with a
#' condition of class \code{dmar_effect_zero}. With \code{specified_N}
#' the power at that size is returned, which is \code{alpha_level}. For an
#' effect that is not 0 the power rises with the sample size, and a
#' \code{desired_power} that no sample size up to 10,000,000 attains is
#' refused at once with a condition of class \code{dmar_computation_failed}.
#'
#' @return
#' \item{ss}{Either the necessary sample size or the specified sample size, depending if one is interested in determining the necessary sample size given a desired degree of statistical power or if one is interested in the determining the value of statistical power given a specified sample size, respectively}
#' \item{actual_power}{Actual power of the situation described}
#' \item{noncentral_t_parm}{Value of the noncentral distribution for the appropriate \emph{t}-distribution}
#' \item{true_f}{Cohen's \emph{f} at the planning values, the effect size for the noncentral \emph{t}-distribution; this is the square root of \code{true_f2}, because \code{true_f2} is the effect size using an \emph{F}-distribution. When the effect is given through \code{true_reg_coef_j} it carries the sign of the coefficient.}
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
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
#' Maxwell, S. E. (2000). Sample size and multiple regression analysis.
#'   \emph{Psychological Methods, 5}(4), 434--458.
#'   \doi{10.1037/1082-989X.5.4.434}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 4 on individual
#'   comparisons of means and Chapter 6 on trend analysis.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
#'   planning for more accurate statistical power: A method adjusting
#'   sample effect sizes for publication bias and uncertainty.
#'   \emph{Psychological Science, 28}(11), 1547--1562.
#'   \doi{10.1177/0956797617723724}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{ss_power_R2}},
#' \code{\link{ss_power_R2_change}}, \code{\link{ci_nc_F}}
#'
#' @examples
#' Cor.Mat <- rbind(
#'   c(1.00, 0.53, 0.58, 0.60, 0.46, 0.66),
#'   c(0.53, 1.00, 0.35, 0.07, 0.14, 0.43),
#'   c(0.58, 0.35, 1.00, 0.18, 0.29, 0.50),
#'   c(0.60, 0.07, 0.18, 1.00, 0.30, 0.26),
#'   c(0.46, 0.14, 0.29, 0.30, 1.00, 0.30),
#'   c(0.66, 0.43, 0.50, 0.26, 0.30, 1.00)
#' )
#'
#' true_cor_XX <- Cor.Mat[2:6, 2:6]
#' true_cor_YX <- Cor.Mat[1, 2:6]
#'
#' # Method 1
#' ss_power_reg_coef(true_R2_Y_X = 0.7826786, true_R2_Y_X_without_j = 0.7363697, p = 5,
#'                   alpha_level = .05, directional = FALSE, desired_power = .80)
#'
#' # Method 2
#' ss_power_reg_coef(alpha_level = .05, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
#'                   which_predictor = 5,
#'                   directional = FALSE, desired_power = .80)
#'
#' # Method 3
#' # Here, true_reg_coef_j is the standardized regression coefficient. Had true_reg_coef_j
#' # been the unstandardized regression coefficient, true_sd_X and true_sd_Y
#' # would have been the standard deviation for the X variable of
#' # interest and Y, respectively.
#' ss_power_reg_coef(true_R2_Y_X = 0.7826786, true_R2_j_X_without_j = 0.3652136,
#'                   true_reg_coef_j = 0.2700964,
#'                   p = 5, alpha_level = .05, true_sd_X = 1, true_sd_Y = 1, directional = FALSE,
#'                   desired_power = .80)
#'
#' # Method 4
#' ss_power_reg_coef(alpha_level = .05, true_f2 = 0.2130898, p = 5,
#'                   directional = FALSE, desired_power = .80)
#'
#' # Power given a specified N and squared multiple correlation coefficients.
#' ss_power_reg_coef(true_R2_Y_X = 0.7826786, true_R2_Y_X_without_j = 0.7363697, specified_N = 25,
#'                   p = 5, alpha_level = .05, directional = FALSE)
#'
#' # Power given a specified N and effect size.
#' ss_power_reg_coef(alpha_level = .05, true_f2 = 0.2130898, p = 5, specified_N = 25,
#'                   directional = FALSE)
#'
#' # Reproducing Maxwell's (2000, p. 445) Example
#' Cor.Mat.Maxwell <- rbind(
#'   c(1.00, 0.35, 0.20, 0.20, 0.20, 0.20),
#'   c(0.35, 1.00, 0.40, 0.40, 0.40, 0.40),
#'   c(0.20, 0.40, 1.00, 0.45, 0.45, 0.45),
#'   c(0.20, 0.40, 0.45, 1.00, 0.45, 0.45),
#'   c(0.20, 0.40, 0.45, 0.45, 1.00, 0.45),
#'   c(0.20, 0.40, 0.45, 0.45, 0.45, 1.00)
#' )
#'
#' RHO.XX.Maxwell <- Cor.Mat.Maxwell[2:6, 2:6]
#' Rho.YX.Maxwell <- Cor.Mat.Maxwell[1, 2:6]
#' R2.Maxwell <- Rho.YX.Maxwell %*% solve(RHO.XX.Maxwell) %*% Rho.YX.Maxwell
#'
#' RHO.XX.Maxwell.no.1 <- Cor.Mat.Maxwell[3:6, 3:6]
#' Rho.YX.Maxwell.no.1 <- Cor.Mat.Maxwell[1, 3:6]
#' R2.Maxwell.no.1 <-
#'   Rho.YX.Maxwell.no.1 %*% solve(RHO.XX.Maxwell.no.1) %*% Rho.YX.Maxwell.no.1
#'
#' # This procedure arrives at N = 111, whereas Maxwell (2000, p. 445)
#' # reports N = 113. The two differ because of the noncentrality
#' # parameterization, not rounding: this function uses the fixed-predictor
#' # noncentrality sqrt(N) * f (see Details), while the tabled value rests on
#' # Cohen's (1988) convention. Neither is a random-predictor result; under
#' # random predictors, where the predictor values vary across replications,
#' # the sample size needed for the same power is larger still.
#' ss_power_reg_coef(true_R2_Y_X = R2.Maxwell, true_R2_Y_X_without_j = R2.Maxwell.no.1, p = 5,
#'                   alpha_level = .05, directional = FALSE, desired_power = .80)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export

ss_power_reg_coef <- function(true_R2_Y_X = NULL, true_R2_Y_X_without_j = NULL, p = NULL, desired_power = .85, alpha_level = .05, directional = FALSE, true_reg_coef_j = NULL, true_sd_X = NULL, true_sd_Y = NULL, true_R2_j_X_without_j = NULL, true_cor_XX = NULL, true_cor_YX = NULL, which_predictor = NULL, true_f2 = NULL, specified_N = NULL, print_iter = FALSE) {
  .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("ss_power_reg_coef", "p")
  # 'p' counts the predictors of the full model, so it is a single whole
  # number of at least 1, checked on every route before any search. Through
  # the first builds of 1.1.0 a fractional p was used as given (p = 2.5
  # returned necessary_N = 47.5), and p = 0, p = -1, Inf, and a vector
  # reached the search or stopped with R's own message.
  if (!is.null(p)) .reg_coef_check_p(p)

  # Validate 'desired_power' before the search: power tops out below 1, so a
  # 'desired_power' at or above 1 would drive the 'while (Dif > 0)' loop below
  # forever. Only enforced when a sample size is being sought.
  if (is.null(specified_N) || !is.null(desired_power)) {
    .dmar_check_desired_power(desired_power)
  }

  # 'f' is the noncentrality-defining effect size; each input method below sets
  # it. Initialize to NULL so a missing or inconsistent effect size
  # specification fails loudly rather than surfacing as 'object f not found'.
  f <- NULL

  # Method 1 to define the effect size.
  if (!is.null(true_R2_Y_X) && !is.null(true_R2_Y_X_without_j)) {
    if (!is.null(true_cor_XX) || !is.null(true_cor_YX)) .dmar_refuse_argument(
      "Since \'true_R2_Y_X,\' and \'true_R2_Y_X_without_j\', have been specified, do not specify \'true_cor_XX\' or \'true_cor_YX\'.",
      arg = c("true_R2_Y_X_without_j", "true_cor_XX", "true_cor_YX"))
    if (!is.null(which_predictor)) .dmar_refuse_argument(
      "Since you have specified \'true_R2_Y_X\' and \'true_R2_Y_X_without_j\', do not specify \'which_predictor\'.",
      arg = c("true_R2_Y_X", "true_R2_Y_X_without_j", "which_predictor"))
    f <- sqrt((true_R2_Y_X - true_R2_Y_X_without_j) / (1 - true_R2_Y_X))
  }


  # Method 2 to define the effect size.
  if (is.null(true_R2_Y_X) && is.null(true_R2_j_X_without_j) && is.null(true_R2_Y_X_without_j) && !is.null(true_cor_YX) && !is.null(true_cor_XX)) {
    if (is.null(which_predictor)) .dmar_refuse_argument(
      "Since you have specified \'true_cor_XX\' and \'true_cor_YX\', you must also specify the \'which_predictor\'.",
      arg = c("true_cor_XX", "true_cor_YX", "which_predictor"))
    if (!is.null(p)) .dmar_refuse_argument(
      "There is no need to specify \'p\' in this situation.", arg = "p")
    if (!is.numeric(true_cor_XX) || length(dim(true_cor_XX)) != 2L ||
        any(!is.finite(true_cor_XX)) || !is.numeric(true_cor_YX) ||
        any(!is.finite(true_cor_YX)) ||
        (dim(true_cor_XX)[1] != dim(true_cor_XX)[2]) ||
        (dim(true_cor_XX)[1] != length(true_cor_YX))) .dmar_refuse_argument(
      c("'true_cor_XX' must be a square numeric matrix, the correlation ",
        "matrix of the p predictors, and 'true_cor_YX' a numeric vector ",
        "of the p correlations of the predictors with the outcome, with ",
        "no missing values."),
      arg = c("true_cor_XX", "true_cor_YX"))
    p <- dim(true_cor_XX)[1]
    which_predictor <- .reg_coef_check_which(which_predictor, p)
    true_cor_YX <- as.numeric(true_cor_YX)
    .reg_coef_check_rho(true_cor_XX)
    true_R2_Y_X <- .reg_coef_R2(true_cor_XX, true_cor_YX)
    # The predictors and the outcome have a joint correlation matrix only
    # when the squared multiple correlation is below 1; otherwise the
    # correlations in 'true_cor_YX' cannot occur with 'true_cor_XX'.
    if (true_R2_Y_X >= 1) .dmar_refuse_argument(c(
      "'true_cor_XX' and 'true_cor_YX' are not correlations that can occur ",
      "together: the squared multiple correlation they imply is ",
      format(signif(true_R2_Y_X, 6)), ", and it must be below 1 (the ",
      "correlation matrix of the outcome and the predictors must be ",
      "positive definite)."), arg = c("true_cor_XX", "true_cor_YX"))
    # With one predictor the restricted model has no predictors, and its
    # squared multiple correlation is 0. Through the first builds of 1.1.0
    # this route solved a 0 by 0 matrix and stopped ("'a' is 0-diml").
    true_R2_Y_X_without_j <- if (p == 1L) 0 else .reg_coef_R2(
      true_cor_XX[-which_predictor, -which_predictor, drop = FALSE],
      true_cor_YX[-which_predictor])
    # For a positive definite 'true_cor_XX' the restricted R2 never exceeds
    # the full one; a difference of a few units in the last place below 0
    # is rounding and is read as 0, the effect of a predictor that adds
    # nothing (the zero-effect refusal below then applies).
    R2_change <- true_R2_Y_X - true_R2_Y_X_without_j
    if (R2_change < 0 && R2_change > -64 * .Machine$double.eps) R2_change <- 0
    f <- sqrt(R2_change / (1 - true_R2_Y_X))
  }

  # Method 3 to define the effect size.
  if (!is.null(true_R2_Y_X) && !is.null(true_reg_coef_j) && !is.null(true_R2_j_X_without_j)) {
    if (!is.null(which_predictor)) .dmar_refuse_argument(
      "Since you have specified \'true_R2_Y_X\', \'true_reg_coef_j\', and \'true_R2_j_X_without_j\', do not specify \'which_predictor\'.",
      arg = c("true_R2_Y_X", "true_reg_coef_j", "true_R2_j_X_without_j"))
    if (is.null(true_sd_X) || is.null(true_sd_Y)) .dmar_refuse_argument(
      "Since you have specified \'true_R2_Y_X\', \'true_reg_coef_j\', and \'true_R2_j_X_without_j\', you must also specify \'true_sd_X\' and \'true_sd_Y\'.",
      arg = c("true_R2_Y_X", "true_reg_coef_j", "true_R2_j_X_without_j"))
    f <- true_reg_coef_j * sqrt((1 - true_R2_j_X_without_j) / (1 - true_R2_Y_X)) * (true_sd_X / true_sd_Y)
  }

  # Method 4, input the effect size, a la Cohen, directly.
  if (!is.null(true_f2)) {
    if (!is.null(which_predictor)) .dmar_refuse_argument(
      "You do not need to specify \'which_predictor\'.",
      arg = "which_predictor")
    if (!is.null(true_sd_X) || !is.null(true_sd_Y)) .dmar_refuse_argument(
      "You do not need to specify \'true_sd_X\' or \'true_sd_Y\'.",
      arg = c("true_sd_X", "true_sd_Y"))
    if (is.null(p)) .dmar_refuse_argument(
      "You need to specify \'p\'.", arg = "p")
    if (!is.null(true_cor_XX) || !is.null(true_cor_YX)) .dmar_refuse_argument(
      "Since the effect size, \'true_f2\' was specified directly, you do not need to specify \'true_cor_XX\' or \'true_cor_YX\'.",
      arg = c("true_f2", "true_cor_XX", "true_cor_YX"))
    if (!is.null(true_R2_Y_X) || !is.null(true_R2_j_X_without_j) || !is.null(true_R2_Y_X_without_j)) .dmar_refuse_argument(
      "Since the effect size, \'true_f2\' was specified directly, you do not need to specify \'true_R2_Y_X\', \'true_R2_j_X_without_j\', or \'true_R2_Y_X_without_j\'.",
      arg = c("true_f2", "true_R2_Y_X", "true_R2_j_X_without_j"))
    if (true_f2 < 0) .dmar_refuse_argument(
      "The effect size, \'true_f2\', cannot be negative.", arg = "true_f2")
    f <- sqrt(true_f2)
  }

  # Confirm one of the effect size input methods actually defined 'f'. A NULL
  # 'f' means no valid effect size specification was supplied; a non-finite 'f'
  # means an inconsistent one (for example true_R2_Y_X_without_j > true_R2_Y_X, which
  # would otherwise produce NaNs partway through the search).
  if (is.null(f) || !is.finite(f)) {
    .dmar_refuse_argument(
      "No valid effect size specification was supplied. Provide one of: ('true_R2_Y_X', 'true_R2_Y_X_without_j'); ('true_cor_XX', 'true_cor_YX', 'which_predictor'); ('true_R2_Y_X', 'true_reg_coef_j', 'true_R2_j_X_without_j', 'true_sd_X', 'true_sd_Y'); or 'true_f2'. Also check that true_R2_Y_X_without_j does not exceed true_R2_Y_X.",
      arg = c("true_R2_Y_X", "true_R2_Y_X_without_j", "true_cor_XX"))
  }
  f <- as.numeric(f)
  # Every route but the matrix one (which sets 'p' from 'true_cor_XX') needs
  # 'p'. Without it the search below began at numeric(0).
  if (is.null(p)) .dmar_refuse_argument(c(
    "'p', the number of predictors in the full model, must be given when ",
    "the effect is given as squared multiple correlations, as 'true_reg_coef_j', ",
    "or as 'true_f2'; it is taken from 'true_cor_XX' only when ",
    "'true_cor_XX' and 'true_cor_YX' are given."), arg = "p")

  # Now that the noncentral t parameter has been defined, contine to sample size determination.
  #################################################

  if (is.null(specified_N)) {
    .power_at <- function(N_now) {
      .power_noncentral_t(ncp = sqrt(N_now) * abs(f), df = N_now - p - 1,
                          alpha_level = alpha_level,
                          directional = directional)
    }
    # Evaluate the smallest admissible size (df >= 1) before incrementing, so a
    # target the minimum already attains returns that minimum, not one above it.
    # Two-sided (nondirectional) power counts both rejection tails; see
    # .power_noncentral_t(). Summing only the upper tail returned alpha/2 at a
    # null effect for a nominal-alpha two-sided test.
    N_i <- p + 1 + 1
    # With an effect of exactly 0 the t statistic has the central t
    # distribution at every N, so the probability that the test rejects is
    # 'alpha_level' for the two-sided test (both tails) and for the
    # directional one (the upper tail, the direction .power_noncentral_t()
    # takes for the effect): the power never rises above 'alpha_level'. A
    # 'desired_power' above it is refused here, at once; through the first
    # builds of 1.1.0 the search ran to N = 10,000,000 (about 40 seconds)
    # and then stopped. A 'desired_power' at or below 'alpha_level' is met
    # at the minimum N, which the search below returns at its first step.
    # The comparison is with 'alpha_level' itself, the exact null rejection
    # probability, not with its computed value, which can differ from it in
    # the last digit.
    if (f == 0 && desired_power > alpha_level) .dmar_refuse(c(
      "The planning values give the predictor an effect of 0 (Cohen's f ",
      "squared is 0), so the null hypothesis is true: the probability ",
      "that the ", if (directional) "directional" else "two-sided",
      " test rejects is 'alpha_level' (", format(alpha_level),
      ") at every sample size and never reaches desired_power = ",
      format(desired_power), ". Give the nonzero effect the study is ",
      "planned to detect, or give 'specified_N' to have the power at ",
      "that sample size (which is 'alpha_level') reported."),
      "effect_zero", fn = "ss_power_reg_coef", alpha_level = alpha_level,
      desired_power = desired_power, directional = directional)
    # At an effect of 0 the power is 'alpha_level' exactly, reported as such.
    Actual_Power <- if (f == 0) alpha_level else .power_at(N_i)
    # Power rises with N for an effect that is not 0, so the power at the
    # largest size the search covers settles at once whether any size it
    # covers attains the target; the search runs only when one does. The
    # search covers sizes up to 10,000,000, or up to the minimum p + 2 when
    # that is larger (where the t statistic has its first degree of
    # freedom; at N = 10,000,000 a larger p would leave none).
    N_max <- max(1e7, N_i)
    if (Actual_Power < desired_power &&
        .power_at(N_max) < desired_power) .dmar_refuse(c(
      "ss_power_reg_coef() found no sample size up to ",
      format(N_max, big.mark = ",", scientific = FALSE), " at which ",
      "the power reaches desired_power = ", format(desired_power),
      "; the effect at the planning values (Cohen's f = ",
      format(signif(f, 3)), ") is too small for that. Check the planning ",
      "values, or ask for less power."),
      "computation_failed", fn = "ss_power_reg_coef",
      desired_power = desired_power)
    repeat {
      if (print_iter == TRUE) cat(c(Current.Power = Actual_Power, Current.NC.t.Parm = (sqrt(N_i) * abs(f)), Current.N = N_i), "\n")

      if (Actual_Power >= desired_power) break
      N_i <- N_i + 1
      if (N_i > 1e7) .dmar_refuse(
        "Failed to reach 'desired_power' within a reasonable sample size.",
        "computation_failed")
      Actual_Power <- .power_at(N_i)
    }
    VALUE <- data.frame(
      term = c("necessary_N", "actual_power", "noncentral_t_parm", "true_f"),
      value = c(N_i, Actual_Power, (sqrt(N_i) * abs(f)), f)
    )
  }

  if (!is.null(specified_N)) {
    # Reject a fractional N or one that leaves no residual degrees of freedom,
    # rather than returning NaN power.
    specified_N <- .check_whole_n(specified_N, "specified_N", p + 2L)
    actual_power_specified_N <- .power_noncentral_t(ncp = sqrt(specified_N) * abs(f),
                                                    df = specified_N - p - 1,
                                                    alpha_level = alpha_level,
                                                    directional = directional)
    VALUE <- data.frame(
      term = c("specified_N", "actual_power", "noncentral_t_parm", "true_f"),
      value = c(specified_N, actual_power_specified_N, (sqrt(specified_N) * abs(f)), f)
    )
  }


  class(VALUE) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  return(VALUE)
}

# The check of 'p' in ss_power_reg_coef() (and so in ss_power_rc()): a
# single finite whole number of at least 1, the number of predictors in
# the full model. NA, NaN, a string, and an empty vector are refused before
# this by .dmar_check_numbers(). Not exported.
.reg_coef_check_p <- function(p) {
  if (is.numeric(p) && length(p) == 1L && is.finite(p) && p == round(p) &&
      p >= 1) {
    return(invisible(p))
  }
  .dmar_refuse_argument(c(
    "'p' must be a single whole number of at least 1, the number of ",
    "predictors in the full model; ", .dmar_given(p), "."),
    arg = "p", value = p)
}

# The check of 'which_predictor' on the route through 'true_cor_XX' and
# 'true_cor_YX': a single whole number from 1 to p, the position of the
# predictor of interest. A fractional value was truncated by indexing, and
# one above p dropped no predictor, so that the effect was 0. Returns it as
# an integer. Not exported.
.reg_coef_check_which <- function(which_predictor, p) {
  if (is.numeric(which_predictor) && length(which_predictor) == 1L &&
      is.finite(which_predictor) &&
      which_predictor == round(which_predictor) &&
      which_predictor >= 1 && which_predictor <= p) {
    return(as.integer(which_predictor))
  }
  .dmar_refuse_argument(c(
    "'which_predictor' must be a single whole number from 1 to ", p,
    ", the position in 'true_cor_XX' of the predictor of interest; ",
    .dmar_given(which_predictor), "."),
    arg = "which_predictor", value = which_predictor)
}

# The check of 'true_cor_XX' as a correlation matrix: symmetric (to the
# tolerance of isSymmetric()), with 1 on the diagonal (to 1e-8), and
# positive definite (chol() succeeds). A symmetric matrix that can be
# inverted but is not positive definite (an off-diagonal 1.2) is not a
# correlation matrix, and through the first builds of 1.1.0 it gave R's
# "NaNs produced" warning and then the general refusal of an effect size
# specification. Not exported.
.reg_coef_check_rho <- function(true_cor_XX) {
  m <- unname(as.matrix(true_cor_XX))
  pd <- isSymmetric(m) && all(abs(diag(m) - 1) <= 1e-8) &&
    !is.null(tryCatch(chol(m), error = function(e) NULL))
  if (pd) return(invisible(true_cor_XX))
  .dmar_refuse_argument(c(
    "'true_cor_XX' is not a correlation matrix: it must be symmetric, ",
    "with 1 on the diagonal, and positive definite (no predictor a linear ",
    "combination of the others, and no set of correlations that cannot ",
    "occur together). Give a positive definite correlation matrix."),
    arg = "true_cor_XX")
}

# The squared multiple correlation of the outcome on the predictors whose
# correlation matrix is 'true_cor_XX' and whose correlations with the
# outcome are 'true_cor_YX': true_cor_YX' true_cor_XX^{-1} true_cor_YX. A
# matrix that cannot be inverted is not a correlation matrix of predictors
# that each add to the model, and is refused. Not exported.
.reg_coef_R2 <- function(true_cor_XX, true_cor_YX) {
  w <- tryCatch(solve(true_cor_XX, true_cor_YX), error = function(e) NULL)
  if (is.null(w)) .dmar_refuse_argument(c(
    "'true_cor_XX' cannot be inverted, so it is not the correlation matrix ",
    "of predictors that are each not a linear combination of the others. ",
    "Give a positive definite correlation matrix."), arg = "true_cor_XX")
  as.numeric(sum(true_cor_YX * w))
}
