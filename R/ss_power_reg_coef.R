#' Sample Size for a Targeted Regression Coefficient
#'
#' Determine the necessary sample size for a targeted regression coefficient or determine the degree of power given a specified sample size
#'
#' @param rho2_Y_X Population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p} predictor variables (i.e., the \emph{X} variables)
#' @param rho2_Y_X_without_j Population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p}-1 predictor variables, where the one not used is the predictor of interest
#' @param p Number of predictor variables
#' @param desired_power Desired degree of statistical power for the test of targeted regression coefficient
#' @param alpha_level Type I error rate
#' @param directional Whether or not a direction or a nondirectional test is to be used (usually \code{directional=FALSE})
#' @param beta_j Population value of the regression coefficient for the predictor of interest
#' @param sigma_X Population standard deviation for the predictor variable of interest
#' @param sigma_Y Population standard deviation for the outcome variable
#' @param rho2_j_X_without_j Population squared multiple correlation coefficient predicting the predictor variable of interest from the remaining p-1 predictor variables
#' @param rho_XX Population correlation matrix for the \code{p} predictor variables
#' @param rho_YX Population vector of correlation coefficient between the \code{p} predictor variables and the criterion variable
#' @param cohen_f2 Cohen's (1988) definition for an effect size for a targeted regression coefficient: \code{(rho2_Y_X-rho2_Y_X_without_j)/(1-rho2_Y_X)}
#' @param which_predictor Identifies the predictor of interest when \code{rho_XX} and \code{rho_YX} are specified
#' @param specified_N Sample size for which power should be evaluated. This is the \emph{total} sample size.
#' @param print_progress If the progress of the iterative procedure is printed to the screen as the iterations are occurring
#'
#' @details
#' Determines the necessary sample size given a desired level of statistical power. Alternatively,
#' determines the statistical power for a given a specified sample size. There are a number of ways
#' that the specification regarding the size of the regression coefficient can be entered. The most
#' basic, and often the simplest, is to specify \code{rho2_Y_X} and \code{rho2_Y_X_without_j}.
#' See the examples section for several options.
#'
#' Power is computed from a noncentral \emph{t} distribution with noncentrality
#' \eqn{\sqrt{N}\,f}, which treats the predictors as fixed (their values held
#' constant across hypothetical replications). This is the standard
#' fixed-predictor power analysis; under random predictors, where the predictor
#' values themselves vary from sample to sample, the sample size required for a
#' given level of power is somewhat larger.
#'
#' @return
#' \item{ss}{Either the necessary sample size or the specified sample size, depending if one is interested in determining the necessary sample size given a desired degree of statistical power or if one is interested in the determining the value of statistical power given a specified sample size, respectively}
#' \item{actual_power}{Actual power of the situation described}
#' \item{noncentral_t_parm}{Value of the noncentral distribution for the appropriate \emph{t}-distribution}
#' \item{effect_size}{Effect size for the noncentral \emph{t}-distribution; this is the square root of \code{cohen_f2}, because \code{cohen_f2} is the effect size using an \emph{F}-distribution}
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
#' \code{\link{ss_aipe_reg_coef}}, \code{\link{ss_power_R2}}, \code{\link{conf_limits_ncf}}
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
#' rho_XX <- Cor.Mat[2:6, 2:6]
#' rho_YX <- Cor.Mat[1, 2:6]
#'
#' # Method 1
#' ss_power_reg_coef(rho2_Y_X = 0.7826786, rho2_Y_X_without_j = 0.7363697, p = 5,
#'                   alpha_level = .05, directional = FALSE, desired_power = .80)
#'
#' # Method 2
#' ss_power_reg_coef(alpha_level = .05, rho_XX = rho_XX, rho_YX = rho_YX, which_predictor = 5,
#'                   directional = FALSE, desired_power = .80)
#'
#' # Method 3
#' # Here, beta_j is the standardized regression coefficient. Had beta_j
#' # been the unstandardized regression coefficient, sigma_X and sigma_Y
#' # would have been the standard deviation for the X variable of
#' # interest and Y, respectively.
#' ss_power_reg_coef(rho2_Y_X = 0.7826786, rho2_j_X_without_j = 0.3652136, beta_j = 0.2700964,
#'                   p = 5, alpha_level = .05, sigma_X = 1, sigma_Y = 1, directional = FALSE,
#'                   desired_power = .80)
#'
#' # Method 4
#' ss_power_reg_coef(alpha_level = .05, cohen_f2 = 0.2130898, p = 5,
#'                   directional = FALSE, desired_power = .80)
#'
#' # Power given a specified N and squared multiple correlation coefficients.
#' ss_power_reg_coef(rho2_Y_X = 0.7826786, rho2_Y_X_without_j = 0.7363697, specified_N = 25,
#'                   p = 5, alpha_level = .05, directional = FALSE)
#'
#' # Power given a specified N and effect size.
#' ss_power_reg_coef(alpha_level = .05, cohen_f2 = 0.2130898, p = 5, specified_N = 25,
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
#' ss_power_reg_coef(rho2_Y_X = R2.Maxwell, rho2_Y_X_without_j = R2.Maxwell.no.1, p = 5,
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

ss_power_reg_coef <- function(rho2_Y_X = NULL, rho2_Y_X_without_j = NULL, p = NULL, desired_power = .85, alpha_level = .05, directional = FALSE, beta_j = NULL, sigma_X = NULL, sigma_Y = NULL, rho2_j_X_without_j = NULL, rho_XX = NULL, rho_YX = NULL, which_predictor = NULL, cohen_f2 = NULL, specified_N = NULL, print_progress = FALSE) {

  # Validate 'desired_power' before the search: power tops out below 1, so a
  # 'desired_power' at or above 1 would drive the 'while (Dif > 0)' loop below
  # forever. Only enforced when a sample size is being sought.
  if (is.null(specified_N) && (!is.numeric(desired_power) || length(desired_power) != 1L || desired_power <= 0 || desired_power >= 1)) {
    stop("'desired_power' must be a single numeric value in (0, 1).", call. = FALSE)
  }

  # 'f' is the noncentrality-defining effect size; each input method below sets
  # it. Initialize to NULL so a missing or inconsistent effect size
  # specification fails loudly rather than surfacing as 'object f not found'.
  f <- NULL

  # Method 1 to define the effect size.
  if (!is.null(rho2_Y_X) && !is.null(rho2_Y_X_without_j)) {
    if (!is.null(rho_XX) || !is.null(rho_YX)) stop("Since \'rho2_Y_X,\' and \'rho2_Y_X_without_j\', have been specified, do not specify \'rho_XX\' or \'rho_YX\'.")
    if (!is.null(which_predictor)) stop("Since you have specified \'rho2_Y_X\' and \'rho2_Y_X_without_j\', do not specify \'which_predictor\'.")
    f <- sqrt((rho2_Y_X - rho2_Y_X_without_j) / (1 - rho2_Y_X))
  }


  # Method 2 to define the effect size.
  if (is.null(rho2_Y_X) && is.null(rho2_j_X_without_j) && is.null(rho2_Y_X_without_j) && !is.null(rho_YX) && !is.null(rho_XX)) {
    if (is.null(which_predictor)) stop("Since you have specified \'rho_XX\' and \'rho_YX\', you must also specify the \'which_predictor\'.")
    if (!is.null(p)) stop("There is no need to specify \'p\' in this situation.")
    if ((dim(rho_XX)[1] != dim(rho_XX)[2]) || (dim(rho_XX)[1] != length(rho_YX))) stop("There is a problem with \'rho_XX\' and/or \'rho_YX\'")
    p <- dim(rho_XX)[1]
    rho2_Y_X <- (rho_YX %*% solve(rho_XX) %*% rho_YX)
    rho2_Y_X_without_j <- (rho_YX[-which_predictor] %*% solve(rho_XX[-which_predictor, -which_predictor]) %*% rho_YX[-which_predictor])
    f <- sqrt((rho2_Y_X - rho2_Y_X_without_j) / (1 - rho2_Y_X))
  }

  # Method 3 to define the effect size.
  if (!is.null(rho2_Y_X) && !is.null(beta_j) && !is.null(rho2_j_X_without_j)) {
    if (!is.null(which_predictor)) stop("Since you have specified \'rho2_Y_X\', \'beta_j\', and \'rho2_j_X_without_j\', do not specify \'which_predictor\'.")
    if (is.null(sigma_X) || is.null(sigma_Y)) stop("Since you have specified \'rho2_Y_X\', \'beta_j\', and \'rho2_j_X_without_j\', you must also specify \'sigma_X\' and \'sigma_Y\'.")
    f <- beta_j * sqrt((1 - rho2_j_X_without_j) / (1 - rho2_Y_X)) * (sigma_X / sigma_Y)
  }

  # Method 4, input the effect size, a la Cohen, directly.
  if (!is.null(cohen_f2)) {
    if (!is.null(which_predictor)) stop("You do not need to specify \'which_predictor\'.")
    if (!is.null(sigma_X) || !is.null(sigma_Y)) stop("You do not need to specify \'sigma_X\' or \'sigma_Y\'.")
    if (is.null(p)) stop("You need to specify \'p\'.")
    if (!is.null(rho_XX) || !is.null(rho_YX)) stop("Since the effect size, \'cohen_f2\' was specified directly, you do not need to specify \'rho_XX\' or \'rho_YX\'.")
    if (!is.null(rho2_Y_X) || !is.null(rho2_j_X_without_j) || !is.null(rho2_Y_X_without_j)) stop("Since the effect size, \'cohen_f2\' was specified directly, you do not need to specify \'rho2_Y_X\', \'rho2_j_X_without_j\', or \'rho2_Y_X_without_j\'.")
    if (cohen_f2 < 0) stop("The effect size, \'cohen_f2\', cannot be negative.")
    f <- sqrt(cohen_f2)
  }

  # Confirm one of the effect size input methods actually defined 'f'. A NULL
  # 'f' means no valid effect size specification was supplied; a non-finite 'f'
  # means an inconsistent one (for example rho2_Y_X_without_j > rho2_Y_X, which
  # would otherwise produce NaNs partway through the search).
  if (is.null(f) || !is.finite(f)) {
    stop("No valid effect size specification was supplied. Provide one of: ('rho2_Y_X', 'rho2_Y_X_without_j'); ('rho_XX', 'rho_YX', 'which_predictor'); ('rho2_Y_X', 'beta_j', 'rho2_j_X_without_j', 'sigma_X', 'sigma_Y'); or 'cohen_f2'. Also check that rho2_Y_X_without_j does not exceed rho2_Y_X.", call. = FALSE)
  }

  # Now that the noncentral t parameter has been defined, contine to sample size determination.
  #################################################

  if (is.null(specified_N)) {
    # Evaluate the smallest admissible size (df >= 1) before incrementing, so a
    # target the minimum already attains returns that minimum, not one above it.
    # Two-sided (nondirectional) power counts both rejection tails; see
    # .power_noncentral_t(). Summing only the upper tail returned alpha/2 at a
    # null effect for a nominal-alpha two-sided test.
    N_i <- p + 1 + 1
    repeat {
      Actual_Power <- .power_noncentral_t(ncp = sqrt(N_i) * abs(f),
                                          df = N_i - p - 1,
                                          alpha_level = alpha_level,
                                          directional = directional)

      if (print_progress == TRUE) cat(c(Current.Power = Actual_Power, Current.NC.t.Parm = (sqrt(N_i) * abs(f)), Current.N = N_i), "\n")

      if (Actual_Power >= desired_power) break
      N_i <- N_i + 1
      if (N_i > 1e7) stop("Failed to reach 'desired_power' within a reasonable sample size.", call. = FALSE)
    }
    VALUE <- data.frame(
      term = c("necessary_N", "actual_power", "noncentral_t_parm", "effect_size"),
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
      term = c("specified_N", "actual_power", "noncentral_t_parm", "effect_size"),
      value = c(specified_N, actual_power_specified_N, (sqrt(specified_N) * abs(f)), f)
    )
  }


  class(VALUE) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  return(VALUE)
}
