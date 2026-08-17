#' Sample Size for a Targeted Regression Coefficient
#'
#' Determine the necessary sample size for a targeted regression coefficient or determine the degree of power given a specified sample size.
#'
#' @param rho2_Y_X Population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \emph{p} predictor variables (i.e., the \emph{X} variables)
#' @param rho2_Y_X_without_j Population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p}-1 predictor variables, where the one not used is the predictor of interest
#' @param p Number of predictor variables
#' @param desired_power Desired degree of statistical power for the test of targeted regression coefficient
#' @param alpha_level Type I error rate
#' @param directional Whether or not a direction or a nondirectional test is to be used (usually \code{directional=FALSE})
#' @param beta_j Population value of the regression coefficient for the predictor of interest
#' @param sigma_X Population standard deviation for the predictor variable of interest
#' @param sigma_Y Population standard deviation for the outcome variable
#' @param Rho2_j_X_without_j Population squared multiple correlation coefficient predicting the predictor variable of interest from the remaining \code{p}-1 predictor variables
#' @param rho_XX Population correlation matrix for the \emph{p} predictor variables
#' @param rho_YX Population vector of correlation coefficient between the \code{p} predictor variables and the criterion variable
#' @param which_predictor Identifies the predictor of interest when \code{rho_XX} and \code{rho_YX} are specified
#' @param cohen_f2 Cohen's (1988) definition for an effect size for a targeted regression coefficient: \code{(rho2_Y_X - rho2_Y_X_without_j) / (1 - rho2_Y_X)}
#' @param specified_N Sample size for which power should be evaluated. This is the \emph{total} sample size.
#' @param print_progress If the progress of the iterative procedure is printed to the screen as the iterations are occurring
#'
#' @details
#' Determines the necessary sample size given a desired level of statistical power. Alternatively, determines the statistical power for a given a specified sample size.
#' There are a number of ways that the specification regarding the size of the regression coefficient can be entered. The most basic, and often the simplest,
#' is to specify \code{rho2_Y_X} and \code{rho2_Y_X_without_j}. See the examples section for several options.
#'
#' @return A tidy \code{data.frame} with a \code{term} column and a numeric
#'   \code{value} column, forwarded unchanged from
#'   \code{\link{ss_power_reg_coef}}: rows for \code{necessary_N} (the
#'   necessary total sample size) or \code{specified_N} (when
#'   \code{specified_N} is supplied), \code{actual_power},
#'   \code{noncentral_t_parm} (the noncentrality of the \emph{t} distribution),
#'   and \code{effect_size} (the square root of \code{cohen_f2}, since
#'   \code{cohen_f2} is the effect size on the \emph{F} scale). The result
#'   carries the \code{dmar_ss_power} class, so \code{\link[generics]{tidy}}
#'   and \code{\link[generics]{glance}} summarize it in broom convention.
#'
#' @references
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
#'   planning for more accurate statistical power: A method adjusting
#'   sample effect sizes for publication bias and uncertainty.
#'   \emph{Psychological Science, 28}(11), 1547--1562.
#'   \doi{10.1177/0956797617723724}
#'
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
#' Maxwell, S. E. (2000). Sample size and multiple regression analysis. \emph{Psychological Methods, 5}(4), 434--458.
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
#'   c(0.66, 0.43, 0.50, 0.26, 0.30, 1.00))
#'
#' rho_XX <- Cor.Mat[2:6, 2:6]
#' rho_YX <- Cor.Mat[1, 2:6]
#'
#' # Method 1
#' ss_power_rc(rho2_Y_X = 0.7826786, rho2_Y_X_without_j = 0.7363697, p = 5,
#'             alpha_level = .05, directional = FALSE, desired_power = .80)
#'
#' # Method 2
#' ss_power_rc(alpha_level = .05, rho_XX = rho_XX, rho_YX = rho_YX, which_predictor = 5,
#'             directional = FALSE, desired_power = .80)
#'
#' # Method 3
#' # Here, beta_j is the standardized regression coefficient. Had beta_j
#' # been the unstandardized regression coefficient, sigma_X and sigma_Y
#' # would have been the standard deviation for the X variable of interest
#' # and Y, respectively.
#' ss_power_rc(rho2_Y_X = 0.7826786, Rho2_j_X_without_j = 0.3652136, beta_j = 0.2700964, p = 5,
#'             alpha_level = .05, sigma_X = 1, sigma_Y = 1,
#'             directional = FALSE, desired_power = .80)
#' # Method 4
#' ss_power_rc(alpha_level = .05, cohen_f2 = 0.2130898, p = 5,
#'             directional = FALSE, desired_power = .80)
#' # Power given a specified N and squared multiple correlation coefficients.
#' ss_power_rc(rho2_Y_X = 0.7826786, rho2_Y_X_without_j = 0.7363697, specified_N = 25, p = 5,
#'             alpha_level = .05, directional = FALSE)
#'
#' # Power given a specified N and effect size.
#' ss_power_rc(alpha_level = .05, cohen_f2 = 0.2130898, p = 5, specified_N = 25, directional = FALSE)
#'
#' # Reproducing Maxwell's (2000, p. 445) Example
#' Cor.Mat.Maxwell <- rbind(
#'     c(1.00, 0.35,  0.20, 0.20, 0.20, 0.20),
#'     c(0.35, 1.00,  0.40, 0.40, 0.40, 0.40),
#'     c(0.20, 0.40,  1.00, 0.45, 0.45, 0.45),
#'     c(0.20, 0.40,  0.45, 1.00, 0.45, 0.45),
#'     c(0.20, 0.40,  0.45, 0.45, 1.00, 0.45),
#'     c(0.20, 0.40,  0.45, 0.45, 0.45, 1.00)
#' )
#'
#' RHO.XX.Maxwell <- Cor.Mat.Maxwell[2:6, 2:6]
#' Rho.YX.Maxwell <- Cor.Mat.Maxwell[1, 2:6]
#' R2.Maxwell <- Rho.YX.Maxwell %*% solve(RHO.XX.Maxwell) %*% Rho.YX.Maxwell
#'
#' RHO.XX.Maxwell.no.1 <- Cor.Mat.Maxwell[3:6, 3:6]
#' Rho.YX.Maxwell.no.1 <- Cor.Mat.Maxwell[1, 3:6]
#' R2.Maxwell.no.1 <- Rho.YX.Maxwell.no.1 %*% solve(RHO.XX.Maxwell.no.1) %*% Rho.YX.Maxwell.no.1
#'
#'
#' # Note that Maxwell arrives at N=113, whereas this procedure arrives at 111.
#' # This seems to be the case becuase of rounding error in calculations
#' # and tables (Cohen, 1988) used. The present procedure is correct and
#' # contains no rounding error in the application of the method.
#' ss_power_rc(rho2_Y_X = R2.Maxwell, rho2_Y_X_without_j = R2.Maxwell.no.1, p = 5,
#'             alpha_level = .05, directional = FALSE, desired_power = .80)
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


ss_power_rc <- function(rho2_Y_X = NULL, rho2_Y_X_without_j = NULL, p = NULL,
                        desired_power = .85, alpha_level = .05, directional = FALSE, beta_j = NULL,
                        sigma_X = NULL, sigma_Y = NULL, Rho2_j_X_without_j = NULL, rho_XX = NULL, rho_YX = NULL,
                        which_predictor = NULL, cohen_f2 = NULL, specified_N = NULL, print_progress = FALSE) {
  ss_power_reg_coef(
    rho2_Y_X = rho2_Y_X, rho2_Y_X_without_j = rho2_Y_X_without_j, p = p,
    desired_power = desired_power, alpha_level = alpha_level, directional = directional, beta_j = beta_j,
    sigma_X = sigma_X, sigma_Y = sigma_Y, rho2_j_X_without_j = Rho2_j_X_without_j, rho_XX = rho_XX, rho_YX = rho_YX,
    which_predictor = which_predictor, cohen_f2 = cohen_f2, specified_N = specified_N, print_progress = print_progress
  )
}
