#' Sensitivity Analysis for the Power of a Targeted Regression Coefficient
#'
#' @description
#' Monte Carlo sensitivity analysis for the statistical power of the
#' \emph{t}-test of a targeted regression coefficient. Given a planned
#' (\code{estimated_*}) covariance structure and a true (\code{true_*})
#' covariance structure, the function draws \code{G} replications,
#' fits the multiple regression, and reports the empirical proportion of
#' replications on which the \emph{t}-test of the targeted coefficient
#' rejects, together with the realized distribution of \eqn{\hat b_j},
#' its standard error, and the test statistic. \code{ss_power_reg_coef_sensitivity()}
#' is the power-oriented sibling of \code{\link{ss_aipe_reg_coef_sensitivity}}
#' (which is CI-width oriented).
#'
#' @param true_var_Y Population variance of the dependent variable (\emph{Y})
#' @param true_cov_YX Population covariance vector between the \code{p} predictor variables and the dependent variable (\emph{Y})
#' @param true_cov_XX Population covariance matrix of the \code{p} predictor variables
#' @param estimated_var_Y Estimated variance of the dependent variable (\emph{Y}) used in sample size planning. Defaults to \code{true_var_Y}.
#' @param estimated_cov_YX Estimated covariance vector between the predictor variables and the dependent variable used in sample size planning. Defaults to \code{true_cov_YX}.
#' @param estimated_cov_XX Estimated covariance matrix of the predictor variables used in sample size planning. Defaults to \code{true_cov_XX}.
#' @param specified_N Directly specified sample size; if supplied, sample size planning is skipped.
#' @param which_predictor Index identifying which of the \code{p} predictors is the targeted predictor for the power test.
#' @param desired_power Desired degree of statistical power used for planning
#' @param alpha_level Type I error rate
#' @param directional Whether a one-sided or two-sided test is used
#' @param standardize Whether each replication's data should be standardized prior to fitting (giving a standardized regression coefficient)
#' @param G Number of Monte Carlo replications
#' @param print_iter Whether to print the iteration number during the simulation
#' @param save Whether to write the per-replication results to a CSV file
#' @param filename Name of the CSV file written when \code{save = TRUE}
#'
#' @details
#' When the estimated and true covariance structures are identical, the
#' function performs a Monte Carlo evaluation of the planning procedure
#' (no misspecification); when they differ, it performs a sensitivity
#' analysis on the consequences of misspecifying the population
#' covariance structure for the targeted coefficient's power. The planning
#' step calls \code{\link{ss_power_reg_coef}} with the estimated
#' covariance structure; the simulation step generates data from the true
#' covariance structure.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis. The \code{term}
#' entries are: \code{total_N} (the sample size evaluated),
#' \code{empirical_power}, \code{analytic_power} (computed from
#' \code{\link{ss_power_reg_coef}}), the mean / median / SD of the
#' realized \eqn{\hat b_j} (\code{mean_b_j}, \code{median_b_j},
#' \code{sd_b_j}), of its standard error (\code{mean_se_b_j},
#' \code{median_se_b_j}, \code{sd_se_b_j}), of the test statistic
#' (\code{mean_t}, \code{median_t}, \code{sd_t}), and of the squared
#' multiple correlation coefficient (\code{mean_R2}, \code{median_R2},
#' \code{sd_R2}), \code{t_crit} (the critical value), and the input
#' echoes \code{p}, \code{which_predictor}, \code{true_b_j} and
#' \code{estimated_b_j} (the population and planning values of the
#' targeted coefficient implied by the supplied covariance structures),
#' \code{desired_power} (NA when \code{specified_N} was supplied
#' instead), and \code{alpha_level}. The result carries the
#' \code{dmar_ss_power_sensitivity} class, so \code{\link[generics]{tidy}}
#' reports the planned sample size beside the empirical and analytic power,
#' and \code{\link[generics]{glance}} adds the simulated estimator
#' distribution beside the echoed inputs.
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
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons of
#'   means and Chapter 6 on trend analysis.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_power_reg_coef}}, \code{\link{ss_aipe_reg_coef_sensitivity}}
#'
#' @examples
#' # Targeted coefficient power sensitivity with two predictors. The
#' # default G = 1000 replications is used in practice; G is reduced here
#' # so the example runs quickly.
#' set.seed(113)
#' Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' cov_YX  <- c(0.4, 0.3)
#' ss_power_reg_coef_sensitivity(
#'   true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
#'   which_predictor = 1, desired_power = 0.80,
#'   G = 100, print_iter = FALSE
#' )
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
#' @importFrom MASS mvrnorm
#' @importFrom utils write.csv

ss_power_reg_coef_sensitivity <- function(true_var_Y = NULL, true_cov_YX = NULL, true_cov_XX = NULL,
                                          estimated_var_Y = NULL, estimated_cov_YX = NULL,
                                          estimated_cov_XX = NULL, specified_N = NULL,
                                          which_predictor = 1, desired_power = .85,
                                          alpha_level = .05, directional = FALSE,
                                          standardize = FALSE, G = 1000,
                                          print_iter = TRUE, save = FALSE,
                                          filename = "ss_power_reg_coef_sensitivity_result.csv") {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("The package \'MASS\' is needed; please install it.")

  if (is.null(true_cov_XX)) stop("You must specify \'true_cov_XX\' (the covariance matrix of the predictors).")
  if (is.null(true_cov_YX)) stop("You must specify \'true_cov_YX\' (the covariance vector of the predictors with Y).")
  if (is.null(true_var_Y))  true_var_Y <- 1

  if ((sum(round(true_cov_XX, 5) == round(t(true_cov_XX), 5))) !=
      (dim(true_cov_XX)[1] * dim(true_cov_XX)[2]))
    stop("\'true_cov_XX\' must be symmetric.")

  p <- dim(true_cov_XX)[1]

  if (!is.numeric(which_predictor) || length(which_predictor) != 1L ||
      which_predictor < 1 || which_predictor > p ||
      which_predictor != as.integer(which_predictor)) {
    stop(sprintf("\'which_predictor\' must be an integer between 1 and p = %d.", p))
  }
  if (is.null(specified_N) && (desired_power <= 0 || desired_power >= 1)) {
    stop("\'desired_power\' must be in (0, 1).")
  }
  if (alpha_level <= 0 || alpha_level >= 1)
    stop("\'alpha_level\' must be in (0, 1).")
  if (!is.numeric(G) || length(G) != 1L || G < 1 || G != as.integer(G))
    stop("\'G\' must be a positive integer.")

  if (is.null(estimated_var_Y))  estimated_var_Y  <- true_var_Y
  if (is.null(estimated_cov_XX)) estimated_cov_XX <- true_cov_XX
  if (is.null(estimated_cov_YX)) estimated_cov_YX <- true_cov_YX

  if ((sum(round(estimated_cov_XX, 5) == round(t(estimated_cov_XX), 5))) !=
      (dim(estimated_cov_XX)[1] * dim(estimated_cov_XX)[2]))
    stop("\'estimated_cov_XX\' must be symmetric.")

  Estimated_Sigma <- cbind(c(estimated_var_Y, estimated_cov_YX),
                           rbind(estimated_cov_YX, estimated_cov_XX))
  True_Sigma <- cbind(c(true_var_Y, true_cov_YX),
                      rbind(true_cov_YX, true_cov_XX))

  Est_Rho2_Y_X  <- (estimated_cov_YX %*% solve(estimated_cov_XX) %*% estimated_cov_YX) / estimated_var_Y
  True_Rho2_Y_X <- (true_cov_YX      %*% solve(true_cov_XX)      %*% true_cov_YX)      / true_var_Y

  if (Est_Rho2_Y_X  > 1) stop("The estimated covariance structure implies R^2 > 1.")
  if (True_Rho2_Y_X > 1) stop("The true covariance structure implies R^2 > 1.")

  True_b_j <- (solve(true_cov_XX) %*% true_cov_YX)[which_predictor]

  if (is.null(specified_N)) {
    Est_cor <- cov2cor(Estimated_Sigma)
    N <- ss_power_reg_coef(
      rho_XX = Est_cor[2:(p + 1), 2:(p + 1)],
      rho_YX = Est_cor[1, 2:(p + 1)],
      which_predictor = which_predictor,
      desired_power = desired_power,
      alpha_level = alpha_level,
      directional = directional
    )$value[1]
  } else {
    # A fractional N, or one leaving no residual degrees of freedom, would give
    # a degenerate simulated fit (NaN power); reject it.
    N <- .check_whole_n(specified_N, "specified_N", p + 2L)
  }

  MU <- rep(0, p + 1)
  Results <- matrix(NA_real_, G, 4)
  colnames(Results) <- c("b_j", "se_b_j", "t_stat", "r_2")

  t_crit <- if (directional) {
    qt(1 - alpha_level, df = N - p - 1)
  } else {
    qt(1 - alpha_level / 2, df = N - p - 1)
  }
  # For directional tests with True_b_j == 0 the rejection rule is
  # degenerate; fall back to the two-sided critical value so the empirical
  # Type I rate equals alpha_level.
  t_crit_two_sided <- qt(1 - alpha_level / 2, df = N - p - 1)

  for (i in seq_len(G)) {
    if (print_iter) cat(i, "\n")
    DATA <- MASS::mvrnorm(N, mu = MU, Sigma = True_Sigma)
    if (standardize) DATA <- scale(DATA)
    fit  <- lm(DATA[, 1] ~ DATA[, -1, drop = FALSE])
    sm   <- summary(fit)
    co   <- coef(sm)[which_predictor + 1, ]
    Results[i, 1] <- co[1]
    Results[i, 2] <- co[2]
    Results[i, 3] <- co[1] / co[2]
    Results[i, 4] <- sm$r.squared
  }

  Results_df <- as.data.frame(Results)
  if (save) {
    message("Simulation results will be saved to a .csv file; overwriting an existing file of the same name.")
    utils::write.csv(Results_df, filename, row.names = FALSE)
  }

  # Rejection rule. Directional tests reject in the direction implied by
  # the sign of True_b_j (the standard convention: ss_power_reg_coef itself
  # treats the effect size as |f| and computes power for the upper tail of
  # the noncentral t, so the test direction is "the assumed direction of
  # the true effect"). With True_b_j exactly 0 there is no direction, so
  # the rule falls back to the two-sided critical value.
  rej <- if (directional) {
    if (True_b_j > 0) {
      Results_df$t_stat >  t_crit
    } else if (True_b_j < 0) {
      Results_df$t_stat < -t_crit
    } else {
      abs(Results_df$t_stat) > t_crit_two_sided
    }
  } else {
    abs(Results_df$t_stat) > t_crit
  }
  empirical_power <- mean(rej, na.rm = TRUE)

  analytic_power <- ss_power_reg_coef(
    rho_XX = cov2cor(True_Sigma)[2:(p + 1), 2:(p + 1)],
    rho_YX = cov2cor(True_Sigma)[1, 2:(p + 1)],
    which_predictor = which_predictor,
    specified_N = N,
    alpha_level = alpha_level,
    directional = directional
  )$value[2]

  Estimated_b_j <- (solve(estimated_cov_XX) %*% estimated_cov_YX)[which_predictor]

  Summary <- data.frame(
    term = c("total_N", "empirical_power", "analytic_power",
             "mean_b_j", "median_b_j", "sd_b_j",
             "mean_se_b_j", "median_se_b_j", "sd_se_b_j",
             "mean_t", "median_t", "sd_t",
             "mean_R2", "median_R2", "sd_R2", "t_crit",
             "p", "which_predictor",
             "true_b_j", "estimated_b_j",
             "desired_power", "alpha_level"),
    value = c(N, empirical_power, analytic_power,
              mean(Results_df$b_j, na.rm = TRUE),
              median(Results_df$b_j, na.rm = TRUE),
              sqrt(var(Results_df$b_j, na.rm = TRUE)),
              mean(Results_df$se_b_j, na.rm = TRUE),
              median(Results_df$se_b_j, na.rm = TRUE),
              sqrt(var(Results_df$se_b_j, na.rm = TRUE)),
              mean(Results_df$t_stat, na.rm = TRUE),
              median(Results_df$t_stat, na.rm = TRUE),
              sqrt(var(Results_df$t_stat, na.rm = TRUE)),
              mean(Results_df$r_2, na.rm = TRUE),
              median(Results_df$r_2, na.rm = TRUE),
              sqrt(var(Results_df$r_2, na.rm = TRUE)),
              t_crit,
              p, which_predictor,
              True_b_j, Estimated_b_j,
              if (is.null(specified_N)) desired_power else NA_real_,
              alpha_level)
  )

  class(Summary) <- c("dmar_ss_power_sensitivity", class(Summary))
  return(.as_dmar_tbl(Summary))
}
