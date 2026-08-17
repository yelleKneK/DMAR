#' Sensitivity Analysis for Sample Size Planning to Make the Omnibus Test of \eqn{R^2} Sufficiently Powerful
#'
#' @description
#' Monte Carlo sensitivity analysis for the power of the omnibus
#' \emph{F}-test of the squared multiple correlation coefficient (\eqn{R^2}).
#' Given an \code{estimated_R2} used for sample size planning and a
#' \code{true_R2} that actually obtains in the population (the two need not
#' agree), the function draws \code{G} replications, fits the regression,
#' compares \eqn{F} to its critical value, and reports the realized empirical
#' power and a summary of the realized \eqn{R^2} and \eqn{F} distributions.
#' The simulation honors the same \code{random_predictors} /
#' \code{generate_random_predictors} crossing as
#' \code{\link{ss_aipe_R2_sensitivity}}, so the user can examine the effect
#' of planning under one regression model (fixed or random predictors) but
#' actually realizing the other.
#'
#' @param true_R2 Value of the population squared multiple correlation coefficient
#' @param estimated_R2 Value of the squared multiple correlation coefficient used for sample size planning. Either \code{estimated_R2} or \code{specified_N} must be supplied (not both).
#' @param desired_power Desired degree of statistical power used for planning
#' @param p Number of predictors
#' @param alpha_level Type I error rate
#' @param random_predictors Whether the sample size planning step treats predictors as random (\code{TRUE}, the default) or fixed (\code{FALSE})
#' @param specified_N Sample size at which the realized power should be computed; alternative to specifying \code{estimated_R2}
#' @param generate_random_predictors Whether the internal simulation should generate predictors as random (\code{TRUE}, the default) or fixed (\code{FALSE})
#' @param rho_yx Correlation between the dependent variable (\emph{Y}) and each of the \emph{X} variables
#' @param rho_xx Correlation among the \emph{X} variables (off-diagonal of the predictor correlation matrix)
#' @param G Number of Monte Carlo replications
#' @param print_iter Whether to print the iteration number during the simulation
#' @param save Whether to write the per-replication results to a CSV file
#' @param filename Name of the CSV file written when \code{save = TRUE}
#' @param \dots Additional arguments forwarded to internal helpers
#'
#' @details
#' When \code{estimated_R2} equals \code{true_R2}, the function performs a
#' straight Monte Carlo evaluation of the planning procedure (no
#' misspecification). Pass \code{specified_N} to evaluate realized power at
#' a specified sample size; in that case \code{estimated_R2} must not be
#' supplied. The crossing of \code{random_predictors} (used in planning)
#' with \code{generate_random_predictors} (used in the simulation) lets
#' the user inspect the consequences of planning under one regression
#' model but realizing the other. See Gatsonis and Sampson (1989) for
#' the comparison of fixed and random predictor power for the omnibus
#' test.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries are: \code{total_N} (the
#' sample size evaluated), \code{empirical_power} (the proportion of
#' replications on which \eqn{F} exceeded the critical value),
#' \code{analytic_power} (computed from \code{\link{ss_power_R2}} under
#' the same model as planning), \code{mean_R2} / \code{median_R2} /
#' \code{sd_R2} and \code{mean_F} / \code{median_F} / \code{sd_F}
#' (summaries of the realized \eqn{R^2} and \eqn{F}), \code{F_crit}
#' (the critical value), and the input echoes \code{p},
#' \code{true_R2}, \code{estimated_R2} and \code{desired_power} (both
#' NA when \code{specified_N} was supplied instead), and
#' \code{alpha_level}. The result carries the
#' \code{dmar_ss_power_sensitivity} class, so \code{\link[generics]{tidy}}
#' reports the planned sample size beside the empirical and analytic power,
#' and \code{\link[generics]{glance}} adds the simulated \eqn{R^2} and
#' \eqn{F} distribution beside the echoed inputs.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact
#'   power and sample size calculations. \emph{Psychological Bulletin,
#'   106}(3), 516--524.
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Lee, Y. S. (1971). Some results on the sampling distribution of the
#'   multiple correlation coefficient. \emph{Journal of the Royal
#'   Statistical Society, Series B, 33}(1), 117--130.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_power_R2}}, \code{\link{ss_aipe_R2_sensitivity}},
#' \code{\link{ci_R2}}
#'
#' @examples
#' set.seed(113)
#' # Realized power when planning under the fixed-predictor model but the
#' # data are actually generated with random predictors. G is small here
#' # for illustration; use G = 10,000 in practice.
#' ss_power_R2_sensitivity(true_R2 = 0.30, estimated_R2 = 0.30,
#'                         desired_power = 0.80, p = 5,
#'                         random_predictors = FALSE,
#'                         generate_random_predictors = TRUE,
#'                         G = 200, print_iter = FALSE)
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

ss_power_R2_sensitivity <- function(true_R2 = NULL, estimated_R2 = NULL,
                                    desired_power = .85, p = NULL,
                                    alpha_level = .05,
                                    random_predictors = TRUE,
                                    specified_N = NULL,
                                    generate_random_predictors = TRUE,
                                    rho_yx = .3, rho_xx = .3, G = 10000,
                                    print_iter = TRUE, save = FALSE,
                                    filename = "ss_power_r2_sensitivity_result.csv",
                                    ...) {
  if (is.null(true_R2)) stop("You must specify \'true_R2\'.")
  if (true_R2 >= 1 || true_R2 <= 0) stop("\'true_R2\' must be between zero and one.")
  if (is.null(p)) stop("You must specify \'p\' (the number of predictors).")
  if (alpha_level <= 0 || alpha_level >= 1) stop("\'alpha_level\' must be in (0, 1).")
  if (is.null(specified_N) && (desired_power <= 0 || desired_power >= 1)) {
    stop("\'desired_power\' must be in (0, 1).")
  }
  if (!is.numeric(G) || length(G) != 1L || G < 1 || G != as.integer(G))
    stop("\'G\' must be a positive integer.")

  if (is.null(estimated_R2) && is.null(specified_N))
    stop("You must specify either \'estimated_R2\' or \'specified_N\'.", call. = FALSE)
  if (!is.null(estimated_R2) && !is.null(specified_N))
    stop("You must specify either \'estimated_R2\' or \'specified_N\', not both.", call. = FALSE)

  prev_warn <- getOption("warn")
  on.exit(options(warn = prev_warn), add = TRUE)
  options(warn = -1)

  if (!is.null(estimated_R2)) {
    if (estimated_R2 >= 1 || estimated_R2 <= 0) stop("\'estimated_R2\' must be between zero and one.")
    N <- ss_power_R2(population_R2 = estimated_R2,
                     alpha_level = alpha_level,
                     desired_power = desired_power,
                     p = p,
                     random_predictors = random_predictors)$value[1]
  } else {
    # A fractional N, or one leaving no residual degrees of freedom, would give
    # a degenerate simulated fit (NaN power, saturated R^2); reject it.
    N <- .check_whole_n(specified_N, "specified_N", p + 2L)
  }

  # Build the joint Y, X covariance so that the population R^2 is exactly
  # true_R2, with the specified rho_yx and rho_xx.
  MU       <- rep(0, p + 1)
  sigma_YX <- rbind(rep(rho_yx, p))
  sigma_XY <- t(sigma_YX)
  Sigma_XX <- matrix(rep(rho_xx, p^2), nrow = p, ncol = p)
  diag(Sigma_XX) <- 1
  Numerator_P_Square <- sigma_YX %*% solve(Sigma_XX) %*% sigma_XY
  sigma_Y  <- Numerator_P_Square / true_R2
  Sigma    <- rbind(c(sigma_Y, sigma_YX), cbind(sigma_XY, Sigma_XX))

  f_crit <- qf(1 - alpha_level, df1 = p, df2 = N - p - 1)

  Results <- matrix(NA_real_, G, 2)
  colnames(Results) <- c("r2", "f_stat")

  if (generate_random_predictors) {
    for (i in seq_len(G)) {
      if (print_iter) cat(i, "\n")
      DATA <- MASS::mvrnorm(N, mu = MU, Sigma = Sigma)
      fit  <- lm(DATA[, 1] ~ DATA[, -1, drop = FALSE])
      sm   <- summary(fit)
      Results[i, 1] <- sm$r.squared
      Results[i, 2] <- sm$fstatistic[1]
    }
  } else {
    DATA_Pop <- MASS::mvrnorm(N, mu = MU, Sigma = Sigma, empirical = TRUE)[, -1, drop = FALSE]
    BETA     <- cbind(c(sigma_YX %*% solve(Sigma_XX)))
    True_Y   <- DATA_Pop %*% BETA
    for (i in seq_len(G)) {
      if (print_iter) cat(i, "\n")
      Obs_Y <- True_Y + rnorm(N, 0, sqrt(sigma_Y * (1 - true_R2)))
      fit   <- lm(Obs_Y ~ DATA_Pop)
      sm    <- summary(fit)
      Results[i, 1] <- sm$r.squared
      Results[i, 2] <- sm$fstatistic[1]
    }
  }

  Results_df <- as.data.frame(Results)
  if (save) {
    message("Simulation results will be saved to a .csv file; overwriting an existing file of the same name.")
    utils::write.csv(Results_df, filename, row.names = FALSE)
  }

  empirical_power <- mean(Results_df$f_stat >= f_crit, na.rm = TRUE)
  analytic_power  <- ss_power_R2(population_R2 = true_R2, p = p,
                                 alpha_level = alpha_level,
                                 specified_N = N,
                                 random_predictors = random_predictors)$value[2]

  Summary <- data.frame(
    term = c("total_N", "empirical_power", "analytic_power",
             "mean_R2", "median_R2", "sd_R2",
             "mean_F", "median_F", "sd_F", "F_crit",
             "p", "true_R2", "estimated_R2",
             "desired_power", "alpha_level"),
    value = c(N, empirical_power, analytic_power,
              mean(Results_df$r2, na.rm = TRUE),
              median(Results_df$r2, na.rm = TRUE),
              sqrt(var(Results_df$r2, na.rm = TRUE)),
              mean(Results_df$f_stat, na.rm = TRUE),
              median(Results_df$f_stat, na.rm = TRUE),
              sqrt(var(Results_df$f_stat, na.rm = TRUE)),
              f_crit,
              p, true_R2,
              if (is.null(estimated_R2)) NA_real_ else estimated_R2,
              if (is.null(estimated_R2)) NA_real_ else desired_power,
              alpha_level)
  )

  class(Summary) <- c("dmar_ss_power_sensitivity", class(Summary))
  return(.as_dmar_tbl(Summary))
}
