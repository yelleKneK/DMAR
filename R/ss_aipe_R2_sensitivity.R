#' Sensitivity Analysis for Sample Size Planning With the Goal of Accuracy in Parameter Estimation (I.e., a Narrow Observed Confidence Interval)
#'
#' @description
#' Given \code{estimated_R2} and \code{true_R2}, one can perform a sensitivity analysis to determine the
#' effect of a misspecified population squared multiple correlation coefficient using the Accuracy in
#' Parameter Estimation (AIPE) approach to sample size planning. The function evaluates the effect of
#' a misspecified \code{true_R2} on the width of obtained confidence intervals.
#'
#' @param true_R2 Value of the population squared multiple correlation coefficient
#' @param estimated_R2 Value of the estimated (for sample size planning) squared multiple correlation coefficient
#' @param w Full confidence interval width of interest
#' @param p Number of predictors
#' @param random_predictors Whether or not the sample size procedure and the simulation itself should be based on random (set to \code{TRUE}) or fixed predictors (set to \code{FALSE})
#' @param specified_N Selected sample size to use in order to determine distributional properties at a given value of sample size
#' @param assurance Parameter to ensure confidence interval width with a specified degree of certainty
#' @param conf_level Confidence interval coverage (symmetric coverage)
#' @param generate_random_predictors Specify whether the simulation should be based on random (default) or fixed regressors.
#' @param rho_yx Value of the correlation between \emph{y} (dependent variable) and each of the \emph{x} variables (independent variables)
#' @param rho_xx Value of the correlation among the \emph{x} variables (independent variables)
#' @param G Number of generations (i.e., replications) of the simulation
#' @param print_iter Should the iteration number (between 1 and \code{G}) during the run of the function
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#' @param ... for modifying parameters of functions this function calls upon
#'
#' @details
#' When \code{estimated_R2}=\code{true_R2}, the results are that of a simulation study when all assumptions
#' are satisfied. Rather than specifying \code{estimated_R2}, one can specify \code{specified_N} to determine
#' the results of a particular sample size (when doing this \code{estimated_R2} cannot be specified).
#'
#' The sample size estimation procedure technically assumes multivariate normal variables (\code{p}+1) with
#' fixed predictors (\code{x}/independent variables), yet the function assumes random multivariate normal
#' predictors (having a \code{p}+1 multivariate distribution). As Gatsonis and Sampson (1989) note in the
#' context of statistical power analysis (recall this function is used in the context of precision), there
#' is little difference in the outcome.
#'
#' In the behavioral, educational, and social sciences, predictor variables are almost always random, and
#' thus \code{random_predictors} should generally be used. \code{random_predictors=TRUE} specifies how both
#' the sample size planning procedure and the confidence intervals are calculated based on the random
#' predictors/regressors. The internal simulation generates random or fixed predictors/regressors based on
#' whether variables predictor variables are random or fixed. However, when \code{random_predictors=FALSE},
#' only the sample size planning procedure and the confidence intervals are calculated based on the
#' parameter. The parameter \code{generate_random_predictors} (where the default is \code{TRUE} so that
#' random predictors/regressors are generated) allows random or fixed predictor variables to be generated.
#' Because the sample size planning procedure and the internal simulation are both specified, for purposes
#' of sensitivity analysis random/fixed can be crossed to examine the effects of specifying sample size
#' based on one but using it on data based on the other.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries are: \code{mean_lower_limit},
#' \code{median_lower_limit}, \code{sd_lower_limit},
#' \code{mean_upper_limit}, \code{median_upper_limit},
#' \code{sd_upper_limit} (summaries of the realized confidence limits);
#' \code{mean_R2}, \code{median_R2}, \code{sd_R2} (summaries of the
#' observed \eqn{R^2}); \code{mean_ci_width_lower},
#' \code{median_ci_width_lower}, \code{sd_ci_width_lower},
#' \code{mean_ci_width_upper}, \code{median_ci_width_upper},
#' \code{sd_ci_width_upper} (summaries of the one-sided widths, measured
#' from the observed \eqn{R^2} to each limit); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the full
#' interval widths); \code{pct_ci_less_w} (proportion of intervals with
#' width at or below the planning target \code{w});
#' \code{pct_ci_miss_low} and \code{pct_ci_miss_high} (tail-specific
#' empirical non-coverage of \code{true_R2});
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); \code{num_probs_with_cis} (number of replications on
#' which a confidence interval could not be obtained); and the input
#' echoes \code{total_N} (the sample size evaluated), \code{p},
#' \code{true_R2}, \code{estimated_R2} (NA when \code{specified_N} was
#' supplied instead), \code{width}, \code{conf_level}, and
#' \code{assurance} (present only when an assurance was supplied). The
#' proportion rows are on the 0 to 1 scale, not percentages.
#'
#' @references
#' Algina, J. & Olejnik, S. (2000). Determining sample size for accurate estimation of the squared
#' multiple correlation coefficient. \emph{Multivariate Behavioral Research, 35}, 119--137.
#'   \doi{10.1207/s15327906mbr3501_5}
#'
#' Gatsonis, C. & Sampson, A. R. (1989). Multiple Correlation: Exact power and sample size calculations.
#' \emph{Psychological Bulletin, 106}(3), 516--524.
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple correlation coefficient:
#' Accuracy in parameter estimation via narrow confidence intervals, \emph{Multivariate Behavioral
#' Research, 43}(4), 524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
#'   interval estimation, power calculations, sample size estimation, and
#'   hypothesis testing in multiple regression. \emph{Behavior Research
#'   Methods, Instruments, & Computers, 24}(4), 581--582.
#'   \doi{10.3758/BF03203611}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ci_R2}}, \code{\link{conf_limits_nct}}, \code{\link{ss_aipe_R2}}
#'
#' @examples
#' # Change 'G' to some large number (e.g., G=10,000)
#' set.seed(113)
#' ss_aipe_R2_sensitivity(true_R2 = .5, estimated_R2 = .4, w = .10, p = 5,
#'                        conf_level = 0.95, G = 25, print_iter = FALSE)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom utils write.csv


ss_aipe_R2_sensitivity <- function(true_R2 = NULL, estimated_R2 = NULL, w = NULL, p = NULL,
                                   random_predictors = TRUE, specified_N = NULL, assurance = NULL, conf_level = .95,
                                   generate_random_predictors = TRUE, rho_yx = .3, rho_xx = .3, G = 10000, print_iter = TRUE,
                                   save = FALSE, filename = "ss_aipe_r2_sensitivity_result.csv", ...) {
  # requireNamespace("MASS", quietly = TRUE)

  if (true_R2 >= 1 || true_R2 <= 0) stop("The values of \'true_R2\' (i.e., the squared multiple correlation coefficient (R^2)) must be between zero and one.")
  if (w == 0 || w >= 1) stop("The width is not specified correctly.")
  if (w == 0 || w >= 1) stop("The width is not specified correctly.")

  if (is.null(estimated_R2) && is.null(specified_N)) stop("You must specify either \'estimated_R2\' or \'specified_N\'.", call. = FALSE)

  if (!is.null(assurance)) {
    if (assurance < 0 || assurance > 1) stop("You must specify either \'assurance\' to be a value between 0 and 1.", call. = FALSE)
  }

  prev_warn <- getOption("warn")
  on.exit(options(warn = prev_warn), add = TRUE)
  options(warn = -1)

  if (!is.null(estimated_R2)) {
    if (estimated_R2 >= 1 || estimated_R2 <= 0) stop("The values of \'estimated_R2\' (i.e., the squared multiple correlation coefficient (R^2)) must be between zero and one.")
    N <- ss_aipe_R2(population_R2 = estimated_R2, conf_level = conf_level, width = w, which_width = "Full", p = p, assurance = assurance, random_predictors = random_predictors)[1, 2]
  } else {
    N <- specified_N
  }
  #############################################################################################################################

  # Means (arbitrary)
  MU <- rep(0, p + 1)

  # Correlation between Y and the X variables (arbitrary for plausible scenarios)
  sigma_YX <- rbind(rep(rho_yx, p))
  sigma_XY <- t(sigma_YX)

  # Correlation among the predictors (arbitrary for plausible scenarios).
  Sigma_XX <- matrix(rep(rho_xx, p^2), nrow = p, ncol = p)
  diag(Sigma_XX) <- 1

  # Defines the numerator so that the desired P^2 (Rho Squared; Population multiple correlation coefficient) can be obtained.
  Numerator_P_Square <- (sigma_YX %*% solve(Sigma_XX) %*% sigma_XY)

  # Define the variance of Y so that the pop. mult. cor. coef. is as specified.
  sigma_Y <- Numerator_P_Square / true_R2

  Sigma <- rbind(c(sigma_Y, sigma_YX), cbind(sigma_XY, Sigma_XX))
  #############################################################################################################################


  R_Square_Results <- matrix(NA, G, 3)
  colnames(R_Square_Results) <- c("lower_limit", "observed_r2", "upper_limit")

  if (generate_random_predictors == TRUE) {
    for (i in 1:G)
    {
      if (print_iter == TRUE) cat(c(i), "\n")
      DATA <- MASS::mvrnorm(N, mu = MU, Sigma = Sigma)

      Regression_Results <- lm(DATA[, 1] ~ DATA[, -1])
      Summary_Regression_Results <- summary(Regression_Results)

      R_Square_Results[i, 2] <- Summary_Regression_Results$r.squared

      CI_Limits_R2 <- try(ci_R2(R2 = R_Square_Results[i, 2], conf_level = conf_level, N = N, p = p, random_predictors = random_predictors))

      R_Square_Results[i, 1] <- CI_Limits_R2[1, 2]
      R_Square_Results[i, 3] <- CI_Limits_R2[3, 2]
    }
  }

  if (generate_random_predictors == FALSE) {
    DATA_Pop_Cov_Structure <- MASS::mvrnorm(N, mu = MU, Sigma = Sigma, empirical = TRUE)[, -1]
    BETA <- cbind(c(sigma_YX %*% solve(Sigma_XX)))
    True_Y <- DATA_Pop_Cov_Structure %*% BETA

    for (i in 1:G)
    {
      if (print_iter == TRUE) cat(c(i), "\n")

      # So, only Y is random from sample to sample.
      Obs_Y <- True_Y + rnorm(N, 0, sqrt(sigma_Y * (1 - true_R2)))

      Regression_Results <- lm(Obs_Y ~ DATA_Pop_Cov_Structure)

      Summary_Regression_Results <- summary(Regression_Results)

      R_Square_Results[i, 2] <- Summary_Regression_Results$r.squared
      # print(Summary_Regression_Results$r.squared)
      CI_Limits_R2 <- try(ci_R2(R2 = R_Square_Results[i, 2], conf_level = conf_level, N = N, p = p, random_predictors = random_predictors))

      R_Square_Results[i, 1] <- CI_Limits_R2[1, 2]
      R_Square_Results[i, 3] <- CI_Limits_R2[3, 2]
    }
  }


  # Summary Section
  # Tail-specific and overall empirical non-coverage of true_R2 (type I error).
  Miss_Low  <- mean(true_R2 <= R_Square_Results[, 1], na.rm = TRUE)
  Miss_High <- mean(true_R2 >= R_Square_Results[, 3], na.rm = TRUE)
  Total_Type_I_Error <- mean(
    (true_R2 <= R_Square_Results[, 1]) | (true_R2 >= R_Square_Results[, 3]),
    na.rm = TRUE
  )

  Lower_Width_CI <- R_Square_Results[, 2] - R_Square_Results[, 1]
  Upper_Width_CI <- R_Square_Results[, 3] - R_Square_Results[, 2]
  Width_CI <- Lower_Width_CI + Upper_Width_CI
  #############################################################

  Results <- data.frame(
    lower_limit = R_Square_Results[, 1], r2 = R_Square_Results[, 2], upper_limit = R_Square_Results[, 3],
    lower_width_ci = Lower_Width_CI, upper_width_ci = Upper_Width_CI, width_ci = Width_CI
  )

  if (save) {
    result_file <- filename
    message("Simulation results will be saved to a .csv file; overwriting a file of the same name if it exists in the directory.")
    utils::write.csv(Results, result_file, row.names = FALSE)
  }

  Num_Probs_with_CIs <- G - length(na.omit(Results$width_ci))

  Summary <- data.frame(
    term = c(
      "mean_lower_limit", "median_lower_limit", "sd_lower_limit", "mean_upper_limit", "median_upper_limit",
      "sd_upper_limit", "mean_R2", "median_R2", "sd_R2", "mean_ci_width_lower", "median_ci_width_lower",
      "sd_ci_width_lower", "mean_ci_width_upper", "median_ci_width_upper", "sd_ci_width_upper",
      "mean_ci_width", "median_ci_width", "sd_ci_width", "pct_ci_less_w",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error", "num_probs_with_cis",
      "total_N", "p", "true_R2", "estimated_R2", "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(Results$lower_limit, na.rm = TRUE), median(Results$lower_limit, na.rm = TRUE),
      sqrt(var(Results$lower_limit, na.rm = TRUE)), mean(Results$upper_limit, na.rm = TRUE),
      median(Results$upper_limit, na.rm = TRUE), sqrt(var(Results$upper_limit, na.rm = TRUE)),
      mean(Results$r2, na.rm = TRUE), median(Results$r2, na.rm = TRUE), sqrt(var(Results$r2, na.rm = TRUE)),
      mean(Results$lower_width_ci, na.rm = TRUE), median(Results$lower_width_ci, na.rm = TRUE),
      sqrt(var(Results$lower_width_ci, na.rm = TRUE)), mean(Results$upper_width_ci, na.rm = TRUE),
      median(Results$upper_width_ci, na.rm = TRUE), sqrt(var(Results$upper_width_ci, na.rm = TRUE)),
      mean(Results$width_ci, na.rm = TRUE), median(Results$width_ci, na.rm = TRUE),
      sqrt(var(Results$width_ci, na.rm = TRUE)), mean(Width_CI <= w, na.rm = TRUE),
      Miss_Low, Miss_High, Total_Type_I_Error, Num_Probs_with_CIs,
      N, p, true_R2,
      if (is.null(estimated_R2)) NA_real_ else estimated_R2,
      w, conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  return(.as_dmar_tbl(Summary, conf_level = conf_level))
}
