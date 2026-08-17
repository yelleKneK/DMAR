#' Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the (Standardized and Unstandardized) Regression Coefficient
#'
#' @description
#' This function performs a sensitivity analysis when planning sample size from the Accuracy in Parameter
#' Estimation Perspective for the standardized or unstandardized regression coefficient.
#'
#' @param true_var_Y Population variance of the dependent variable (\emph{Y})
#' @param true_cov_YX Population covariances vector between the \code{p} predictor variables and the dependent variable (\emph{Y})
#' @param true_cov_XX Population covariance matrix of the \code{p} predictor variables
#' @param estimated_var_Y Estimated variance of the dependent variable (\emph{Y})
#' @param estimated_cov_YX Estimated covariances vector between the \code{p} predictor variables and the dependent variable (\code{Y})
#' @param estimated_cov_XX Estimated Population covariance matrix of the \code{p} predictor variables
#' @param specified_N Directly specified sample size (instead of planning one from the estimated covariance structure)
#' @param which_predictor Identifies which of the \emph{p} predictors is of interest
#' @param w desired Confidence interval width for the regression coefficient of interest
#' @param noncentral Specify with a \code{TRUE} or \code{FALSE} statement whether or not the noncentral approach to sample size planning should be used
#' @param standardize Specify with a \code{TRUE} or \code{FALSE} statement whether or not the regression coefficient will be standardized
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param assurance Degree of certainty that the obtained confidence interval will be sufficiently narrow
#' @param G The number of generations (i.e., replications) of the simulation within the function
#' @param print_iter Specify with a \code{TRUE}/\code{FALSE} statement if the iteration number should be printed as the simulation within the function runs
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#'
#' @details
#' Direct specification of \code{true_cov_YX} and \code{true_cov_XX} is necessary, even if one is
#' interested in a single regression coefficient, so that the covariance/correlation structure can be
#' specified when the simulation within the function runs.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries are: \code{mean_b_j},
#' \code{median_b_j}, \code{sd_b_j} (summaries of the realized
#' regression-coefficient point estimates); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the realized
#' interval widths); \code{pct_ci_less_w} (proportion of intervals at or
#' below the planning target \code{w}); \code{pct_ci_miss_low} and
#' \code{pct_ci_miss_high} (tail-specific empirical non-coverage of the
#' population coefficient); \code{total_type_I_error} (overall empirical
#' non-coverage, the sum of the two tails); \code{mean_R2},
#' \code{median_R2}, \code{sd_R2} (summaries of the realized squared
#' multiple correlation coefficient); and the input echoes
#' \code{total_N} (the sample size evaluated), \code{p},
#' \code{which_predictor}, \code{true_b_j} and \code{estimated_b_j} (the
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
#' Note that when the true and estimated covariance structures agree
#' (\code{true_cov_YX} equals \code{estimated_cov_YX} and \code{true_cov_XX}
#' equals \code{estimated_cov_XX}), the results are not literally from a
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
#' Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' cov_YX <- c(0.4, 0.3)
#' ss_aipe_reg_coef_sensitivity(
#'   true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
#'   estimated_var_Y = 1, estimated_cov_YX = cov_YX, estimated_cov_XX = Sigma_X,
#'   which_predictor = 1, w = 0.20, conf_level = 0.95,
#'   G = 100, print_iter = FALSE
#' )
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_reg_coef_sensitivity <- function(true_var_Y = NULL, true_cov_YX = NULL, true_cov_XX = NULL,
                                         estimated_var_Y = NULL, estimated_cov_YX = NULL, estimated_cov_XX = NULL, specified_N = NULL,
                                         which_predictor = 1, w = NULL, noncentral = FALSE, standardize = FALSE, conf_level = .95,
                                         assurance = NULL, G = 1000, print_iter = TRUE, save = FALSE,
                                         filename = "ss_aipe_reg_coef_sensitivity_result.csv") {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

  if (noncentral == TRUE && is.null(true_var_Y)) true_var_Y <- 1

  if (is.null(w)) stop("You must specify \'w\' (i.e., a confidence interval width).")
  width <- w
  if (is.null(conf_level)) stop("You must specify a confidence level (i.e., 1 - Type I error rate).")
  if (is.null(G)) stop("You must specify 'G/' (i.e., the number of generations of the simulation).")

  if (is.null(true_cov_XX)) stop("You must specify 'true_cov_XX' (i.e., the covariance matrix of the predictors).")
  if (is.null(true_cov_YX)) stop("You must specify 'true_cov_YX' (i.e., the covariance vector of the predictors with the dependent variable).")

  if ((sum(round(true_cov_XX, 5) == round(t(true_cov_XX), 5))) != (dim(true_cov_XX)[1] * dim(true_cov_XX)[2])) stop("The correlation matrix, \'true_cov_XX\' should be symmetric.")


  p <- dim(true_cov_XX)[1]

  if (!is.null(estimated_cov_XX)) {
    if ((sum(round(estimated_cov_XX, 5) == round(t(estimated_cov_XX), 5))) != (dim(estimated_cov_XX)[1] * dim(estimated_cov_XX)[2])) stop("The covariance matrix, \'estimated_cov_XX\' should be symmetric.")
  }

  if (is.null(estimated_var_Y)) estimated_var_Y <- true_var_Y
  if (is.null(estimated_cov_XX)) estimated_cov_XX <- true_cov_XX
  if (is.null(estimated_cov_YX)) estimated_cov_YX <- true_cov_YX

  Estimated_Sigma <- cbind(c(estimated_var_Y, estimated_cov_YX), rbind(estimated_cov_YX, estimated_cov_XX))
  sigma_Y <- sqrt(Estimated_Sigma[1, 1])
  sigma_X <- sqrt(Estimated_Sigma[(1 + which_predictor), (1 + which_predictor)])
  Estimated_Rho2_Y_X <- (estimated_cov_YX %*% solve(estimated_cov_XX) %*% estimated_cov_YX) / (sigma_Y^2)
  Estimated_Rho2_j_X_without_j <- 1 - ((solve(estimated_cov_XX)[which_predictor, which_predictor] * estimated_cov_XX)[which_predictor, which_predictor])^(-1)
  Estimated_b_j <- (solve(estimated_cov_XX) %*% estimated_cov_YX)[which_predictor]

  # Covariance structure.
  True_Sigma <- cbind(c(true_var_Y, true_cov_YX), rbind(true_cov_YX, true_cov_XX))

  True_Rho2_Y_X <- (True_Sigma[1, -1] %*% solve(True_Sigma[-1, -1]) %*% True_Sigma[-1, 1]) / (True_Sigma[1, 1])
  True_Rho2_j_X_without_j <- 1 - ((solve(true_cov_XX)[which_predictor, which_predictor] * true_cov_XX)[which_predictor, which_predictor])^(-1)
  True_b_j <- (solve(True_Sigma[-1, -1]) %*% True_Sigma[-1, 1])[which_predictor]

  if (True_Rho2_Y_X > 1) stop("You have specified an impossible correlational structure of \'true_cov_XX\' and/or \'true_cov_YX\' (the multiple R square is above 1).")
  if (Estimated_Rho2_Y_X > 1) stop("You have specified an impossible correlational structure of \'estimated_cov_XX\' and/or \'estimated_cov_YX\' (the multiple R square is above 1).")

  # See if this needs to be modified
  if (is.null(specified_N)) {
    Estimated_Sigma_as_Cor <- cov2cor(Estimated_Sigma)
    N <- ss_aipe_reg_coef(
      width = width, rho_XX = Estimated_Sigma_as_Cor[2:(p + 1), 2:(p + 1)], rho_YX = Estimated_Sigma_as_Cor[1, 2:(p + 1)], which_predictor = which_predictor,
      conf_level = conf_level, noncentral = noncentral, assurance = assurance, sigma_Y = Estimated_Sigma[1, 1]^.5, sigma_X = (Estimated_Sigma[(1 + which_predictor), (1 + which_predictor)])^.5
    )[1, 2]
  } else {
    N <- specified_N
  }

  # Means (arbitrary)
  MU <- rep(0, p + 1)


  # Begin simulation.
  Results <- matrix(NA, G, 6)
  colnames(Results) <- c("b_j", "ll_ci_beta_j", "ul_ci_beta_j", "r_2", "se_b_j", "t_for_b_j")
  for (i in 1:G)
  {
    if (print_iter == TRUE) cat(c(i), "\n")
    DATA <- MASS::mvrnorm(N, mu = MU, Sigma = True_Sigma)

    if (standardize == TRUE) DATA <- scale(DATA)

    Regression_Results <- lm(DATA[, 1] ~ DATA[, -1])
    Summary_Results <- summary(Regression_Results)

    b_j <- coef(Summary_Results)[(which_predictor + 1), 1]
    SE_b_j <- coef(Summary_Results)[(which_predictor + 1), 2]

    if (noncentral == FALSE) CI_Lims <- ci_reg_coef(b_j = b_j, SE_b_j = SE_b_j, s_Y = (var(DATA[, 1]))^.5, s_X = (var(DATA[, 1 + which_predictor]))^.5, N = dim(DATA)[1], p = (dim(DATA)[2] - 1), R2_Y_X = NULL, R2_j_X_without_j = NULL, conf_level = conf_level, R2_Y_X_without_j = NULL, t_value = NULL, alpha_lower = NULL, alpha_upper = NULL, noncentral = FALSE)
    if (noncentral == TRUE) CI_Lims <- ci_reg_coef(b_j = b_j, SE_b_j = SE_b_j, s_Y = (var(DATA[, 1]))^.5, s_X = (var(DATA[, 1 + which_predictor]))^.5, N = dim(DATA)[1], p = (dim(DATA)[2] - 1), R2_Y_X = NULL, R2_j_X_without_j = NULL, conf_level = conf_level, R2_Y_X_without_j = NULL, t_value = NULL, alpha_lower = NULL, alpha_upper = NULL, noncentral = TRUE)


    Results[i, 1] <- b_j
    Results[i, 2] <- CI_Lims[which(CI_Lims$term == "lower_limit"), 2]
    Results[i, 3] <- CI_Lims[which(CI_Lims$term == "upper_limit"), 2]
    Results[i, 4] <- Summary_Results$r.squared
    Results[i, 5] <- SE_b_j
    Results[i, 6] <- b_j / SE_b_j
  }

  # End Simulation.

  Results <- as.data.frame(Results)

  if (save) {
    result_file <- filename
    # print("Simulation results will be saved to a .csv file")
    suppressWarnings(file_exist <- try(utils::read.csv(result_file), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      utils::write.table(Results, result_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      cat("A file in the local directory has the same name as the file where simulation", "\n", "results will be saved to. Simulation results will be appended to this file.", "\n", sep = "")
    } else {
      utils::write.table(Results, result_file, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  Summary <- data.frame(
    term = c(
      "mean_b_j", "median_b_j", "sd_b_j", "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "mean_R2", "median_R2", "sd_R2",
      "total_N", "p", "which_predictor",
      "true_b_j", "estimated_b_j",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(Results[, 1]), median(Results[, 1]), (var(Results[, 1]))^.5,
      mean(Results[, 3] - Results[, 2]), median(Results[, 3] - Results[, 2]),
      (var(Results[, 3] - Results[, 2]))^.5, mean((Results[, 3] - Results[, 2]) <= w),
      mean(True_b_j < Results[, 2]), mean(True_b_j > Results[, 3]),
      (mean((True_b_j < Results[, 2]) | (True_b_j > Results[, 3]))),
      mean(Results[, 4]), median(Results[, 4]), (var(Results[, 4]))^.5,
      N, p, which_predictor,
      True_b_j, Estimated_b_j,
      w, conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  return(.as_dmar_tbl(Summary, conf_level = conf_level))
}
