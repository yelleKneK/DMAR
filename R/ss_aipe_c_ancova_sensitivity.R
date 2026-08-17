#' @title Sensitivity Analysis for Sample Size Planning for the (Unstandardized) Contrast in Randomized ANCOVA
#' From the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE)
#' Perspective for the (unstandardized) contrast in randomized ANCOVA design.
#'
#' @param true_error_var_ancova population error variance of the ANCOVA model
#' @param est_error_var_ancova estimated error variance of the ANCOVA model
#' @param true_error_var_anova population error variance of the ANOVA model (i.e., excluding the covariate)
#' @param est_error_var_anova estimated error variance of the ANOVA model (i.e., excluding the covariate)
#' @param rho population correlation coefficient of the response and the covariate
#' @param est_rho estimated correlation coefficient of the response and the covariate
#' @param G number of generations (i.e., replications) of the simulation
#' @param mu_y vector that contains the response's population mean of each group
#' @param sigma_y the population standard deviation of the response
#' @param mu_x the population mean of the covariate
#' @param sigma_x the population standard deviation of the covariate
#' @param c_weights the contrast weights
#' @param width the desired full width of the obtained confidence interval
#' @param conf_level the desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#'
#' @details
#' The arguments \code{mu_y}, \code{mu_x}, \code{sigma_y}, and \code{sigma_x} are used to generate random data in the simulations
#' for the sensitivity analysis. The value of \code{sigma_y} should be the same as the square root of \code{true_error_var_anova}.
#'
#' So far this function is based on one-covariate randomized ANCOVA design only. The argument \code{mu_x} should be a single number,
#' because it is assumed that the population mean of the covariate is equal across groups in randomized ANCOVA.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries are: \code{mean_psi},
#' \code{median_psi}, \code{sd_psi} (summaries of the realized
#' unstandardized contrast); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the realized
#' interval widths); \code{pct_ci_less_w} (proportion of intervals
#' narrower than the planning target \code{width});
#' \code{pct_ci_miss_low} and \code{pct_ci_miss_high} (tail-specific
#' empirical non-coverage of the population contrast);
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); \code{mean_se_ratio} (mean ratio of the contrast
#' standard error that ignores the covariate-imbalance term to the full
#' ANCOVA standard error); and the input echoes \code{n_per_group},
#' \code{total_N}, \code{true_psi} (the population contrast implied by
#' \code{mu_y} and \code{c_weights}), \code{est_error_var_ancova} (as
#' supplied or as resolved from \code{est_error_var_anova} and
#' \code{est_rho}), \code{rho}, \code{width}, \code{conf_level}, and
#' \code{assurance} (present only when an assurance was supplied). The
#' proportion rows are on the 0 to 1 scale, not percentages. The
#' per-replication vectors (\code{psi_obs}, \code{se_psi},
#' \code{se_psi_restricted}, \code{width_obs}) are not returned; they
#' are written to the CSV named by \code{filename} when
#' \code{save = TRUE}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Monte Carlo sensitivity sweep; G is small here so the example runs quickly.
#' # Raise G (e.g., G = 1000 or more) for a stable sensitivity analysis.
#' set.seed(113)
#' ss_aipe_c_ancova_sensitivity(true_error_var_ancova=30,
#'                              est_error_var_ancova=30, rho=.2, mu_y=c(10,12,15,13), mu_x=2,
#'                              G=50, sigma_x=1.3, sigma_y=2, c_weights=c(1,0,-1,0), width=3)
#'
#' ss_aipe_c_ancova_sensitivity(true_error_var_anova=36,
#'                              est_error_var_anova=36, rho=.2, est_rho=.2, G=50,
#'                              mu_y=c(10,12,15,13), mu_x=2, sigma_x=1.3, sigma_y=6,
#'                              c_weights=c(1,0,-1,0), width=3, assurance=NULL)
#'
#' @references
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_c_ancova_sensitivity <- function(true_error_var_ancova = NULL, est_error_var_ancova = NULL, true_error_var_anova = NULL, est_error_var_anova = NULL,
                                         rho, est_rho = NULL, G = 10000, mu_y, sigma_y, mu_x, sigma_x, c_weights, width, conf_level = .95, assurance = NULL,
                                         save = FALSE, filename = "ss_aipe_c_ancova_sensitivity_result.csv") {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

  if (!is.null(sigma_y) && !is.null(true_error_var_anova)) {
    if (sigma_y != sqrt(true_error_var_anova)) stop("'sigma_y' and 'true_error_var_anova' should be the same")
  }

  if (is.null(est_error_var_ancova)) {
    if (is.null(est_error_var_anova) || is.null(est_rho)) {
      stop("Please specify the either 'est_error_var_ancova', or both 'est_error_var_anova' and 'est_rho'")
    } else {
      est_error_var_ancova <- est_error_var_anova * (1 - est_rho^2)
    }
  }
  if (length(mu_y) != length(c_weights)) stop("'mu_y' must be a vector that contains the mean of each group")

  n <- ss_aipe_c_ancova(
    error_var_ancova = est_error_var_ancova, c_weights = c_weights, width = width, conf_level = conf_level,
    assurance = assurance
  )[1, 2]
  J <- length(c_weights)
  Psi_pop <- sum(c_weights * mu_y)

  Type_I_Error <- rep(NA, G)
  Type_I_Error_Upper <- rep(NA, G)
  Type_I_Error_Lower <- rep(NA, G)
  se_Psi <- rep(NA, G)
  se_Psi_res <- rep(NA, G)
  se_res_vs_se_full <- rep(NA, G)
  Psi_obs <- rep(NA, G)
  width_obs <- rep(NA, G)
  w_vs_omega <- rep(NA, G)

  for (g in 1:G)
  {
    random_data <- simulate_ancova_data(mu_y = mu_y, mu_x = mu_x,
                                        sigma_y = sigma_y, sigma_x = sigma_x,
                                        rho = rho, a = J, n = n)
    y_vector <- random_data$y
    x_vector <- random_data$x
    group    <- random_data$group

    # fit the ANCOVA model with the random data
    fitted_model <- lm(y_vector ~ x_vector + group)

    # calculate adjusted Y means
    beta <- fitted_model$coefficients[2]
    y_bar_adj <- rep(NA, J)
    x_bar     <- rep(NA, J)
    grand_x   <- mean(x_vector)
    for (j in 1:J)
    {
      in_j         <- group == levels(group)[j]
      x_bar[j]     <- mean(x_vector[in_j])
      y_bar_adj[j] <- mean(y_vector[in_j]) - beta * (x_bar[j] - grand_x)
    }

    # calculate the ANCOVA contrast
    Psi_obs[g] <- sum(c_weights * y_bar_adj)

    # extract the MSwithin from the ANCOVA table
    error_var_ancova <- anova(fitted_model)[3, 3]

    # extract SSwithin_x from ANOVA on the covariate
    x_anova <- lm(x_vector ~ group)
    SSwithin_x <- anova(x_anova)[2, 2]

    f_x_numerater <- (sum(c_weights * x_bar))^2
    f_x_denominator <- SSwithin_x
    sample_size_weighted <- sum(c_weights^2) / n

    se_Psi2 <- error_var_ancova * (sample_size_weighted + f_x_numerater / f_x_denominator)
    se_Psi[g] <- (sqrt(se_Psi2))
    se_Psi_res[g] <- sqrt(error_var_ancova * sample_size_weighted)
    se_res_vs_se_full[g] <- se_Psi_res[g] / se_Psi[g]

    # calculate the confidence interval
    alpha <- 1 - conf_level
    nu <- n * J - J - 1
    t_value <- qt(1 - alpha / 2, df = nu)

    ci_Psi <- list(upper = Psi_obs[g] + t_value * se_Psi[g], lower = Psi_obs[g] - t_value * se_Psi[g])
    width_obs[g] <- ci_Psi$upper - ci_Psi$lower

    Type_I_Error_Upper[g] <- Psi_pop > ci_Psi$upper
    Type_I_Error_Lower[g] <- Psi_pop < ci_Psi$lower
    Type_I_Error[g] <- Type_I_Error_Upper[g] | Type_I_Error_Lower[g]

    if (width_obs[g] < width) {
      w_vs_omega[g] <- TRUE
    } else {
      w_vs_omega[g] <- FALSE
    }
  }

  Results <- data.frame(
    psi_obs = Psi_obs, se_psi = se_Psi, se_psi_restricted = se_Psi_res, ratio = se_res_vs_se_full,
    width_obs = width_obs, type_I_error = Type_I_Error, type_I_error_upper = Type_I_Error_Upper,
    type_I_error_lower = Type_I_Error_Lower
  )

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

  # Summary
  out <- data.frame(
    term = c(
      "mean_psi", "median_psi", "sd_psi",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_w",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "mean_se_ratio",
      "n_per_group", "total_N",
      "true_psi", "est_error_var_ancova", "rho",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(Psi_obs), median(Psi_obs), sd(Psi_obs),
      mean(width_obs), median(width_obs), sd(width_obs),
      mean(w_vs_omega),
      mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
      mean(se_res_vs_se_full),
      n, n * J,
      Psi_pop, est_error_var_ancova, rho,
      width, conf_level,
      if (!is.null(assurance)) assurance
    )
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
