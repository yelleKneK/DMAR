#' @title Sensitivity Analysis for Sample Size Planning for the (Unstandardized) Contrast in Randomized ANCOVA
#' From the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE)
#' Perspective for the (unstandardized) contrast in randomized ANCOVA design.
#'
#' @param true_var_error_ancova population error variance of the ANCOVA model
#' @param supposed_var_error_ancova the supposed error variance of the ANCOVA model, the planning value given to \code{\link{ss_aipe_c_ancova}} to choose the sample size
#' @param true_var_error_anova population error variance of the ANOVA model (i.e., excluding the covariate)
#' @param supposed_var_error_anova the supposed error variance of the ANOVA model (i.e., excluding the covariate), a planning value used with \code{supposed_cor}
#' @param true_cor population correlation coefficient of the response and the covariate
#' @param supposed_cor the supposed correlation coefficient of the response and the covariate, a planning value used with \code{supposed_var_error_anova}
#' @param G number of generations (i.e., replications) of the simulation
#' @param true_mean_Y vector that contains the response's population mean of each group
#' @param true_sd_Y the population standard deviation of the response
#' @param true_mean_X the population mean of the covariate
#' @param true_sd_X the population standard deviation of the covariate
#' @param c_weights the contrast weights
#' @param width the desired full width of the obtained confidence interval
#' @param conf_level the desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param filename an optional path for a comma separated file recording
#'   every replication (the realized contrast, its full and
#'   covariate-ignoring standard errors and their ratio, the interval
#'   width, and the three non-coverage indicators): nothing is written
#'   when \code{filename} is \code{NULL} (the default), a new file with a
#'   header row is created otherwise, an existing file at that path is
#'   appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' The arguments \code{true_mean_Y}, \code{true_mean_X}, \code{true_sd_Y}, and \code{true_sd_X} are used to generate random data in the simulations
#' for the sensitivity analysis. The value of \code{true_sd_Y} should be the same as the square root of \code{true_var_error_anova}.
#'
#' So far this function is based on one-covariate randomized ANCOVA design only. The argument \code{true_mean_X} should be a single number,
#' because it is assumed that the population mean of the covariate is equal across groups in randomized ANCOVA.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries are: \code{mean_contrast},
#' \code{median_contrast}, \code{sd_contrast} (summaries of the realized
#' unstandardized contrast); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the realized
#' interval widths); \code{pct_ci_less_width} (proportion of intervals
#' narrower than the planning target \code{width});
#' \code{pct_ci_miss_low} and \code{pct_ci_miss_high} (tail-specific
#' empirical non-coverage of the population contrast);
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); \code{mean_se_ratio} (mean ratio of the contrast
#' standard error that ignores the covariate-imbalance term to the full
#' ANCOVA standard error); and the input echoes \code{n_per_group},
#' \code{total_N}, \code{true_contrast} (the population contrast implied by
#' \code{true_mean_Y} and \code{c_weights}), \code{supposed_var_error_ancova} (as
#' supplied or as resolved from \code{supposed_var_error_anova} and
#' \code{supposed_cor}), \code{true_cor}, \code{width}, \code{conf_level}, and
#' \code{assurance} (present only when an assurance was supplied). The
#' proportion rows are on the 0 to 1 scale, not percentages. The
#' per-replication vectors (\code{contrast}, \code{se_contrast},
#' \code{se_contrast_restricted}, \code{width_obs}) are not returned; they
#' are written to the comma separated file named by \code{filename} when
#' one is supplied.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Monte Carlo sensitivity sweep; G is small here so the example runs quickly.
#' # Raise G (e.g., G = 1000 or more) for a stable sensitivity analysis.
#' set.seed(113)
#' ss_aipe_c_ancova_sensitivity(true_var_error_ancova = 30,
#'                              supposed_var_error_ancova = 30,
#'                              true_cor = .2, true_mean_Y = c(10, 12, 15, 13),
#'                              true_mean_X = 2, true_sd_X = 1.3, true_sd_Y = 2,
#'                              c_weights = c(1, 0, -1, 0), width = 3, G = 25)
#'
#' ss_aipe_c_ancova_sensitivity(true_var_error_anova = 36,
#'                              supposed_var_error_anova = 36,
#'                              true_cor = .2, supposed_cor = .2,
#'                              true_mean_Y = c(10, 12, 15, 13), true_mean_X = 2,
#'                              true_sd_X = 1.3, true_sd_Y = 6,
#'                              c_weights = c(1, 0, -1, 0), width = 3,
#'                              assurance = NULL, G = 25)
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
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export


ss_aipe_c_ancova_sensitivity <- function(true_var_error_ancova = NULL, supposed_var_error_ancova = NULL, true_var_error_anova = NULL, supposed_var_error_anova = NULL,
                                         true_cor, supposed_cor = NULL, G = 10000, true_mean_Y, true_sd_Y, true_mean_X, true_sd_X, c_weights, width, conf_level = .95, assurance = NULL,
                                         filename = NULL,
                                         seed = NULL) {
  if (!requireNamespace("MASS", quietly = TRUE)) .dmar_refuse(
    "The package 'MASS' is needed; please install the package and try again.",
    "package_required", package = "MASS")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  if (!is.null(true_sd_Y) && !is.null(true_var_error_anova)) {
    if (true_sd_Y != sqrt(true_var_error_anova)) .dmar_refuse_argument(
      "'true_sd_Y' and 'true_var_error_anova' should be the same",
      arg = c("true_sd_Y", "true_var_error_anova"))
  }

  if (is.null(supposed_var_error_ancova)) {
    if (is.null(supposed_var_error_anova) || is.null(supposed_cor)) {
      .dmar_refuse_argument(
        "Please specify the either 'supposed_var_error_ancova', or both 'supposed_var_error_anova' and 'supposed_cor'",
        arg = c("supposed_var_error_ancova", "supposed_var_error_anova", "supposed_cor"))
    } else {
      supposed_var_error_ancova <- supposed_var_error_anova * (1 - supposed_cor^2)
    }
  }
  # All-zero weights define no contrast; they gave a summary of zeros.
  .dmar_check_contrast_nonzero(c_weights)
  if (length(true_mean_Y) != length(c_weights)) .dmar_refuse_argument(
    "'true_mean_Y' must be a vector that contains the mean of each group",
    arg = "true_mean_Y")

  n <- ss_aipe_c_ancova(
    true_var_error_ancova = supposed_var_error_ancova, c_weights = c_weights, width = width, conf_level = conf_level,
    assurance = assurance
  )[1, 2]
  J <- length(c_weights)
  Psi_pop <- sum(c_weights * true_mean_Y)

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
    random_data <- simulate_ancova_data(true_mean_Y = true_mean_Y, true_mean_X = true_mean_X,
                                        true_sd_Y = true_sd_Y, true_sd_X = true_sd_X,
                                        true_cor = true_cor, a = J, n_per_group = n)
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
    contrast = Psi_obs, se_contrast = se_Psi, se_contrast_restricted = se_Psi_res, ratio = se_res_vs_se_full,
    width_obs = width_obs, type_I_error = Type_I_Error, type_I_error_upper = Type_I_Error_Upper,
    type_I_error_lower = Type_I_Error_Lower
  )

  if (!is.null(filename)) {
    if (file.exists(filename) && file.size(filename) > 0) {
      message("The file '", filename, "' already exists; the simulation results are appended to it.")
      utils::write.table(Results, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    } else {
      utils::write.table(Results, filename, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  # Summary
  out <- data.frame(
    term = c(
      "mean_contrast", "median_contrast", "sd_contrast",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "pct_ci_less_width",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "mean_se_ratio",
      "n_per_group", "total_N",
      "true_contrast", "supposed_var_error_ancova", "true_cor",
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
      Psi_pop, supposed_var_error_ancova, true_cor,
      width, conf_level,
      if (!is.null(assurance)) assurance
    )
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(width_obs), conf_level)
}
