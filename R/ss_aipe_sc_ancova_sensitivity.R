#' Sensitivity Analysis for the Sample Size Planning Method for Standardized ANCOVA Contrast
#'
#' @description
#' Sensitivity analysis for the sample size planning method with the goal to obtain sufficiently
#' narrow confidence intervals for standardized ANCOVA complex contrasts.
#'
#' @param true_psi the population standardized ANCOVA contrast
#' @param estimated_psi the estimated standardized ANCOVA contrast
#' @param c_weights the contrast weights
#' @param desired_width the desired full width of the obtained confidence interval
#' @param n_per_group selected sample size to use in order to determine distributional properties of a given value of sample size
#' @param mu_x the population mean for the covariate
#' @param sigma_x the population standard deviation of the covariate
#' @param rho the population correlation coefficient between the response and the covariate
#' @param divisor which error standard deviation to be used in standardizing the contrast; the value can be either \code{"s_ancova"} or \code{"s_anova"}
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity)
#' @param conf_level the desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param G number of generations (i.e., replications) of the simulation
#' @param print_iter to print the current value of the iterations
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#' @param \dots allows one to potentially include parameter values for inner functions
#'
#' @details
#' The sample size planning method this function is based on is developed in the context of simple (i.e., one-response-one-covariate)
#' ANCOVA model and randomized design (i.e., same population covariate mean across groups).
#'
#' An ANCOVA contrast can be standardized in at least two ways: (a) divided by the error standard deviation of the
#' ANOVA model, (b) divided by the error standard deviation of the ANCOVA model. This function can be used to analyze
#' both types of standardized ANCOVA contrasts.
#'
#' The population mean and standard deviation of the covariate does not affect the sample size planning procedure;
#' they can be specified as any values that are considered as reasonable by the user.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across \code{G}
#' replications. The \code{term} entries are: \code{mean_psi},
#' \code{median_psi}, \code{sd_psi} (summaries of the realized
#' standardized ANCOVA contrast); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the full
#' interval widths); \code{mean_ci_width_lower} and
#' \code{mean_ci_width_upper} (mean one-sided widths, measured from the
#' observed contrast to each limit); \code{pct_ci_less_w} (proportion of
#' intervals at or below the target width); \code{pct_ci_miss_low} and
#' \code{pct_ci_miss_high} (tail-specific empirical non-coverage of
#' \code{true_psi}); \code{total_type_I_error} (overall empirical
#' non-coverage, the sum of the two tails); and the input echoes
#' \code{n_per_group}, \code{total_N}, \code{true_psi},
#' \code{estimated_psi} (NA when \code{n_per_group} was supplied
#' instead), \code{rho}, \code{width}, \code{conf_level}, and
#' \code{assurance} (present only when an assurance was supplied). The
#' proportion and Type I error rows are proportions on the 0 to 1
#' scale, not percentages.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals.
#'   \emph{Psychological Methods, 11}(4), 363--385.
#'   \doi{10.1037/1082-989X.11.4.363}
#'
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
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_sc_ancova}}, \code{\link{ss_aipe_sc_sensitivity}}
#'
#' @examples
#' # Sensitivity analysis for a standardized ANCOVA contrast across
#' # three groups, contrast (-1, 0, 1), a covariate-outcome correlation
#' # of 0.4, and a planning target width of 0.5. Sizes are kept small
#' # here so the Monte Carlo sweep runs quickly; raise G for a stable
#' # estimate in practice.
#' set.seed(113)
#' ss_aipe_sc_ancova_sensitivity(
#'   true_psi = 0.5, estimated_psi = 0.5,
#'   c_weights = c(-1, 0, 1),
#'   desired_width = 0.5, rho = 0.4,
#'   conf_level = 0.95, G = 50, print_iter = FALSE
#' )
#'
#' @keywords design
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sc_ancova_sensitivity <- function(true_psi = NULL, estimated_psi = NULL, c_weights,
                                          desired_width = NULL, n_per_group = NULL, mu_x = 0, sigma_x = 1, rho, divisor = "s_ancova",
                                          assurance = NULL, conf_level = .95, G = 10000, print_iter = TRUE,
                                          save = FALSE, filename = "ss_aipe_sc_ancova_sensitivity_result.csv", ...) {
  prev_warn <- getOption("warn")
  on.exit(options(warn = prev_warn), add = TRUE)
  options(warn = -1)
  if (divisor != "s_ancova" && divisor != "s_anova") stop("The argument 'divisor' must be either 's_ancova' or 's_anova'")

  if (divisor == "s_ancova") {
    ##  standardized ANCOVA contrast using s_ancova as divisor

    if (is.null(estimated_psi) && is.null(n_per_group)) stop("You must specify either \'estimated_psi\' or \'n_per_group\' (i.e., the sample size per group).", call. = FALSE)
    if (!is.null(estimated_psi) && !is.null(n_per_group)) stop("You must specify either \'estimated_psi\' or \'n_per_group\' (i.e., the sample size per group), but not both.", call. = FALSE)

    if (abs(sum(c_weights)) > 1e-8) stop("The sum of the contrast weights must be zero")
    if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")

    J <- length(c_weights)
    width <- desired_width
    sigma_ancova <- 2
    sigma_anova <- sqrt(sigma_ancova^2 / (1 - rho^2))
    mu_y <- rep(NA, J)

    if (!is.null(estimated_psi)) {
      n <- ss_aipe_sc_ancova(
        psi_standardized = estimated_psi, c_weights = c_weights, width = width,
        conf_level = conf_level, assurance = assurance, ...
      )[1, 2]
    }
    else {
      n <- n_per_group
    }

    group <- gl(J, n)
    y_bar <- rep(NA, J)
    y_bar_adj <- rep(NA, J)
    cov_means <- rep(NA, J)

    G <- G
    Full_Width <- rep(NA, G)
    Width_from_psi_obs_Lower <- rep(NA, G)
    Width_from_psi_obs_Upper <- rep(NA, G)
    Type_I_Error_Upper <- rep(NA, G)
    Type_I_Error_Lower <- rep(NA, G)
    Type_I_Error <- rep(NA, G)
    Lower_Limit <- rep(NA, G)
    Upper_Limit <- rep(NA, G)
    psi_obs <- rep(NA, G)
    psi_pop <- true_psi

    # create mu_y, the vector that contains the population mean of each group
    # J-1 out of the J means are free to be any value, and are set to 1 here in this function
    # The mean for the 1st group whose contrast weight is not zero is determined by the rest J-1 means
    Psi <- psi_pop * sigma_ancova
    signal <- 1
    for (j in 1:J)
    {
      if (c_weights[j] != 0 && signal == 1) {
        mu_y[j] <- 0
        signal <- 0
        k <- j
        j <- j + 1
      }
      else {
        mu_y[j] <- 1
        j <- j + 1
      }
    }
    mu_y_k <- (Psi - sum(c_weights * mu_y)) / c_weights[k]
    mu_y[k] <- mu_y_k

    for (g in 1:G)
    {
      if (print_iter) cat(c(g), "\n")

      # generate random samples (long-format data.frame: group, y, x)
      random_data <- simulate_ancova_data(mu_y = mu_y, mu_x = mu_x,
                                          sigma_y = sigma_anova, sigma_x = sigma_x,
                                          rho = rho, a = J, n = n)
      y_vector <- random_data$y
      x_vector <- random_data$x
      group    <- random_data$group

      # fit the ANCOVA model with the random data
      fitted_model <- lm(y_vector ~ x_vector + group)

      # extract the error variance from ANCOVA table
      error_var_ancova <- anova(fitted_model)[3, 3]
      s_ancova <- sqrt(error_var_ancova)

      # obtain the regression coefficient and the adjusted means
      beta <- fitted_model$coefficients[2]

      for (j in 1:J)
      {
        in_j <- group == levels(group)[j]
        y_bar[j]     <- mean(y_vector[in_j])
        cov_means[j] <- mean(x_vector[in_j])
        y_bar_adj[j] <- y_bar[j] - beta * (cov_means[j] - mean(x_vector))
      }

      # obtain SSwithin_x from ANOVA on the covariate
      x_anova <- lm(x_vector ~ group)
      SSwithin_x <- anova(x_anova)[2, 2]

      # calculate the observed psi and CI for psi
      psi_obs[g] <- sum(c_weights * y_bar_adj) / s_ancova

      ci_psi <- ci_sc_ancova(
        adj_means = y_bar_adj, s_ancova = s_ancova, c_weights = c_weights, n = n, cov_means = cov_means,
        SSwithin_x = SSwithin_x, conf_level = conf_level
      )

      psi_limit_lower <- ci_psi$value[ci_psi$term == "lower_limit"]
      psi_limit_upper <- ci_psi$value[ci_psi$term == "upper_limit"]

      # evaluate the obtained confidence interval
      Full_Width[g] <- abs(psi_limit_upper - psi_limit_lower)

      Width_from_psi_obs_Lower[g] <- psi_obs[g] - psi_limit_lower
      Width_from_psi_obs_Upper[g] <- psi_limit_upper - psi_obs[g]

      Type_I_Error_Upper[g] <- psi_pop > psi_limit_upper
      Type_I_Error_Lower[g] <- psi_pop < psi_limit_lower
      Type_I_Error[g] <- Type_I_Error_Upper[g] | Type_I_Error_Lower[g]

      Lower_Limit[g] <- psi_limit_lower
      Upper_Limit[g] <- psi_limit_upper
    }

    Results <- data.frame(
      psi_obs = psi_obs, full_width = Full_Width, width_lower = Width_from_psi_obs_Lower,
      width_upper = Width_from_psi_obs_Upper, type_I_error_upper = Type_I_Error_Upper,
      type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, lower_limit = Lower_Limit,
      upper_limit = Upper_Limit
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

    Summary <- data.frame(
      term = c(
        "mean_psi", "median_psi", "sd_psi",
        "mean_ci_width", "median_ci_width", "sd_ci_width",
        "mean_ci_width_lower", "mean_ci_width_upper",
        "pct_ci_less_w",
        "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
        "n_per_group", "total_N",
        "true_psi", "estimated_psi", "rho",
        "width", "conf_level",
        if (!is.null(assurance)) "assurance"
      ),
      value = c(
        mean(psi_obs), median(psi_obs), sd(psi_obs),
        mean(Full_Width), median(Full_Width), sd(Full_Width),
        mean(Width_from_psi_obs_Lower), mean(Width_from_psi_obs_Upper),
        mean(Full_Width <= desired_width),
        mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
        n, n * J,
        true_psi,
        if (is.null(estimated_psi)) NA_real_ else estimated_psi,
        rho,
        if (is.null(desired_width)) NA_real_ else desired_width,
        conf_level,
        if (!is.null(assurance)) assurance
      )
    )

    return(.as_dmar_tbl(Summary, conf_level = conf_level))
  }
  #################################################################################################
  #################################################################################################

  if (divisor == "s_anova")
  ## use s_anova to standardize the ANCOVA contrast

    {
      if (is.null(estimated_psi) && is.null(n_per_group)) stop("You must specify either \'estimated_psi\' or \'n_per_group\' (i.e., the sample size per group ).", call. = FALSE)
      if (!is.null(estimated_psi) && !is.null(n_per_group)) stop("You must specify either \'estimated_psi\' or \'n_per_group\' (i.e., the sample size per group), but not both.", call. = FALSE)

      if (abs(sum(c_weights)) > 1e-8) stop("The sum of the contrast weights must be zero")
      if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")

      width <- desired_width
      sigma_anova <- 2
      sigma_ancova <- sqrt(sigma_anova^2 * (1 - rho^2))

      if (!is.null(estimated_psi)) {
        n <- ss_aipe_sc_ancova(
          psi = estimated_psi * sigma_anova, sigma_ancova = sigma_ancova, sigma_anova = sigma_anova, c_weights = c_weights,
          divisor = "s_anova", width = width, conf_level = conf_level, assurance = assurance, ...
        )[1, 2]
      }
      else {
        n <- n_per_group
      }

      J <- length(c_weights)
      group <- gl(J, n)
      y_bar <- rep(NA, J)
      y_bar_adj <- rep(NA, J)
      cov_means <- rep(NA, J)
      mu_y <- rep(NA, J)

      G <- G
      Full_Width <- rep(NA, G)
      Width_from_psi_obs_Lower <- rep(NA, G)
      Width_from_psi_obs_Upper <- rep(NA, G)
      Type_I_Error_Upper <- rep(NA, G)
      Type_I_Error_Lower <- rep(NA, G)
      Type_I_Error <- rep(NA, G)
      Lower_Limit <- rep(NA, G)
      Upper_Limit <- rep(NA, G)
      psi_obs <- rep(NA, G)
      psi_pop <- true_psi

      # create mu_y, the vector that contains the population mean of each group
      # J-1 out of the J means are free to be any value, and are set to 1 here in this function
      # The mean for the 1st group whose contrast weight is not zero is determined by the rest J-1 means
      Psi <- sigma_anova * psi_pop
      signal <- 1
      for (j in 1:J)
      {
        if (c_weights[j] != 0 && signal == 1) {
          mu_y[j] <- 0
          signal <- 0
          k <- j
          j <- j + 1
        }
        else {
          mu_y[j] <- 1
          j <- j + 1
        }
      }
      mu_y_k <- (Psi - sum(c_weights * mu_y)) / c_weights[k]
      mu_y[k] <- mu_y_k

      for (g in 1:G)
      {
        if (print_iter) cat(c(g), "\n")

        # generate random samples (long-format data.frame: group, y, x)
        random_data <- simulate_ancova_data(mu_y = mu_y, mu_x = mu_x,
                                            sigma_y = sigma_anova, sigma_x = sigma_x,
                                            rho = rho, a = J, n = n)
        y_vector <- random_data$y
        x_vector <- random_data$x
        group    <- random_data$group

        # fit the ANCOVA model with the random data
        fitted_model <- lm(y_vector ~ x_vector + group)

        # extract the error variance from ANCOVA table
        error_var_ancova <- anova(fitted_model)[3, 3]
        s_ancova <- sqrt(error_var_ancova)

        # obtain the regression coefficient and the adjusted means
        beta <- fitted_model$coefficients[2]

        for (j in 1:J)
        {
          in_j <- group == levels(group)[j]
          y_bar[j]     <- mean(y_vector[in_j])
          cov_means[j] <- mean(x_vector[in_j])
          y_bar_adj[j] <- y_bar[j] - beta * (cov_means[j] - mean(x_vector))
        }

        # extract the error variance from the ANOVA on y
        anova_y <- lm(y_vector ~ group)
        error_var_anova <- anova(anova_y)[2, 3]
        s_anova <- sqrt(error_var_anova)

        # obtain SSwithin_x from ANOVA on the covariate
        x_anova <- lm(x_vector ~ group)
        SSwithin_x <- anova(x_anova)[2, 2]

        # calculate the observed psi and CI for psi based on s_anova
        psi_obs[g] <- sum(c_weights * y_bar_adj) / s_anova

        ci_psi <- ci_sc_ancova(
          psi = sum(c_weights * y_bar_adj), s_anova = s_anova, s_ancova = s_ancova, standardizer = "s_anova",
          c_weights = c_weights, n = n, cov_means = cov_means, SSwithin_x = SSwithin_x, conf_level = conf_level
        )

        psi_limit_lower <- ci_psi$value[ci_psi$term == "lower_limit"]
        psi_limit_upper <- ci_psi$value[ci_psi$term == "upper_limit"]

        # evaluate the obtained confidence interval
        Full_Width[g] <- abs(psi_limit_upper - psi_limit_lower)

        Width_from_psi_obs_Lower[g] <- psi_obs[g] - psi_limit_lower
        Width_from_psi_obs_Upper[g] <- psi_limit_upper - psi_obs[g]

        Type_I_Error_Upper[g] <- psi_pop > psi_limit_upper
        Type_I_Error_Lower[g] <- psi_pop < psi_limit_lower
        Type_I_Error[g] <- Type_I_Error_Upper[g] | Type_I_Error_Lower[g]

        Lower_Limit[g] <- psi_limit_lower
        Upper_Limit[g] <- psi_limit_upper
      }

      Results <- data.frame(
        psi_obs = psi_obs, full_width = Full_Width, width_lower = Width_from_psi_obs_Lower,
        width_upper = Width_from_psi_obs_Upper, type_I_error_upper = Type_I_Error_Upper,
        type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, lower_limit = Lower_Limit,
        upper_limit = Upper_Limit
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

      Summary <- data.frame(
        term = c(
          "mean_psi", "median_psi", "sd_psi",
          "mean_ci_width", "median_ci_width", "sd_ci_width",
          "mean_ci_width_lower", "mean_ci_width_upper",
          "pct_ci_less_w",
          "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
          "n_per_group", "total_N",
          "true_psi", "estimated_psi", "rho",
          "width", "conf_level",
          if (!is.null(assurance)) "assurance"
        ),
        value = c(
          mean(psi_obs), median(psi_obs), sd(psi_obs),
          mean(Full_Width), median(Full_Width), sd(Full_Width),
          mean(Width_from_psi_obs_Lower), mean(Width_from_psi_obs_Upper),
          mean(Full_Width <= desired_width),
          mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
          n, n * J,
          true_psi,
          if (is.null(estimated_psi)) NA_real_ else estimated_psi,
          rho,
          if (is.null(desired_width)) NA_real_ else desired_width,
          conf_level,
          if (!is.null(assurance)) assurance
        )
      )

      return(.as_dmar_tbl(Summary, conf_level = conf_level))
    }
}
