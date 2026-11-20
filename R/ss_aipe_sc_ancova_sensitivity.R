#' Sensitivity Analysis for the Sample Size Planning Method for Standardized ANCOVA Contrast
#'
#' @description
#' Sensitivity analysis for the sample size planning method with the goal to obtain sufficiently
#' narrow confidence intervals for standardized ANCOVA complex contrasts.
#'
#' @param true_std_contrast the population standardized ANCOVA contrast
#' @param supposed_std_contrast the supposed standardized ANCOVA contrast, the planning value given to \code{\link{ss_aipe_sc_ancova}} to choose the sample size
#' @param c_weights the contrast weights
#' @param width the desired full width of the obtained confidence interval.
#'   Required when \code{supposed_std_contrast} is given, because the
#'   sample size is planned from it; a call without it is refused with a
#'   condition of class \code{dmar_argument_invalid}.
#' @param specified_n_per_group selected sample size per group to use in order to determine distributional properties of a given value of sample size
#' @param true_mean_X the population mean for the covariate
#' @param true_sd_X the population standard deviation of the covariate
#' @param true_cor the population correlation coefficient between the response and the covariate
#' @param divisor which error standard deviation to be used in standardizing the contrast; the value can be either \code{"s_ancova"} or \code{"s_anova"}
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty. When the sample size is planned from \code{supposed_std_contrast}, it is passed to \code{\link{ss_aipe_sc_ancova}} and must be \code{NULL} or greater than .5 and less than 1; with \code{specified_n_per_group} it is only reported in the result.
#' @param assurance_method Passed to \code{\link{ss_aipe_sc_ancova}},
#'   which plans the sample size from \code{supposed_std_contrast} with
#'   \code{assurance}: \code{"analytic"} (the default) for the published
#'   plan of Lai and Kelley (2012), which does not compute the probability
#'   it attains and can attain less than \code{assurance}, and
#'   \code{"simulation"} for the smallest sample size per group at which
#'   the probability of an interval no wider than \code{width}, estimated
#'   by the planner from its default \code{G} of 10000 simulated samples
#'   with the random number stream that \code{seed} sets, is at least
#'   \code{assurance}. Not used with \code{specified_n_per_group} or
#'   without \code{assurance}.
#' @param conf_level the desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param G number of generations (i.e., replications) of the simulation
#' @param print_iter to print the current value of the iterations
#' @param filename Optional path of a CSV file to receive the per-replication
#'   results (the observed standardized contrast, \code{std_contrast}, the
#'   full and one-sided
#'   interval widths, the tail and overall misses, and the confidence
#'   limits), appended when the file already exists and created otherwise;
#'   the default \code{NULL} writes nothing, and a throwaway run that wants
#'   the file should point it at \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
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
#' replications. The \code{term} entries are: \code{mean_std_contrast},
#' \code{median_std_contrast}, \code{sd_std_contrast} (summaries of the realized
#' standardized ANCOVA contrast); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the full
#' interval widths); \code{mean_ci_width_lower} and
#' \code{mean_ci_width_upper} (mean one-sided widths, measured from the
#' observed contrast to each limit); \code{pct_ci_less_width} (proportion of
#' intervals at or below the target width); \code{pct_ci_miss_low} and
#' \code{pct_ci_miss_high} (tail-specific empirical non-coverage of
#' \code{true_std_contrast}); \code{total_type_I_error} (overall empirical
#' non-coverage, the sum of the two tails); and the input echoes
#' \code{n_per_group} (the size simulated at, planned or specified), \code{total_N}, \code{true_std_contrast},
#' \code{supposed_std_contrast} (NA when \code{specified_n_per_group} was supplied
#' instead), \code{true_cor}, \code{width}, \code{conf_level}, and
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
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
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
#' # of 0.4, and a planning target width of 0.5. G = 50 keeps the
#' # example quick; a reported sensitivity analysis deserves the default of
#' # G = 10000 replications.
#' set.seed(113)
#' ss_aipe_sc_ancova_sensitivity(
#'   true_std_contrast = 0.5, supposed_std_contrast = 0.5,
#'   c_weights = c(-1, 0, 1),
#'   width = 0.5, true_cor = 0.4,
#'   conf_level = 0.95, G = 50, print_iter = FALSE
#' )
#'
#' @keywords design
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export


ss_aipe_sc_ancova_sensitivity <- function(true_std_contrast = NULL, supposed_std_contrast = NULL, c_weights,
                                          width = NULL, specified_n_per_group = NULL, true_mean_X = 0, true_sd_X = 1, true_cor, divisor = "s_ancova",
                                          assurance = NULL, conf_level = .95, G = 10000, print_iter = TRUE,
                                          filename = NULL,
                                          assurance_method = c("analytic", "simulation"),
                                          seed = NULL) {
  .dmar_check_conf_level(conf_level)
  .dmar_check_assurance_single(assurance)
  .dmar_check_numbers("ss_aipe_sc_ancova_sensitivity", c("true_std_contrast", "supposed_std_contrast"))
  if (divisor != "s_ancova" && divisor != "s_anova") .dmar_refuse_argument(
    "The argument 'divisor' must be either 's_ancova' or 's_anova'",
    arg = "divisor")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  # A sample size is planned from 'width' when 'supposed_std_contrast'
  # is given; without the check a missing width met R's "argument is of
  # length zero" inside the planner.
  if (!is.null(supposed_std_contrast)) {
    .dmar_check_width(width, "the standardized contrast")
  }
  # The planner's method (R/ss_aipe_assurance_internals.R), checked here in
  # this function's name.
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_sc_ancova_sensitivity", simulation = TRUE)
  .dmar_local_seed(seed)
  # The plan from the supposed value (with "simulation", the planner's G,
  # 10000, and the seed set above).
  plan_n <- function(...) {
    ss_aipe_sc_ancova(..., width = width, c_weights = c_weights,
                      conf_level = conf_level, assurance = assurance,
                      assurance_method = assurance_method)[1, 2]
  }

  if (divisor == "s_ancova") {
    ##  standardized ANCOVA contrast using s_ancova as divisor

    if (is.null(supposed_std_contrast) && is.null(specified_n_per_group)) .dmar_refuse_argument(
      "You must specify either \'supposed_std_contrast\' or \'specified_n_per_group\' (i.e., the sample size per group).",
      arg = c("supposed_std_contrast", "specified_n_per_group"))
    if (!is.null(supposed_std_contrast) && !is.null(specified_n_per_group)) .dmar_refuse_argument(
      "You must specify either \'supposed_std_contrast\' or \'specified_n_per_group\' (i.e., the sample size per group), but not both.",
      arg = c("supposed_std_contrast", "specified_n_per_group"))

    if (abs(sum(c_weights)) > 1e-8) .dmar_refuse_argument(
      "'c_weights' must sum to zero.", arg = "c_weights")
    # All-zero weights define no contrast. With 'specified_n_per_group'
    # no planner is called, and they stopped below with "object 'k' not
    # found".
    .dmar_check_contrast_nonzero(c_weights)
    if (sum(c_weights[c_weights > 0]) > 1) .dmar_refuse_argument(
      "The positive values of 'c_weights' must not sum to more than 1; give the weights as fractions, such as c(1/2, 1/2, -1).", arg = "c_weights")

    J <- length(c_weights)
    sigma_ancova <- 2
    sigma_anova <- sqrt(sigma_ancova^2 / (1 - true_cor^2))
    mu_y <- rep(NA, J)

    if (!is.null(supposed_std_contrast)) {
      n <- plan_n(true_std_contrast = supposed_std_contrast)
    }
    else {
      n <- specified_n_per_group
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
    psi_pop <- true_std_contrast

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
      random_data <- simulate_ancova_data(true_mean_Y = mu_y, true_mean_X = true_mean_X,
                                          true_sd_Y = sigma_anova, true_sd_X = true_sd_X,
                                          true_cor = true_cor, a = J, n_per_group = n)
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

      # calculate the observed standardized contrast and its CI
      psi_obs[g] <- sum(c_weights * y_bar_adj) / s_ancova

      ci_psi <- ci_sc_ancova(
        adj_means = y_bar_adj, s_ancova = s_ancova, c_weights = c_weights, n_per_group = n, cov_means = cov_means,
        SSwithin_x = SSwithin_x, conf_level = conf_level
      )

      psi_limit_lower <- ci_psi$value[ci_psi$term == "ci_lower"]
      psi_limit_upper <- ci_psi$value[ci_psi$term == "ci_upper"]

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
      std_contrast = psi_obs, full_width = Full_Width, width_lower = Width_from_psi_obs_Lower,
      width_upper = Width_from_psi_obs_Upper, type_I_error_upper = Type_I_Error_Upper,
      type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, ci_lower = Lower_Limit,
      ci_upper = Upper_Limit
    )

    if (!is.null(filename)) {
      suppressWarnings(file_exist <- try(utils::read.csv(filename), silent = TRUE))
      if (!is.null(dim(file_exist))) {
        utils::write.table(Results, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        cat("The file '", filename, "' already exists; the simulation results were appended to it.\n", sep = "")
      } else {
        utils::write.table(Results, filename, sep = ",", row.names = FALSE, append = FALSE)
      }
    }

    Summary <- data.frame(
      term = c(
        "mean_std_contrast", "median_std_contrast", "sd_std_contrast",
        "mean_ci_width", "median_ci_width", "sd_ci_width",
        "mean_ci_width_lower", "mean_ci_width_upper",
        "pct_ci_less_width",
        "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
        "n_per_group", "total_N",
        "true_std_contrast", "supposed_std_contrast", "true_cor",
        "width", "conf_level",
        if (!is.null(assurance)) "assurance"
      ),
      value = c(
        mean(psi_obs), median(psi_obs), sd(psi_obs),
        mean(Full_Width), median(Full_Width), sd(Full_Width),
        mean(Width_from_psi_obs_Lower), mean(Width_from_psi_obs_Upper),
        mean(Full_Width <= width),
        mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
        n, n * J,
        true_std_contrast,
        if (is.null(supposed_std_contrast)) NA_real_ else supposed_std_contrast,
        true_cor,
        if (is.null(width)) NA_real_ else width,
        conf_level,
        if (!is.null(assurance)) assurance
      )
    )

    return(.ss_sens_record(.as_dmar_tbl(Summary, conf_level = conf_level), G,
                         !is.na(Full_Width), conf_level))
  }
  #################################################################################################
  #################################################################################################

  if (divisor == "s_anova")
  ## use s_anova to standardize the ANCOVA contrast

    {
      if (is.null(supposed_std_contrast) && is.null(specified_n_per_group)) .dmar_refuse_argument(
        "You must specify either \'supposed_std_contrast\' or \'specified_n_per_group\' (i.e., the sample size per group ).",
        arg = c("supposed_std_contrast", "specified_n_per_group"))
      if (!is.null(supposed_std_contrast) && !is.null(specified_n_per_group)) .dmar_refuse_argument(
        "You must specify either \'supposed_std_contrast\' or \'specified_n_per_group\' (i.e., the sample size per group), but not both.",
        arg = c("supposed_std_contrast", "specified_n_per_group"))

      if (abs(sum(c_weights)) > 1e-8) .dmar_refuse_argument(
        "'c_weights' must sum to zero.", arg = "c_weights")
      .dmar_check_contrast_nonzero(c_weights)
      if (sum(c_weights[c_weights > 0]) > 1) .dmar_refuse_argument(
        "The positive values of 'c_weights' must not sum to more than 1; give the weights as fractions, such as c(1/2, 1/2, -1).", arg = "c_weights")

      sigma_anova <- 2
      sigma_ancova <- sqrt(sigma_anova^2 * (1 - true_cor^2))

      if (!is.null(supposed_std_contrast)) {
        n <- plan_n(true_contrast = supposed_std_contrast * sigma_anova,
                    true_sd_ancova = sigma_ancova, true_sd_anova = sigma_anova,
                    divisor = "s_anova")
      }
      else {
        n <- specified_n_per_group
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
      psi_pop <- true_std_contrast

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
        random_data <- simulate_ancova_data(true_mean_Y = mu_y, true_mean_X = true_mean_X,
                                            true_sd_Y = sigma_anova, true_sd_X = true_sd_X,
                                            true_cor = true_cor, a = J, n_per_group = n)
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

        # calculate the observed standardized contrast and its CI based on s_anova
        psi_obs[g] <- sum(c_weights * y_bar_adj) / s_anova

        ci_psi <- ci_sc_ancova(
          contrast = sum(c_weights * y_bar_adj), s_anova = s_anova, s_ancova = s_ancova, standardizer = "s_anova",
          c_weights = c_weights, n_per_group = n, cov_means = cov_means, SSwithin_x = SSwithin_x, conf_level = conf_level
        )

        psi_limit_lower <- ci_psi$value[ci_psi$term == "ci_lower"]
        psi_limit_upper <- ci_psi$value[ci_psi$term == "ci_upper"]

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
        std_contrast = psi_obs, full_width = Full_Width, width_lower = Width_from_psi_obs_Lower,
        width_upper = Width_from_psi_obs_Upper, type_I_error_upper = Type_I_Error_Upper,
        type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, ci_lower = Lower_Limit,
        ci_upper = Upper_Limit
      )

      if (!is.null(filename)) {
        suppressWarnings(file_exist <- try(utils::read.csv(filename), silent = TRUE))
        if (!is.null(dim(file_exist))) {
          utils::write.table(Results, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          cat("The file '", filename, "' already exists; the simulation results were appended to it.\n", sep = "")
        } else {
          utils::write.table(Results, filename, sep = ",", row.names = FALSE, append = FALSE)
        }
      }

      Summary <- data.frame(
        term = c(
          "mean_std_contrast", "median_std_contrast", "sd_std_contrast",
          "mean_ci_width", "median_ci_width", "sd_ci_width",
          "mean_ci_width_lower", "mean_ci_width_upper",
          "pct_ci_less_width",
          "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
          "n_per_group", "total_N",
          "true_std_contrast", "supposed_std_contrast", "true_cor",
          "width", "conf_level",
          if (!is.null(assurance)) "assurance"
        ),
        value = c(
          mean(psi_obs), median(psi_obs), sd(psi_obs),
          mean(Full_Width), median(Full_Width), sd(Full_Width),
          mean(Width_from_psi_obs_Lower), mean(Width_from_psi_obs_Upper),
          mean(Full_Width <= width),
          mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
          n, n * J,
          true_std_contrast,
          if (is.null(supposed_std_contrast)) NA_real_ else supposed_std_contrast,
          true_cor,
          if (is.null(width)) NA_real_ else width,
          conf_level,
          if (!is.null(assurance)) assurance
        )
      )

      return(.ss_sens_record(.as_dmar_tbl(Summary, conf_level = conf_level), G,
                         !is.na(Full_Width), conf_level))
    }
}
