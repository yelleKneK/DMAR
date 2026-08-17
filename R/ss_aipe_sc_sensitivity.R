#' @title Sensitivity Analysis for Sample Size Planning for the Standardized ANOVA Contrast From
#' the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE)
#' Perspective for the standardized ANOVA contrast.
#'
#' @param true_psi population standardized contrast
#' @param estimated_psi estimated standardized contrast
#' @param c_weights the contrast weights
#' @param desired_width the desired full width of the obtained confidence interval
#' @param n_per_group selected sample size to use in order to determine distributional properties of at a given value of sample size
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param conf_level the desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param G number of generations (i.e., replications) of the simulation
#' @param print_iter to print the current value of the iterations
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#' @param \dots allows one to potentially include parameter values for inner functions
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across the \code{G}
#' replications. The \code{term} entries are: \code{mean_psi},
#' \code{median_psi}, \code{sd_psi} (summaries of the realized
#' standardized contrast); \code{mean_ci_width}, \code{median_ci_width},
#' \code{sd_ci_width} (summaries of the full interval widths);
#' \code{mean_ci_width_lower} and \code{mean_ci_width_upper} (mean
#' one-sided widths, measured from the observed contrast to each limit);
#' \code{pct_ci_less_w} (proportion of intervals at or below the target
#' width); \code{pct_ci_miss_low} and \code{pct_ci_miss_high}
#' (tail-specific empirical non-coverage of \code{true_psi});
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); and the input echoes \code{n_per_group},
#' \code{total_N}, \code{true_psi}, \code{estimated_psi} (NA when
#' \code{n_per_group} was supplied instead), \code{width},
#' \code{conf_level}, and \code{assurance} (present only when an
#' assurance was supplied). The proportion and Type I error rows are proportions
#' on the 0 to 1 scale, not percentages.
#'
#' @references
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators.
#' \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application,
#' and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: Accuracy in Parameter
#' Estimation via narrow confidence intervals. \emph{Psychological Methods, 11}(4), 363--385.
#'   \doi{10.1037/1082-989X.11.4.363}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#' ANCOVA and ANOVA contrasts: Sample size planning via narrow
#' confidence intervals.
#' \emph{British Journal of Mathematical and Statistical Psychology, 65},
#' 350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_sc}}, \code{\link{ss_aipe_c}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # Sensitivity analysis for a standardized three-group ANOVA contrast
#' # (-1, 0, 1) at psi = 0.5 and target full width 0.40. G is kept small
#' # here so the example runs quickly; raise it for a stable sweep.
#' set.seed(113)
#' ss_aipe_sc_sensitivity(
#'   true_psi = 0.5, estimated_psi = 0.5,
#'   c_weights = c(-1, 0, 1),
#'   desired_width = 0.40,
#'   conf_level = 0.95, G = 50, print_iter = FALSE
#' )
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sc_sensitivity <- function(true_psi = NULL, estimated_psi = NULL, c_weights, desired_width = NULL,
                                   n_per_group = NULL, assurance = NULL, conf_level = .95, G = 10000,
                                   print_iter = TRUE, save = FALSE, filename = "ss_aipe_sc_sensitivity_result.csv", ...) {
  if (is.null(estimated_psi) && is.null(n_per_group)) stop("You must specify either \'estimated_psi\' or \'n_per_group\' (i.e., the sample size per group ).", call. = FALSE)
  if (!is.null(estimated_psi) && !is.null(n_per_group)) stop("You must specify either \'estimated_psi\' or \'n_per_group\' (i.e., the per group sample size), but not both.", call. = FALSE)

  if (abs(sum(c_weights)) > 1e-8) stop("The sum of the coefficients must be zero")
  if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")

  if (!is.null(estimated_psi)) {
    n <- ss_aipe_sc(psi_standardized = estimated_psi, conf_level = conf_level, c_weights = c_weights, width = desired_width, assurance = assurance, Tolerance = 1e-7)[1, 2]
  }
  else {
    n <- n_per_group
  }

  J <- length(c_weights)

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

  # Create population group means for group 1 to group J
  signal <- TRUE
  mu_c <- rep(NA, J)
  for (cw in 1:J)
  {
    if (c_weights[cw] > 0 && signal) {
      mu_c[cw] <- psi_pop / c_weights[cw]
      signal <- FALSE
    }

    else {
      mu_c[cw] <- 0
    }
  }

  for (i in 1:G)
  {
    if (print_iter) cat(c(i), "\n")

    # generate random samples
    group_data <- array(NA, dim = c(n, J))
    x_bar <- rep(NA, J)
    sd_group <- rep(NA, J)

    for (p in 1:J)
    {
      group_data[, p] <- rnorm(n, mean = mu_c[p], sd = 1)
      x_bar[p] <- mean(group_data[, p])
      sd_group[p] <- sd(group_data[, p])
    }

    # calculate the observed psi_standardized and CI for psi_standardized
    s_pooled <- mean(sd_group)
    psi_obs[i] <- sum(c_weights * x_bar) / s_pooled
    lambda <- psi_obs[i] / sqrt(sum(c_weights^2) / n)
    lambda_limits <- conf_limits_nct(ncp = lambda, df = n * J - J, conf_level = conf_level)
    psi_limit_upper <- lambda_limits[2, 2] * sqrt(sum(c_weights^2) / n)
    psi_limit_lower <- lambda_limits[1, 2] * sqrt(sum(c_weights^2) / n)

    Full_Width[i] <- psi_limit_upper - psi_limit_lower
    Width_from_psi_obs_Lower[i] <- psi_obs[i] - psi_limit_lower
    Width_from_psi_obs_Upper[i] <- psi_limit_upper - psi_obs[i]

    Type_I_Error_Upper[i] <- psi_pop > psi_limit_upper
    Type_I_Error_Lower[i] <- psi_pop < psi_limit_lower
    Type_I_Error[i] <- Type_I_Error_Upper[i] | Type_I_Error_Lower[i]

    Lower_Limit[i] <- psi_limit_lower
    Upper_Limit[i] <- psi_limit_upper
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
      "true_psi", "estimated_psi",
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
      if (is.null(desired_width)) NA_real_ else desired_width,
      conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  return(.as_dmar_tbl(Summary, conf_level = conf_level))
}
