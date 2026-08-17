#' Sensitivity Analysis for Sample Size Given the Accuracy in Parameter Estimation Approach for the Standardized Mean Difference
#'
#' @description
#' Performs sensitivity analysis for sample size determination for the standardized mean difference
#' given a population and a standardized mean difference. Allows one to determine the effect of being
#' wrong when estimating the population standardized mean difference in terms of the
#' width of the obtained (two-sided) confidence intervals.
#'
#' @param true_delta population standardized mean difference
#' @param estimated_delta estimated standardized mean difference; can be \code{true_delta} to perform standard simulations
#' @param desired_width describe full width for the confidence interval around the population standardized mean difference
#' @param n_per_group selected sample size to use in order to determine distributional properties of at a given value of sample size
#' @param assurance parameter to ensure confidence interval width with a specified degree of certainty (must be \code{NULL} or between zero and unity)
#' @param conf_level the desired degree of confidence (i.e., 1-Type I error rate)
#' @param G number of generations (i.e., replications) of the simulation
#' @param print_iter to print the current value of the iterations
#' @param save option to save simulation results. It can be saved with \code{save = TRUE} outside of the printed results
#' @param filename the name of the file that simulation results will be saved to
#' @param ... for modifying parameters of functions this function calls
#'
#' @details
#' For sensitivity analysis when planning sample size given the desire to obtain narrow confidence intervals
#' for the population standardized mean difference. Given a population value and an estimated value, one can determine
#' the effects of incorrectly specifying the population standardized mean difference (\code{true_delta}) on the
#' obtained widths of the confidence intervals. Also, one can evaluate the percent of the confidence intervals
#' that are less than the desired width (especially when modifying the \code{assurance} parameter); see \code{ss_aipe_smd})
#' Alternatively, one can specify \code{n_per_group} to determine the results at a particular sample size
#' (when doing this \code{estimated_delta} cannot be specified).
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across the \code{G}
#' replications. The \code{term} entries are: \code{mean_smd},
#' \code{median_smd}, \code{sd_smd} (summaries of the realized
#' standardized mean difference); \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width} (summaries of the full
#' interval widths); \code{mean_ci_width_lower} and
#' \code{mean_ci_width_upper} (mean one-sided widths, measured from the
#' observed standardized mean difference to each limit);
#' \code{pct_ci_less_w} (proportion of intervals at or below the target
#' width); \code{pct_ci_miss_low} and \code{pct_ci_miss_high}
#' (tail-specific empirical non-coverage of \code{true_delta});
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); and the input echoes \code{n_per_group},
#' \code{total_N}, \code{true_delta}, \code{estimated_delta} (NA when
#' \code{n_per_group} was supplied instead), \code{width},
#' \code{conf_level}, and \code{assurance} (present only when an
#' assurance was supplied). The proportion and Type I error rows are
#' proportions on the 0 to 1 scale, not percentages.
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
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
#' difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_smd}}
#'
#' @examples
#' # Since 'true_delta' equals 'estimated_delta', this usage
#' # returns the results of a correctly specified situation.
#' # Note that 'G' should be large (50 is used to make the example run easily)
#' set.seed(113)
#' Res.1 <- ss_aipe_smd_sensitivity(true_delta=.5, estimated_delta=.5, desired_width=.30,
#'                                  assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)
#'
#' # Objects contained in the 'summary'.
#' Res.1$term
#'
#' # True standardized mean difference is .4, but specified at .5.
#' # Change 'G' to some large number (e.g., G=5,000)
#' Res.2 <- ss_aipe_smd_sensitivity(true_delta=.4, estimated_delta=.5, desired_width=.30,
#'                                  assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)
#'
#' # The effect of the misspecification on mean confidence intervals is:
#' Res.2[1,]
#'
#' # True standardized mean difference is .5, but specified at .4.
#' Res.3 <- ss_aipe_smd_sensitivity(true_delta=.5, estimated_delta=.4, desired_width=.30,
#'                                  assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)
#'
#' # The effect of the misspecification on mean confidence intervals is:
#' Res.3[1,]
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_smd_sensitivity <- function(true_delta = NULL, estimated_delta = NULL, desired_width = NULL, n_per_group = NULL,
                                    assurance = NULL, conf_level = .95, G = 1000, print_iter = FALSE,
                                    save = FALSE, filename = "ss_aipe_smd_sensitivity_result.csv", ...) {
  if (is.null(estimated_delta) && is.null(n_per_group)) stop("You must specify either \'estimated_delta\' or \'n_per_group\' (i.e., the per group sample size).", call. = FALSE)
  if (!is.null(estimated_delta) && !is.null(n_per_group)) stop("You must specify either \'estimated_delta\' or \'n_per_group\' (i.e., the per group sample size), but not both.", call. = FALSE)

  if (!is.null(estimated_delta)) {
    n <- ss_aipe_smd(
      delta = estimated_delta, conf_level = conf_level, width = desired_width,
      assurance = assurance
    )[1, 2]
  }
  else {
    n <- n_per_group
  }

  G <- G

  Full_Width <- rep(NA, G)
  Width_from_d_Lower <- rep(NA, G)
  Width_from_d_Upper <- rep(NA, G)
  Type_I_Error_Upper <- rep(NA, G)
  Type_I_Error_Lower <- rep(NA, G)
  Type_I_Error <- rep(NA, G)
  Low_Limit <- rep(NA, G)
  Upper_Limit <- rep(NA, G)
  d <- rep(NA, G)

  delta <- true_delta
  for (i in 1:G)
  {
    if (print_iter == TRUE) cat(c(i), "\n")
    d[i] <- smd(group_1 = rnorm(n, delta, 1), group_2 = rnorm(n, 0, 1))[1, 2]

    ci_delta <- ci_smd(smd = d[i], n_1 = n, n_2 = n, conf_level = conf_level)

    Full_Width[i] <- ci_delta[3, 2] - ci_delta[1, 2]

    Width_from_d_Lower[i] <- d[i] - ci_delta[1, 2]
    Width_from_d_Upper[i] <- ci_delta[3, 2] - d[i]

    Type_I_Error_Upper[i] <- delta > ci_delta[3, 2]
    Type_I_Error_Lower[i] <- delta < ci_delta[1, 2]
    Type_I_Error[i] <- Type_I_Error_Upper[i] | Type_I_Error_Lower[i]

    Low_Limit[i] <- ci_delta[1, 2]
    Upper_Limit[i] <- ci_delta[3, 2]
  }

  result <- data.frame(
    d = d, full_width = Full_Width, width_lower = Width_from_d_Lower,
    width_upper = Width_from_d_Upper, type_I_error_upper = Type_I_Error_Upper,
    type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, lower_limit = Low_Limit,
    upper_limit = Upper_Limit
  )

  if (save) {
    result_file <- filename
    # print("Simulation results will be saved to a .csv file")
    suppressWarnings(file_exist <- try(utils::read.csv(result_file), silent = TRUE))
    if (!is.null(dim(file_exist))) {
      utils::write.table(result, result_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      cat("A file in the local directory has the same name as the file where simulation", "\n", "results will be saved to. Simulation results will be appended to this file.", "\n", sep = "")
    } else {
      utils::write.table(result, result_file, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  out <- data.frame(
    term = c(
      "mean_smd", "median_smd", "sd_smd",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "mean_ci_width_lower", "mean_ci_width_upper",
      "pct_ci_less_w",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "n_per_group", "total_N",
      "true_delta", "estimated_delta",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(d), median(d), sd(d),
      mean(Full_Width), median(Full_Width), sd(Full_Width),
      mean(Width_from_d_Lower), mean(Width_from_d_Upper),
      mean(Full_Width <= desired_width),
      mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
      n, 2 * n,
      true_delta,
      if (is.null(estimated_delta)) NA_real_ else estimated_delta,
      if (is.null(desired_width)) NA_real_ else desired_width,
      conf_level,
      if (!is.null(assurance)) assurance
    )
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
