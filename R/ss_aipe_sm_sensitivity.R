#' Sensitivity Analysis for Sample Size Planning for the Standardized Mean From the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description
#' Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE)
#' Perspective for the standardized mean.
#'
#' @param true_sm population standardized mean
#' @param estimated_sm estimated standardized mean
#' @param desired_width desired full width of the confidence interval for the population standardized mean
#' @param specified_N selected sample size to use in order to determine distributional properties of a given value of sample size
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity)
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
#' replications. The \code{term} entries are: \code{mean_sm},
#' \code{median_sm}, \code{sd_sm} (summaries of the realized
#' standardized mean); \code{mean_ci_width}, \code{median_ci_width},
#' \code{sd_ci_width} (summaries of the full interval widths);
#' \code{mean_ci_width_lower} and \code{mean_ci_width_upper} (mean
#' one-sided widths, measured from the observed standardized mean to
#' each limit); \code{pct_ci_less_w} (proportion of intervals at or
#' below the target width); \code{pct_ci_miss_low} and
#' \code{pct_ci_miss_high} (tail-specific empirical non-coverage of
#' \code{true_sm}); \code{total_type_I_error} (overall empirical
#' non-coverage, the sum of the two tails); and the input echoes
#' \code{total_N}, \code{true_sm}, \code{estimated_sm} (NA when
#' \code{specified_N} was supplied instead), \code{width},
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
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the
#' standardized mean difference: Bootstrap and parametric confidence intervals,
#' \emph{Educational and Psychological Measurement, 65}, 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application,
#' and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_sm}}
#'
#' @examples
#' # Since 'true_sm' equals 'estimated_sm', this usage
#' # returns the results of a correctly specified situation.
#' # Note that 'G' should be large (10 is used to make the
#' # example run easily)
#' #Res.1 <- ss_aipe_sm_sensitivity(true_sm=10, estimated_sm=10,
#' #desired_width=.5, assurance=.95, conf_level=.95, G=10,
#' #print_iter=FALSE)
#'
#' # Objects contained in the 'Summary'.
#' # Res.1$term
#'
#' # What proportion of the obtained full widths are narrower than the
#' # desired one?
#' # Res.1[which(Res.1$term == 'pct_ci_less_w'),2]
#'
#' # True standardized mean difference is 10, but specified at 12.
#' # Change 'G' to some large number (e.g., G=20)
#' #Res.2 <- ss_aipe_sm_sensitivity(true_sm=10, estimated_sm=12,
#' #desired_width=.5, assurance=NULL, conf_level=.95, G=20)
#'
#' # The effect of the misspecification on mean confidence intervals is:
#' # Res.2[which(Res.2$term == 'mean_ci_width'),2]
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_sm_sensitivity <- function(true_sm = NULL, estimated_sm = NULL, desired_width = NULL, specified_N = NULL,
                                   assurance = NULL, conf_level = .95, G = 10000, print_iter = TRUE,
                                   save = FALSE, filename = "ss_aipe_sm_sensitivity_result.csv", ...) {
  if (is.null(estimated_sm) && is.null(specified_N)) stop("You must specify either 'estimated_sm' or 'specified_N' (i.e., the total sample size ).", call. = FALSE)
  if (!is.null(estimated_sm) && !is.null(specified_N)) stop("You must specify either 'estimated_sm' or 'specified_N' (i.e., the total sample size), but not both.", call. = FALSE)

  # The planner search and the Monte Carlo loop both call conf_limits_nct() many
  # times. When the observed standardized mean is large the noncentrality
  # parameter can exceed the magnitude at which R's pt()/qt() stay accurate,
  # which conf_limits_nct() signals via a warning; surfacing it on every
  # evaluation produces dozens of identical messages. We muffle only that
  # specific warning with withCallingHandlers, count it, and emit a single
  # summary warning at the end (via on.exit so it fires on any return path).
  # Warnings from any other source still surface normally.
  .nct_ncp_count <- 0L
  on.exit({
    if (.nct_ncp_count > 0L) {
      warning(sprintf(
        "During the sample size search and Monte Carlo sensitivity loop, the noncentrality parameter exceeded the accurate range of R's noncentral t functions in %d evaluations (see ?conf_limits_nct). Those evaluations may be inaccurate.",
        .nct_ncp_count
      ), call. = FALSE)
    }
  }, add = TRUE)

  .muffle_nct_ncp <- function(w) {
    if (grepl("noncentrality parameter exceeds", conditionMessage(w), fixed = TRUE)) {
      .nct_ncp_count <<- .nct_ncp_count + 1L
      invokeRestart("muffleWarning")
    }
  }

  if (!is.null(estimated_sm)) {
    n <- withCallingHandlers(
      ss_aipe_sm(sm = estimated_sm, conf_level = conf_level, width = desired_width, assurance = assurance, Tolerance = 1e-7)[1, 2],
      warning = .muffle_nct_ncp
    )
  }
  else {
    n <- specified_N
  }

  G <- G
  Full_Width <- rep(NA, G)
  Width_from_sm_obs_Lower <- rep(NA, G)
  Width_from_sm_obs_Upper <- rep(NA, G)
  Type_I_Error_Upper <- rep(NA, G)
  Type_I_Error_Lower <- rep(NA, G)
  Type_I_Error <- rep(NA, G)
  Lower_Limit <- rep(NA, G)
  Upper_Limit <- rep(NA, G)
  sm_obs <- rep(NA, G)
  sm_pop <- true_sm

  withCallingHandlers(
    {
      for (i in 1:G)
      {
        if (print_iter) cat(c(i), "\n")
        sample_data <- rnorm(n, mean = sm_pop, sd = 1)
        x_bar <- mean(sample_data)
        sm_obs[i] <- x_bar / sd(sample_data)

        lambda <- sm_obs[i] * sqrt(n)
        lambda_limits <- conf_limits_nct(ncp = lambda, df = n - 1, conf_level = conf_level)
        sm_limit_upper <- lambda_limits[2, 2] / sqrt(n)
        sm_limit_lower <- lambda_limits[1, 2] / sqrt(n)

        Full_Width[i] <- sm_limit_upper - sm_limit_lower
        Width_from_sm_obs_Lower[i] <- sm_obs[i] - sm_limit_lower
        Width_from_sm_obs_Upper[i] <- sm_limit_upper - sm_obs[i]

        Type_I_Error_Upper[i] <- sm_pop > sm_limit_upper
        Type_I_Error_Lower[i] <- sm_pop < sm_limit_lower
        Type_I_Error[i] <- Type_I_Error_Upper[i] | Type_I_Error_Lower[i]

        Lower_Limit[i] <- sm_limit_lower
        Upper_Limit[i] <- sm_limit_upper
      }
    },
    warning = .muffle_nct_ncp
  )

  Results <- data.frame(
    sm_obs = sm_obs, full_width = Full_Width, width_lower = Width_from_sm_obs_Lower,
    width_upper = Width_from_sm_obs_Upper, type_I_error_upper = Type_I_Error_Upper,
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
      "mean_sm", "median_sm", "sd_sm",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "mean_ci_width_lower", "mean_ci_width_upper",
      "pct_ci_less_w",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "total_N",
      "true_sm", "estimated_sm",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(sm_obs), median(sm_obs), sd(sm_obs),
      mean(Full_Width), median(Full_Width), sd(Full_Width),
      mean(Width_from_sm_obs_Lower), mean(Width_from_sm_obs_Upper),
      mean(Full_Width <= desired_width),
      mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
      n,
      true_sm,
      if (is.null(estimated_sm)) NA_real_ else estimated_sm,
      if (is.null(desired_width)) NA_real_ else desired_width,
      conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  return(.as_dmar_tbl(Summary, conf_level = conf_level))
}
