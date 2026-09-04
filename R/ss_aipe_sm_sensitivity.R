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
#' @param filename Optional path to a CSV file; when supplied, the per-replication results (the
#'   observed standardized mean, the full and one-sided widths, the tail misses, and the interval
#'   limits) are written there, appended when the file already exists, and a throwaway run should
#'   point it at \code{tempfile(fileext = ".csv")}; the default \code{NULL} writes nothing.
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
#' # With true_sm equal to estimated_sm the planning value is correct, and
#' # the sweep reports what a correctly specified plan delivers: the share of
#' # intervals no wider than the target should sit near the assurance. G = 10
#' # keeps the example quick; a reported sensitivity study deserves the
#' # default G = 10000.
#' set.seed(113)
#' res_correct <- ss_aipe_sm_sensitivity(true_sm = 1.5, estimated_sm = 1.5,
#'   desired_width = 0.5, assurance = 0.95, conf_level = 0.95, G = 10,
#'   print_iter = FALSE)
#' res_correct
#'
#' # The terms the summary reports.
#' res_correct$term
#'
#' # The proportion of realized full widths no wider than the target.
#' res_correct$value[res_correct$term == "pct_ci_less_w"]
#'
#' # The population standardized mean is 1.5 but the plan assumed 2, so the
#' # planner sizes the study for a wider sampling distribution than the data
#' # will show, and the realized intervals come in narrower than the target.
#' set.seed(113)
#' res_misspecified <- ss_aipe_sm_sensitivity(true_sm = 1.5, estimated_sm = 2,
#'   desired_width = 0.5, G = 20, print_iter = FALSE)
#'
#' # The effect of the misspecification on the mean interval width.
#' res_misspecified$value[res_misspecified$term == "mean_ci_width"]
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_sm_sensitivity <- function(true_sm = NULL, estimated_sm = NULL, desired_width = NULL, specified_N = NULL,
                                   assurance = NULL, conf_level = .95, G = 10000, print_iter = TRUE,
                                   filename = NULL) {
  if (is.null(estimated_sm) && is.null(specified_N)) stop("You must specify either 'estimated_sm' or 'specified_N' (i.e., the total sample size ).", call. = FALSE)
  if (!is.null(estimated_sm) && !is.null(specified_N)) stop("You must specify either 'estimated_sm' or 'specified_N' (i.e., the total sample size), but not both.", call. = FALSE)
  .check_filename(filename)

  # The planner search and the Monte Carlo loop both call ci_nc_t() many
  # times. When the observed standardized mean is large the noncentrality
  # parameter can exceed the magnitude at which R's pt()/qt() stay accurate,
  # which ci_nc_t() signals via a warning; surfacing it on every
  # evaluation produces dozens of identical messages. We muffle only that
  # specific warning with withCallingHandlers, count it, and emit a single
  # summary warning at the end (via on.exit so it fires on any return path).
  # Warnings from any other source still surface normally. The count lives in
  # an environment made for it so the handler updates it by ordinary assignment.
  state <- new.env(parent = emptyenv())
  state$nct_ncp_count <- 0L
  on.exit({
    if (state$nct_ncp_count > 0L) {
      warning(sprintf(
        "During the sample size search and Monte Carlo sensitivity loop, the noncentrality parameter exceeded the accurate range of R's noncentral t functions in %d evaluations (see ?ci_nc_t). Those evaluations may be inaccurate.",
        state$nct_ncp_count
      ), call. = FALSE)
    }
  }, add = TRUE)

  .muffle_nct_ncp <- function(w) {
    if (grepl("noncentrality parameter exceeds", conditionMessage(w), fixed = TRUE)) {
      state$nct_ncp_count <- state$nct_ncp_count + 1L
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
        lambda_limits <- ci_nc_t(ncp = lambda, df = n - 1, conf_level = conf_level)
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

  # The family's writer: a new file gets a header, an existing file is
  # appended to without one.
  if (!is.null(filename)) .write_sensitivity_csv(Results, filename)

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
