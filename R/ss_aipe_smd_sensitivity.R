#' Sensitivity Analysis for Sample Size Given the Accuracy in Parameter Estimation Approach for the Standardized Mean Difference
#'
#' @description
#' Performs sensitivity analysis for sample size determination for the standardized mean difference
#' given a population and a standardized mean difference. Allows one to determine the effect of being
#' wrong when estimating the population standardized mean difference in terms of the
#' width of the obtained (two-sided) confidence intervals.
#'
#' @param true_smd population standardized mean difference
#' @param supposed_smd the supposed standardized mean difference: the planning value from which
#'   \code{\link{ss_aipe_smd}} plans the per-group sample size (not an estimate computed from data);
#'   it can equal \code{true_smd} to simulate a correctly specified plan
#' @param width describe full width for the confidence interval around the population standardized mean difference
#' @param specified_n_per_group selected sample size to use in order to determine distributional properties of at a given value of sample size
#' @param assurance parameter to ensure confidence interval width with a specified degree of certainty (must be \code{NULL} or between zero and unity). When the sample size is planned from \code{supposed_smd}, an assurance above .999999998 is refused, as in \code{\link{ss_aipe_smd}}.
#' @param assurance_method passed to \code{\link{ss_aipe_smd}} with
#'   \code{assurance} when the sample size is planned from
#'   \code{supposed_smd}. Only \code{"analytic"} (the default) is offered:
#'   the plan is the smallest per-group sample size, at or above the plan
#'   without \code{assurance}, at which the exact probability that the
#'   interval is no wider than \code{width} is at least \code{assurance}.
#'   \code{"simulation"} is refused with a condition of class
#'   \code{dmar_argument_invalid}, as it is by \code{\link{ss_aipe_smd}},
#'   since that probability is exact, whether or not \code{assurance} or
#'   \code{specified_n_per_group} is given.
#' @param conf_level the desired degree of confidence (i.e., 1-Type I error rate)
#' @param G number of generations (i.e., replications) of the simulation
#' @param print_iter to print the current value of the iterations
#' @param filename an optional path for a comma separated file recording
#'   every replication (the realized standardized mean difference, the
#'   full and the two one-sided interval widths, the three non-coverage
#'   indicators, and the two confidence limits): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a header
#'   row is created otherwise, an existing file at that path is appended
#'   to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' For sensitivity analysis when planning sample size given the desire to obtain narrow confidence intervals
#' for the population standardized mean difference. Given a population value and a supposed (planning) value, one can determine
#' the effects of incorrectly specifying the population standardized mean difference (\code{true_smd}) on the
#' obtained widths of the confidence intervals. Also, one can evaluate the percent of the confidence intervals
#' that are less than the desired width (especially when modifying the \code{assurance} parameter); see \code{ss_aipe_smd})
#' Alternatively, one can specify \code{specified_n_per_group} to determine the results at a particular sample size
#' (when doing this \code{supposed_smd} cannot be specified).
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
#' \code{pct_ci_less_width} (proportion of intervals at or below the target
#' width); \code{pct_ci_miss_low} and \code{pct_ci_miss_high}
#' (tail-specific empirical non-coverage of \code{true_smd});
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); and the input echoes \code{n_per_group} (the size simulated at, planned or specified),
#' \code{total_N}, \code{true_smd}, \code{supposed_smd} (NA when
#' \code{specified_n_per_group} was supplied instead), \code{width},
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
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_smd}}
#'
#' @examples
#' # Since 'true_smd' equals 'supposed_smd', this usage
#' # returns the results of a correctly specified situation.
#' # Note that 'G' should be large (50 is used to make the example run easily)
#' set.seed(113)
#' Res.1 <- ss_aipe_smd_sensitivity(true_smd=.5, supposed_smd=.5, width=.30,
#'                                  assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)
#'
#' # Objects contained in the 'summary'.
#' Res.1$term
#'
#' # True standardized mean difference is .4, but specified at .5.
#' # Change 'G' to some large number (e.g., G=5,000)
#' Res.2 <- ss_aipe_smd_sensitivity(true_smd=.4, supposed_smd=.5, width=.30,
#'                                  assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)
#'
#' # The effect of the misspecification on mean confidence intervals is:
#' Res.2[1,]
#'
#' # True standardized mean difference is .5, but specified at .4.
#' Res.3 <- ss_aipe_smd_sensitivity(true_smd=.5, supposed_smd=.4, width=.30,
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
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export


ss_aipe_smd_sensitivity <- function(true_smd = NULL, supposed_smd = NULL, width = NULL, specified_n_per_group = NULL,
                                    assurance = NULL, assurance_method = "analytic", conf_level = .95, G = 1000, print_iter = FALSE,
                                    filename = NULL,
                                    seed = NULL) {
  .dmar_check_assurance_single(assurance)
  # The planner's method, refused here in this function's name when it is
  # not one ss_aipe_smd() offers (R/ss_aipe_assurance_internals.R).
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_smd_sensitivity", simulation = FALSE,
    both = FALSE)
  .dmar_check_numbers("ss_aipe_smd_sensitivity", "true_smd")
  if (is.null(supposed_smd) && is.null(specified_n_per_group)) .dmar_refuse_argument(
    "You must specify either \'supposed_smd\' or \'specified_n_per_group\' (i.e., the per group sample size).",
    arg = c("supposed_smd", "specified_n_per_group"))
  if (!is.null(supposed_smd) && !is.null(specified_n_per_group)) .dmar_refuse_argument(
    "You must specify either \'supposed_smd\' or \'specified_n_per_group\' (i.e., the per group sample size), but not both.",
    arg = c("supposed_smd", "specified_n_per_group"))
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  if (!is.null(supposed_smd)) {
    n <- ss_aipe_smd(
      true_smd = supposed_smd, conf_level = conf_level, width = width,
      assurance = assurance, assurance_method = assurance_method
    )[1, 2]
  }
  else {
    n <- specified_n_per_group
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

  delta <- true_smd
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
    type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, ci_lower = Low_Limit,
    ci_upper = Upper_Limit
  )

  if (!is.null(filename)) {
    if (file.exists(filename) && file.size(filename) > 0) {
      message("The file '", filename, "' already exists; the simulation results are appended to it.")
      utils::write.table(result, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    } else {
      utils::write.table(result, filename, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  out <- data.frame(
    term = c(
      "mean_smd", "median_smd", "sd_smd",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "mean_ci_width_lower", "mean_ci_width_upper",
      "pct_ci_less_width",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "n_per_group", "total_N",
      "true_smd", "supposed_smd",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(d), median(d), sd(d),
      mean(Full_Width), median(Full_Width), sd(Full_Width),
      mean(Width_from_d_Lower), mean(Width_from_d_Upper),
      mean(Full_Width <= width),
      mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
      n, 2 * n,
      true_smd,
      if (is.null(supposed_smd)) NA_real_ else supposed_smd,
      if (is.null(width)) NA_real_ else width,
      conf_level,
      if (!is.null(assurance)) assurance
    )
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(Full_Width), conf_level)
}
